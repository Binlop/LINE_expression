import requests
from bs4 import BeautifulSoup
import schedule
import time
from urllib.parse import urlencode
import pandas as pd
from LINE_preparation import LINE_Processing
from exons_coverage import SamtoolsCoverage
import subprocess

class BlastProtein:

    NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    def __init__(self, seq: str):
        self.seq = seq
        self.blast_data = []

    def blast(self):
        params = {
            "DATABASE": "nr",
            "ENTREZ_QUERY": "(none)",
            "EXPECT": "10.0",
            "HITLIST_SIZE": "50",
            "PROGRAM": "blastx",
            "QUERY": self.seq,
            "CMD": "Put"
        }
            
        rid = None
        i = [0] 
        message = urlencode(params).encode()

        headers = {"User-Agent": "BiopythonClient"}
        response = requests.put(self.NCBI_BLAST_URL, data=message, headers=headers)

        response.raise_for_status()

        soup = BeautifulSoup(response.content, "html.parser")
        rid = soup.find("input", {"name": "RID"})['value']

        finished = [False]

        schedule.every(10).seconds.do(lambda: self.check_status(rid, i, finished))

        while not finished[0]:
            schedule.run_pending()
            time.sleep(5)

    def check_status(self, rid, i, finished):
        i[0] += 1 
        params = {
            'RID': rid,
            'CMD': 'Get',
            # 'FORMAT_TYPE': 'XML'
        }
        response = requests.get(self.NCBI_BLAST_URL, params=params)
        soup = BeautifulSoup(response.content, "html.parser")
        table_stat = soup.find("table", {"id": "statInfo"})
        table_data = soup.find("table", {"id": "dscTable"})

        if table_data:
            self.extract_data(table_data) 
            finished[0] = True

        elif table_stat:
            rows = table_stat.find_all("tr")

            for row in rows:
                cells = row.find_all("td")
                if len(cells) > 1 and cells[0].get_text(strip=True) == "Status":
                    status_value = cells[1].get_text(strip=True)
                    break

            print(f"Status: {status_value}, запрос {i[0]}")
            if status_value == "Searching":
                pass
                # schedule.every(10).seconds.do(lambda: self.check_status(rid, i, finished))


    def extract_data(self, table_data):
        rows = table_data.find_all("tr")
        for row in rows:
            id_cell = row.find("td", {"class": "l c0"})
            function_cell = row.find("td", {"class": "ellipsis c2"})
            coincidence_cell = row.find("td", {"class": "c6"})
            if id_cell and function_cell and coincidence_cell is not None:
                id_value = id_cell.find("input", {"name":"getSeqGi"})["value"]
                function_value = function_cell.find("a", class_="deflnDesc").get_text(strip=True)
                coincidence_value = coincidence_cell.get_text(strip=True)
                data = {'id': id_value, 'fucntion': function_value, 'coincidence': coincidence_value}
                self.blast_data.append(data)

        return self.blast_data


class Bedtools(SamtoolsCoverage):

    def __init__(self, df: pd.DataFrame, transcript: str):
        self.df = df
        self.transcript = transcript
        self.path_to_file = '/home/nock/VSC_Projects/alter_LINE_expr/human_genome/hg38.fa'
    

    def get_exons_positions_to_gene(self):
        filtered_gene_segments = self.df.loc[self.df['transcript'] == self.transcript]
        only_exons_df = filtered_gene_segments.loc[(self.df['name'].str.contains('exon'))]
        return only_exons_df

    def filter_exons_after_LINE(self, exons: pd.DataFrame, LINE: dict) -> pd.DataFrame:
        exons['start'] = pd.to_numeric(exons['start'], errors='coerce')
        exons['end'] = pd.to_numeric(exons['end'], errors='coerce')

        gene = exons.iloc[0]['gene']
        end_LINE = LINE[gene]['end']
        strand_LINE = LINE[gene]
        gene_stand = '+' if strand_LINE == '-' else '-'
        if gene_stand == '+':
            end_LINE = LINE[gene]['end']
            exons_after_LINE = exons[exons['start'] > end_LINE]
        elif gene_stand == '-':
            start_LINE = LINE[gene]['start']
            exons_after_LINE = exons[exons['end'] < start_LINE]
        else:
            raise ValueError(f'Undefined gene position: {gene_stand}')
        
        return exons_after_LINE
    
    def get_fasta_seq_from_exons(self, exons: list[dict]) -> list[str]:
        exons_seq = []
        for exon in exons:
            chrom = exon['chrom']
            start_exon = exon['start']
            end_exon = exon['end']
            seq = self.run_command_to_subprocces(chrom, start_exon, end_exon)
            seq = self.processing_seq(seq)
            exons_seq.append(seq)

        return exons_seq

    def processing_seq(self, fasta_seq: str):
        return fasta_seq.split('\n')[1].upper()

    def run_command_to_subprocces(self, chr, start, end):
        fasta_file = self.get_bam_filename()
        gene_position = self.get_gene_position(chr=chr, start=start, end=end)
        
        command = ["bedtools", "getfasta", "-fi", fasta_file, "-bed", "stdin"]
        result = subprocess.run(command, input=gene_position, text=True, capture_output=True)

        if result.returncode != 0:
            raise Exception(f"Ошибка при выполнении bedtools getfasta: {result.stderr}")

        return result.stdout

    def get_gene_position(self, chr, start, end):
        return f"{chr}\t{start}\t{end}\n"

class BlastManager:

    def blast(self):
        gene_segments_df = pd.read_csv('../genes/exons_and_intron_positions_in_genome.csv', sep='\t')
        bedtools = Bedtools(df=gene_segments_df, transcript='ENST00000615172.4')
        genes_exons = bedtools.get_exons_positions_to_gene()
        exons_seq = bedtools.get_fasta_seq_from_exons(exons=genes_exons.to_dict(orient='records'))
        unite_seq = ''.join(exons_seq)
        print('ДНК полностью')
        print(unite_seq)

        df = pd.read_csv('../intersect_LINE_and_genome.bed', sep='\t')
        line_processing = LINE_Processing(df=df)
        print('Фильтрация экзонов после LINE')
        exons_after_LINE = bedtools.filter_exons_after_LINE(exons=genes_exons, LINE=line_processing.get_positions_last_LINE_interval())
        exons_after_LINE = exons_after_LINE.to_dict(orient='records')
        exons_seq = bedtools.get_fasta_seq_from_exons(exons=exons_after_LINE)
        print(len(exons_seq))
        unite_seq = ''.join(exons_seq)

        print('ДНК от LINE до конца гена')
        print(unite_seq)
        blast_protein = BlastProtein(seq=unite_seq)
        blast_protein.blast()
        print(blast_protein.blast_data[0])


if __name__ == "__main__":
    blast = BlastManager()
    blast.blast()

    # sequence = "CTTTGCTGGGAAGAAAAGCAGAGGGCCTTTCCTGTCGCTTGCATCATCTTTCCTGCTCTGGTTTTATCTTGGAAACCTTGGGATCGAAGAAATTTTCCCAGTTCGTTTTCTTCTTGAGACAAGAA"
    # blast = BlastProtein(seq=sequence)
    # blast.blast()
    # print(blast.blast_data[0])