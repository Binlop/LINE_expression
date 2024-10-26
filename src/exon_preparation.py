import pandas as pd
from utils import get_value_from_dict

class USCS_Processing:

    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df
        self.gene_segments_coords = []

    def iterate_rows_of_df_pandas(self):
        for index, row in df.iterrows():
            uscs_procc.generate_gene_segments_coords(row)

    def generate_gene_segments_coords(self, row: dict):
        """Создает список из координат в геноме всех экзонов и интронов""" 
        transcript_name = self.get_transcript_name_from_row(row=row)
        start_exons = self.get_exons_coords(row=row)
        end_exons = self.get_exons_coords(row=row, field_name="exonEnds")
        count_exons = get_value_from_dict(vocab=row, key="exonCount")
        gene = get_value_from_dict(vocab=row, key="geneName")
        chr = get_value_from_dict(vocab=row, key="chrom")
        strand = get_value_from_dict(vocab=row, key="strand")
        for i in range(count_exons):
            exon_name = self.get_exon_name_from_transcript(transcript=transcript_name, exon_number=i)
            start_coord_exon = start_exons[i]
            end_coord_exon = end_exons[i]
            self.gene_segments_coords.append({
                'chrom': chr, 
                'start': start_coord_exon, 
                'end': end_coord_exon,
                'name': exon_name, 
                'transcript': transcript_name,
                'gene': gene,
                'strand': strand,
                })
            # if i < count_exons-1:
            #     if self.check_len_between_exon_and_intron(start=int(end_exons[i])+1, end=int(start_exons[i+1])-1):
            #         intron = self.make_intron(chr=chr, start=int(end_exons[i])+1, end=int(start_exons[i+1])-1, transcript_name=transcript_name, gene=gene, number_intron=i, strand=strand)
            #         self.gene_segments_coords.append(intron)
            
    def make_intron(self, chr: str,  start: int, end: int, transcript_name: str, gene: str, number_intron: int, strand: str) -> dict:
        intron_name = self.get_intron_name_from_transcript(transcript=transcript_name, intron_number=number_intron)
        return {'chrom': chr, 'start': start, 'end': end, 'name': intron_name, 'transcript': transcript_name,'gene': gene, 'strand': strand}

    def check_len_between_exon_and_intron(self, start, end):
        if end - start < 10:
            return False
        return True

    def get_exon_name_from_transcript(self, transcript: str, exon_number: int) -> str:
        return f'{transcript}_exon_{exon_number}'

    def get_intron_name_from_transcript(self, transcript: str, intron_number: int) -> str:
        return f'{transcript}_intron_{intron_number}'

    def get_transcript_name_from_row(self, row: dict, field_name: str = "name") -> str:
        transcript_name = get_value_from_dict(vocab=row, key=field_name)
        return transcript_name

    def get_exons_coords(self, row: dict, field_name: str = "exonStarts") -> list:
        coords_exons =  get_value_from_dict(vocab=row, key=field_name)
        coords_exons_as_list = coords_exons.split(',')
        return coords_exons_as_list

    def get_genes_border(self):
        genes_coords = []
        for gene, transcripts in self.df.groupby('geneName'):
            max_exons_count = transcripts['exonCount'].max()
            transcript_with_max_exons = transcripts[transcripts['exonCount'] == max_exons_count].iloc[0]
            start_gene = transcript_with_max_exons['txStart']
            end_gene = transcript_with_max_exons['txEnd']
            protein_id = transcript_with_max_exons['proteinID']
            strand = transcript_with_max_exons['strand']
            genes_coords.append({'name': gene, 'start': start_gene, 'end_gene': end_gene, 'protein_id': protein_id, 'strand': strand})
        return genes_coords

        
if __name__ == "__main__":
    df = pd.read_csv('../genes/genome.csv', sep='\t')
    # df = df[:10]
    df = df[df['chrom'].str.len() < 6]
    uscs_procc = USCS_Processing(df=df)
    genes_coords = uscs_procc.get_genes_border()
    df = pd.DataFrame(genes_coords)
    df.to_csv('gene_borders.csv', sep='\t', index=False, na_rep=None)
    # uscs_procc.iterate_rows_of_df_pandas()
    # df = pd.DataFrame(uscs_procc.gene_segments_coords)
    # df.to_csv("../genes/exons_positions_in_genome.csv", sep='\t', index=False)