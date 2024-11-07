import pandas as pd

from concurrent.futures import ProcessPoolExecutor
import subprocess

import itertools
import os

from tqdm import tqdm
from colorama import Fore, Back, Style

from reading import ParsingDNA
from LINE_preparation import LINE_Processing
from utils import get_value_from_dict

class CoverageLINE_Seed:

    def __init__(self, df: pd.DataFrame = None, df_LINE_seed_coord: pd.DataFrame = None, path_to_bamfile: str = None) -> None:
        self.transcripts_coverage_df = df
        self.df_with_LINE_seed_coords_in_genes = df_LINE_seed_coord
        self.path_to_bamfile = path_to_bamfile
        self.line_processing = LINE_Processing(df=pd.read_csv('../new_intersect_genes_and_LINE.bed', sep='\t'))
        self.transcripts_with_LINE_exon_overlap = self.line_processing.validate_min_distance_between_exon_and_LINE()
        self.LINE_seeds = self.get_LINE_seeds('../LINE/LINE_seed_seq.fa')
        self.cores: int = 2
        self.field_name_with_coverage_before_LINE: str = 'coverage_before_LINE'
        self.field_name_with_coverage_after_LINE: str = 'coverage_after_LINE'

    def get_LINE_seeds(self, path_to_file_with_LINE_seed_coords: str = None) -> dict[str, str]:
        if not os.path.exists(path_to_file_with_LINE_seed_coords):
            raise ValueError(f'Path {path_to_file_with_LINE_seed_coords} to file with LINE seeds coords was incorrect')
        parsing_fasta = ParsingDNA()
        LINE_seedings = parsing_fasta.open_fasta_file(path_to_fasta_file=path_to_file_with_LINE_seed_coords, multi_line_format=False)
        return parsing_fasta.remove_coords_from_read_name_after_bedtools(LINE_seedings) # Удаление позиций из названия рида >ZNG1A::chr9:136119-136179


    def df_to_dict(self, index_name: str = None) -> dict:
        self.transcripts_coverage_df.set_index(index_name)
        return self.transcripts_coverage_df.to_dict(orient='index')
        
    def filter_after_LINE_gt_zero(self):
        """Фильтрация df, чтобы покрытие после вставки LINE было больше 0"""
        self.transcripts_coverage_df = self.transcripts_coverage_df[self.transcripts_coverage_df[self.field_name_with_coverage_after_LINE] > 0]


    def filter_after_LINE_greater_than_pre_LINE(self):
        """Фильтрация df, чтобы покрытие после вставки LINE было больше, чем до LINE"""
        self.transcripts_coverage_df = self.transcripts_coverage_df[self.transcripts_coverage_df[self.field_name_with_coverage_after_LINE] > self.transcripts_coverage_df[self.field_name_with_coverage_before_LINE]]

    def apply_filters_to_df_with_coverage(self):
        self.filter_after_LINE_gt_zero()
        self.filter_after_LINE_greater_than_pre_LINE()

    def coverage_seed(self):

        self.apply_filters_to_df_with_coverage()
        transcripts = list(self.transcripts_coverage_df['name'])

        split_list_with_transcripts = self.split_list_to_several(transcripts)

        print(Fore.CYAN + Style.BRIGHT + "Получение кол-во ридов с LINE транскриптов")

        with ProcessPoolExecutor(max_workers=self.cores) as executor:
            results = list(executor.map(self.multiple_coverage_seed, split_list_with_transcripts))            
        
        print(Style.RESET_ALL)

        transcripts_with_count_LINE_reads = list(itertools.chain(*results)) 
        transcripts_with_count_LINE_reads = [transcript for transcript in transcripts_with_count_LINE_reads if transcript != None]

        return transcripts_with_count_LINE_reads


    def split_list_to_several(self, list_with_all_data: list[str]) -> list[list[str]]:
        len_list = len(list_with_all_data)
        len_one_list = len_list // self.cores
        dfs = []
        x = 0
        for i in range(self.cores):
            start_pos = x
            end_pos = x + len_one_list
            dfs.append(list_with_all_data[start_pos:end_pos])
            x += len_one_list
        
        if x < len_list:
            dfs.append(list_with_all_data[x:len_list])
        
        return dfs

    def get_coverage_seed(self, transcript):
        if not transcript in self.transcripts_with_LINE_exon_overlap:
            transcript = self.transcripts_coverage_df[self.transcripts_coverage_df['name'] == transcript].iloc[0]
            name = transcript['name']
            gene = transcript['gene']
            coverage_after = transcript['coverage_after_LINE']
            seed_seq = get_value_from_dict(self.LINE_seeds, gene)
            reads = self.get_reads(gene)
            count_reads_with_seed = self.get_reads_with_LINE_seedings(reads, seed_seq)
            return {'name': name, 'count_reads_with_LINE': count_reads_with_seed, 
                                                        'count_reads': len(reads.split('\n')), 'gene': gene}

    def multiple_coverage_seed(self, part_transcripts: list[str]):
        transcripts_with_count_LINE_reads = [self.get_coverage_seed(row) for row in tqdm(part_transcripts, total=len(part_transcripts))]
        return transcripts_with_count_LINE_reads

    def get_reads_with_LINE_seedings(self, reads: str, LINE_seed: str = None) -> int:
        """Принимает список генов с ридами, возвращает гены с ридами, в которых включены 90% от последовательности затравок LINE"""
        count_reads = 0
        if len(reads.split('\n')) > 1:
            for read in reads.split('\n'):
                if len(read.split('\t')) > 8:
                    read = read.split('\t')[9]
                    if self.max_matches(seq=read, seed=LINE_seed) is True:
                        count_reads += 1
            return count_reads
        else:
            return 0

    def cut_read_name_and_sequence(self, data: str) -> list[dict]: 
        reads = data.split('\n')
        read_name_and_seq = []
        # available_flags = {99, 147, 83, 163}
        bad_flags = {77, 141}
        for read in reads:
            read_fields = read.split('\t')
            if len(read_fields) > 1 and int(read_fields[1]) not in bad_flags:
                    read_name_and_seq.append({"read_name": read_fields[0], "seq": read_fields[9]})
        return read_name_and_seq

    def max_matches(self, seq, seed):
        max_matches = 0
        for i in range(len(seq) - len(seed) + 1):
            matches = sum(1 for j in range(len(seed)) if seq[i + j].lower() == seed[j].lower())
            max_matches = max(max_matches, matches)
        
        return True if max_matches/30 >= 0.7 else False

    def blast_read_and_seeding(self, seq: str, subseq: str):
        """Сравнивает попарно рид и затравку LINE со смещением каждый раз на 1 нуклеотид к концу рида, определяет наибольший % совпадений рида с затравкой"""
        matches = []
        for i in range(len(seq)):
            matches.append(self.count_mathes(seq[i:(i+len(subseq))], subseq))

        for i in range(len(seq)):
            matches.append(self.count_mathes(seq[i:(i+len(subseq))], ''.join(reversed(subseq))))

        max_match = max(matches)
        if max_match >= 0.7:
            return True
        else: return False

    def count_mathes(self, seq1: str, seq2: str) -> int:
        """Возвращает отношение совпадений к общей длине последовательностей"""
        m = 0
        if len(seq1) == len(seq2):
            for i, nuc in enumerate(seq1):
                if nuc.lower() == seq2[i].lower():
                    m += 1
            return m/len(seq1)
        return 0

    def get_reads(self, gene: str) -> str:
        chr, start, end = self.get_region_to_samtools(gene)
        reads = self.make_command_to_samtools(chr, start, end)
        return reads

    def get_region_to_samtools(self, gene: str) -> tuple:
        LINE_coords = self.df_with_LINE_seed_coords_in_genes[self.df_with_LINE_seed_coords_in_genes['gene'] == gene].iloc[0]
        return LINE_coords['chrom'], LINE_coords['start']-5, LINE_coords['end']+5

    def make_command_to_samtools(self, chr, start, end) -> str:
        samtools = self.generate_samtools()
        bam_filename = samtools.get_bam_filename()
        gene_position = samtools.get_gene_position(chr=chr, start=start, end=end)
        command = ['samtools', 'view', bam_filename, gene_position]
        result = samtools.run_subprocess_command_and_return_result(command)
        return result 
    
    def generate_samtools(self):
        return Samtools(self.path_to_bamfile)


class Samtools:
    
    def __init__(self, path_to_bamfile: str = None) -> None:
        self.path_to_bamfile = path_to_bamfile

    def run_subprocess_command_and_return_result(self, command: list[str]):
        try:
            coverage = subprocess.run(command, check=True, stdout=subprocess.PIPE)
            return coverage.stdout.decode()
        except subprocess.CalledProcessError as e:
            return (f'Command {e.cmd} failed with error {e.returncode}')

    def get_gene_position(self, chr: str, start: int, end: int):
        return f'{chr}:{start}-{end}'
    
    def get_bam_filename(self) -> str:
        return self.path_to_bamfile

class CoverageSeedManager:

    def coverage_seed(args = None):
        path_to_file = '../coverage/SRX_coverage.csv'
        df_with_transcripts_coverage = pd.read_csv(args.input, sep='\t')
        df_with_LINE_seed_coords_in_genes = pd.read_csv('../LINE/LINE_seed_coords_in_genes.bed', sep='\t')

        coverage_seed = CoverageLINE_Seed(df=df_with_transcripts_coverage, df_LINE_seed_coord=df_with_LINE_seed_coords_in_genes, path_to_bamfile=args.bamfile)
        transcripts_with_LINE_seed_coverage = coverage_seed.coverage_seed()
        df = pd.DataFrame(transcripts_with_LINE_seed_coverage)
        filename = path_to_file.split('/')[-1].split('.')[0]
        df.to_csv(args.output, index=False, sep='\t')   

def add_coverage_to_norm_LINE_reads():
    filename = 'SRX'
    df_with_coverage = pd.read_csv(f'../gencode47/normalization_data/coverage/{filename}_upd_cov.csv', sep='\t')
    df_with_coverage.drop(df_with_coverage.columns[[1,3,4, 5]], axis=1, inplace=True)

    df_with_count_LINE_reads = pd.read_csv(f'../gencode47/normalization_data/coverage_seed/{filename}_upd_cov_seed.csv', sep='\t')

    new_order_cols = ['name', 'index', f'{filename}_count_reads_with_LINE',  'count_reads', 'coverage_after_LINE', 'gene']
    df_merge = df_with_count_LINE_reads.merge(df_with_coverage, how='inner', on='name')
    df_merge['index'] = df_merge[f'{filename}_count_reads_with_LINE'] / df_merge['coverage_after_LINE']
    df_filter_by_nan = df_merge.dropna()
    df_filter_by_nan = df_filter_by_nan[new_order_cols]
    
    df_filter_by_nan.to_csv(f'../gencode47/normalization_data/coverage_seed/with_indexes/{filename}_cov_seed_i.csv', sep='\t', index=False)

if __name__ == '__main__':
    # CoverageSeedManager.coverage_seed()
    add_coverage_to_norm_LINE_reads()