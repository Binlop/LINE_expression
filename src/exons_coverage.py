import pandas as pd
import subprocess
import os
from utils import get_value_from_dict
from tqdm import tqdm
from colorama import Fore, Back, Style
from LINE_preparation import LINE_Processing
import time
from concurrent.futures import ProcessPoolExecutor

class ExonCoverage:

    def __init__(self, df: pd.DataFrame, bamfile: str, cores: int = 1) -> None:
        self.df = df
        self.bamfile = bamfile
        self.line_processing = LINE_Processing()
        self.cores = cores


    def coverage(self):
        LINE_genes = self.line_processing.filter_by_LINE_genes()
        only_exons_df = self.df.loc[(self.df['name'].str.contains('exon'))]
        filtered_df = only_exons_df.loc[(self.df['gene'].isin(LINE_genes))]
        # filtered_df = filtered_df[:1001]

        dfs = self.split_df_to_cores(filtered_df)
        start_time = time.time()
        print(Fore.CYAN + Style.BRIGHT + "Получение покрытия у экзонов")
        

        with ProcessPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(self.multiple_coverage, dfs))            

        df_with_coverage = pd.concat(results)          

        print(Style.RESET_ALL)
        end_time = time.time()
        total_time = end_time - start_time
        print(f'Итоговое время выполнения программы: {total_time:.2f} секунд')        
        return df_with_coverage

    def split_df_to_cores(self, filtered_df: pd.DataFrame) -> list[pd.DataFrame]:
        len_df = len(filtered_df)
        len_one_df = len_df // self.cores
        dfs = []
        x = 0
        for i in range(self.cores):
            start_pos = x
            end_pos = x + len_one_df
            dfs.append(filtered_df[start_pos:end_pos])
            x += len_one_df
        
        if x < len_df:
            dfs.append(filtered_df[x:len_df])
        
        return dfs

    def multiple_coverage(self, filtered_df):
        filtered_df.loc[:, 'coverage'] = [self.get_coverage_to_exon(row) for row in tqdm(filtered_df.itertuples(), total=len(filtered_df))]
        return filtered_df

        
    def get_coverage_to_exon(self, exon):
        exon = exon._asdict()
        chr = get_value_from_dict(vocab=exon, key="chrom")
        start = get_value_from_dict(vocab=exon, key="start")
        end = get_value_from_dict(vocab=exon, key="end")
        coverage = SamtoolsCoverage(path_to_bamfile=self.bamfile)
        command = coverage.make_command_to_subprocces(chr, start, end)
        result_command = coverage.run_subprocess_command_and_return_result(command)
        coverage = self.parsing_coverage_from_result_str(result_command)
        return coverage
    
    def parsing_coverage_from_result_str(self, result_command: str) -> float:
        return float(result_command.split('\n')[1].split()[6])



class SamtoolsCoverage:

    def __init__(self, path_to_bamfile: str) -> None:
        self.path_to_bamfile = path_to_bamfile

    def run_subprocess_command_and_return_result(self, command: list[str]):
        try:
            coverage = subprocess.run(command, check=True, stdout=subprocess.PIPE)
            return coverage.stdout.decode()
        except subprocess.CalledProcessError as e:
            return (f'Command {e.cmd} failed with error {e.returncode}')

    def make_command_to_subprocces(self, chr, start, end):
        bam_filename = self.get_bam_filename()
        gene_position = self.get_gene_position(chr=chr, start=start, end=end)
        return ['samtools', 'coverage', '-r', gene_position, bam_filename]
    
    def get_gene_position(self, chr: str, start: int, end: int):
        return f'{chr}:{start}-{end}'
    
    def get_bam_filename(self) -> str:
        return self.path_to_bamfile

class CoverageManager:

    def coverage(args = None):
        df = pd.read_csv('../genes/exons_and_intron_positions_in_genome.csv', sep='\t')
        if not args.bamfile or not os.path.exists(args.bamfile):
            raise ValueError("BAM file does not exist or is not specified.")
        exon_coverage = ExonCoverage(df=df, bamfile = args.bamfile, cores = args.cores)
        filtered_df = exon_coverage.coverage()
        filtered_df.to_csv(args.output, sep='\t', index=False) 

if __name__ == "__main__":
    df = pd.read_csv('../genes/exons_and_intron_positions_in_genome.csv', sep='\t')
    exon_coverage = ExonCoverage(df=df, bamfile='../../alter_LINE_expr/SRR14374304_sorted.bam', cores=4)
    filtered_df = exon_coverage.coverage()
    filtered_df.to_csv('SRX_exons_coverage.csv', sep='\t', index=False)