import pandas as pd
import os
import pickle as pkl

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
import time

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm

from cluster import TranscriptsDelta
from utils import get_files_with_extension

OUTPUT_PATH = "../output_files/synthetic_example/"

class DeSeq(TranscriptsDelta):

    def __init__(self, cores: int):
        super().__init__()
        self.cores = cores
        self.normal_df: pd.DataFrame = None
        self.column_name_with_data = 'count_reads_with_LINE'


    def deseq(self):
        path_to_files = '../gencode47/raw_data/coverage_seed/'
        files = get_files_with_extension(path_to_files, 'csv')
        dataframes = self.load_files_to_dataframes(files, key_field=self.column_name_with_data)
        copies = self.make_dfs_copies(dataframes)
        dataframes_mod = self.drop_columns_from_dataframes(dataframes, columns_to_remove=[2,3])
        merged_df = self.merge_datafames(dataframes_mod)
        print(merged_df.head())
        columns = ['name', 'SRX', 'ER12', 'ER13']
        merged_df.columns = columns
        merged_df.set_index('name', inplace=True)
        df_transposed = merged_df.T
        df_transposed.loc[:, df_transposed.columns != 'name'] = df_transposed.loc[:, df_transposed.columns != 'name'].astype(int)

        res, size_factors = deseq2_norm(df_transposed)

        start = time.time()
        update_coverage = self.update_coverage(res, copies)
        end = time.time()
        total = end - start
        print(f'Итоговое время выполнения программы: {total:.2f} секунд')        

        self.dfs_to_files(update_coverage)


    def make_dfs_copies(self, dfs: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
        copies = {}
        for name, df in dfs.items():
            copies[name] = df.copy()
        return copies

    def normalization_expression(self, genes: pd.DataFrame, metadata: pd.DataFrame) -> DeseqDataSet:
        inference = DefaultInference(n_cpus=8)
        dds = DeseqDataSet(
            counts=genes,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            inference=inference,
            # n_cpus=8, # n_cpus can be specified here or in the inference object
        )

        dds.deseq2()
        return dds

    def get_and_write_stat(self, dds: DeseqDataSet, inference: DefaultInference = None):
        if inference is None:
            inference = DefaultInference(n_cpus=8)
        stat_res = DeseqStats(dds, inference=inference)
        stat_res.summary()
        stat_res.results_df.to_csv(os.path.join(OUTPUT_PATH, "results_my.csv"), sep='\t', index=False)

    def get_metadata(self):
        data = {'condition': ['lymphocytes', 'placenta', 'placenta']}
        df = pd.DataFrame(data, index=['SRX', 'ER12', 'ER13'])
        return df

    def update_coverage(self, df_with_norm_values: pd.DataFrame, original_dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
        for name, row in df_with_norm_values.iterrows():
            df_name = name
            original_df = original_dfs[df_name]
            for transcript, norm_coverage in row.items():
                if transcript != 'name':
                    original_df.loc[original_df['name'] == transcript, f'{name}_{self.column_name_with_data}'] = norm_coverage
                
            original_dfs[df_name] = original_df
            
        return original_dfs
    
    def update_coverage_pool(self, args: tuple[str, pd.DataFrame]) -> pd.DataFrame:
        name, df = args
        for name, row in self.normal_df.iterrows():
            for transcript, norm_coverage in row.items():
                if transcript != 'name':
                    df.loc[df['name'] == transcript, f'coverage_before_LINE_{name}'] = norm_coverage
                            
        return name, df

    def check_path_exist(self, path: str):
        os.makedirs(path, exist_ok=True)

    def dfs_to_files(self, dfs: dict[str, pd.DataFrame]):
        output_path = '../gencode47/raw_data/coverage_seed/'
        self.check_path_exist(output_path)
        for name, df in dfs.items():
            df.to_csv(os.path.join(output_path, f'{name}_upd_cov_seed.csv'), sep='\t', index=False)


if __name__ == '__main__':
    deseq = DeSeq(cores=2)
    deseq.deseq()