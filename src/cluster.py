import pandas as pd
import numpy as np
from utils import get_value_from_dict, get_files_with_extension
from LINE_preparation import LINE_Processing
from statistics import mean, median
import os
from tqdm import tqdm
from colorama import Fore, Back, Style

SEP = '\n-----------------------------------------------------------------------------------------------\n'


class ExonsClusterization:
    """Кластеризация экзонов на основе их позиции относительно LINE"""
    def __init__(self, df: pd.DataFrame) -> None:
        self.exons_coverage_df = df #Датасет с экзонами и их покрытием

    def clusterization_exons(self, LINE_coords_in_genes: dict):
        coverage_before_LINE = {}
        coverage_after_LINE = {}
        exons = self.exons_coverage_df.to_dict(orient='records')
        undefined = []

        print(Fore.CYAN + Style.BRIGHT + "Кластеризация покрытий экзонов относительно LINE")
        for exon in tqdm(exons[:100]):
            coverage = get_value_from_dict(exon, 'coverage')
            transcript_name = get_value_from_dict(exon, 'transcript')
            gene_name = get_value_from_dict(exon, 'gene')
            if gene_name in LINE_coords_in_genes: #Проверка, что такой транскрипт есть в словаре транскриптов с координатами LINE
                strand = get_value_from_dict(exon, 'strand')
                LINE_coords = get_value_from_dict(LINE_coords_in_genes, gene_name)
                LINE_start = get_value_from_dict(LINE_coords, 'start')
                LINE_end = get_value_from_dict(LINE_coords, 'end')
                start_exon = get_value_from_dict(exon, 'start')
                end_exon = get_value_from_dict(exon, 'end')

                if strand == '+':
                    if start_exon > LINE_end:
                        coverage_after_LINE.setdefault(transcript_name, []).append(coverage)
                    elif end_exon < LINE_end:
                        coverage_before_LINE.setdefault(transcript_name, []).append(coverage)
                    else:
                        undefined.append(transcript_name)

                elif strand == '-':

                    if end_exon < LINE_start:
                        coverage_after_LINE.setdefault(transcript_name, []).append(coverage)
                    elif start_exon > LINE_start:
                        coverage_before_LINE.setdefault(transcript_name, []).append(coverage)
                    else:
                        undefined.append(transcript_name)


        print(Style.RESET_ALL)
        # print(f'Покрытие до LINE {coverage_before_LINE}', end=sep)
        # print(f'Покрытие после LINE {coverage_after_LINE}', end=sep)

        coverage_before_LINE = self.get_mean_value(coverage_before_LINE)
        coverage_after_LINE = self.get_mean_value(coverage_after_LINE)
        
        return coverage_before_LINE, coverage_after_LINE, undefined

    def get_mean_value(self, coverage: dict):
        return {key: mean(value) for key, value in coverage.items()}

    def clusters_to_common_dict(self, coverage_before_LINE: dict, coverage_after_LINE: dict) -> list[dict]:
        coverage_before_and_after = []
        for transcript in coverage_before_LINE.keys():
            if transcript in coverage_after_LINE:
                exons = self.exons_coverage_df.loc[self.exons_coverage_df['transcript'] == transcript]
                if exons.empty:
                    raise ValueError(f"No exons found for transcript {transcript}")                   
                coverage_before = coverage_before_LINE[transcript]
                coverage_after = coverage_after_LINE[transcript]

                gene = exons['gene'].iloc[0]
                strand_exon = exons['strand'].iloc[0]
                strand_LINE = '-' if strand_exon == '+' else '+'

                coverage_before_and_after.append({'name': transcript, 'coverage_before_LINE': coverage_before, 'coverage_after_LINE': coverage_after,
                                                    'gene': gene, 'strand_transcript': strand_exon, 'strand_LINE': strand_LINE})
        return coverage_before_and_after

    def get_stat(self, coverage_before_LINE: dict, coverage_after_LINE: dict, undefined: dict):
        """
        Получение информации о кол-ве транскриптов до и после LINE, кол-ве потерянных
        транскриптов в результате фильтрации по стартпу и стопу транскриптов
        """
        print('Кол-во транскриптов до LINE', len(coverage_before_LINE.keys()), end=SEP)
        print('Кол-во транскриптов после LINE', len(coverage_after_LINE.keys()), end=SEP)

        before_keys = set(coverage_before_LINE.keys())
        after_keys = set(coverage_after_LINE.keys())
        x = 0
        for key in before_keys:
            if key in after_keys:
                x += 1

        print(f"Количество транскриптов, которых есть в покрытии до или после: {x}", end=SEP)

        if len(undefined) > 0:
            print('Кол-во потерянных транскриптов: ', len(undefined))
            print('Как пример утерянного транскрипта: ', undefined[0], end=SEP)


class TranscriptsDelta:

    def delta(self):
        files = get_files_with_extension('../coverage_seed/', 'csv')
        dataframes = self.load_files_to_dataframes(files)
        dataframes_mod = self.drop_columns_from_dataframes(dataframes, columns_to_remove=[2,3,4])
        merged_df = self.merge_datafames(dataframes_mod)
        df_with_outlier = self.get_outlier(merged_df)
        print(df_with_outlier.head())
        df_with_mean_ER = self.apply_mean_to_df(df_with_outlier, [2, 3, 4, 5, 6, 7])
        df_with_delta = self.compute_experimental_control_delta(df_with_mean_ER, 'mean_ER_index')
        clusters_df = self.cluster_dataframe_by_delta(df_with_delta)
        placenta_df_with_genes = self.add_gene_name_to_transcript_delta(clusters_df['placenta_genes'])
        placenta_df_with_genes.to_csv('../clusters/placenta_transcripts.csv', sep='\t', index=False)

    def load_files_to_dataframes(self, files: list[str]) -> dict[str, pd.DataFrame]:
        dataframes = {}
        for file in files:
            uniq_file_name = file.split('/')[-1].split('_')[0]
            df = pd.read_csv(file, sep='\t')
            df = df.rename(columns={'index': f'index_{uniq_file_name}'})
            dataframes[uniq_file_name] = df
        return dataframes

    def drop_columns_from_dataframes(self, dataframes: dict[str, pd.DataFrame], columns_to_remove: list[str]) -> dict[str, pd.DataFrame]:
        """Удаление столбцов на основе их номера из датафрейма"""
        dataframes_with_drop_columns = {}
        for name, dataframe in dataframes.items():
            dataframe.drop(dataframe.columns[columns_to_remove], axis=1, inplace=True)
            dataframes_with_drop_columns[name] = dataframe
        
        return dataframes_with_drop_columns
    
    def merge_datafames(self, dataframes: dict[str, pd.DataFrame]) -> pd.DataFrame:
        srx_df = dataframes.pop('SRX', None)
        if srx_df is None:
            raise ValueError(f'Датафрейм с контролем отсутствует: {srx_df}')

        column_names = ['transcript', 'index_SRX']
        for name, df in dataframes.items():
            column_names.append(name)
            srx_df = srx_df.merge(df, on='name', how='inner')
        
        return srx_df

    def apply_mean_to_df(self, df: pd.DataFrame, columns_to_get_mean: list[int]) -> pd.DataFrame:
        df['mean_ER_index'] = df.iloc[:, columns_to_get_mean].mean(axis=1)
        df.drop(df.columns[columns_to_get_mean], axis=1, inplace=True)
        return df
    
    def apply_median_to_df(self, df: pd.DataFrame, columns_to_get_mean: list[int]) -> pd.DataFrame:
        df['median_ER_index'] = df.iloc[:, columns_to_get_mean].median(axis=1)
        df.drop(df.columns[columns_to_get_mean], axis=1, inplace=True)
        return df

    def compute_experimental_control_delta(self, df: pd.DataFrame, column_name_with_exp_value: str) -> pd.DataFrame:
        """Расчет дельта на основе разницы между (средним индексом опытных - индекс контроля)"""
        df["delta"] = df.apply(lambda x: x[column_name_with_exp_value] - x["index_SRX"], axis=1)    
        return df     
    
    def cluster_dataframe_by_delta(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        cluster_level = 0.1
        df['cluster'] = df['delta'].apply(lambda delta: "placenta_genes" if delta > cluster_level else "lymphocytes_genes")
        cluster_dfs = {cluster: df[df['cluster'] == cluster].drop(['cluster'], axis=1) for cluster in df['cluster'].unique()}
        return cluster_dfs

    def get_trascripts_genes(self):
        df = pd.read_csv('../genes/exons_and_intron_positions_in_genome.csv', sep='\t')
        df = df.drop_duplicates(subset='transcript')
        df.set_index('transcript', inplace=True)
        transcripts_with_genes = df.to_dict(orient='index')
        return {transcript: row['gene'] for transcript, row in transcripts_with_genes.items()}
    
    def add_gene_name_to_transcript_delta(self, df: pd.DataFrame):
        genes_transcripts = self.get_trascripts_genes()
        df['gene'] = df['name'].apply(lambda name: get_value_from_dict(genes_transcripts, name))
        return df

    def get_outlier(self, df: pd.DataFrame) -> pd.DataFrame:
        df["outlier"] = df.apply(lambda x: self.check_outlier(x), axis=1)
        return df

    def check_outlier(self, row):
        q1 = np.percentile(row[2:8], 25)
        name = row['name']
        if name == 'ENST00000426083.5':
            print(row)     
        q3 = np.percentile(row[2:8], 75)        
        IQ = q3 - q1
        index_SRX = row["index_SRX"]
        if index_SRX > q3 + IQ*1.5 or index_SRX < q1 - IQ*1.5:
            return True
        else: return False

class ManagerCluster:

    def cluster(args):
        if not args.input or not os.path.exists(args.input):
            raise ValueError("Input file with exons coverage does not exist or is not specified.")

        df = pd.read_csv('../intersect_LINE_and_genome.bed', sep='\t')
        line_processing = LINE_Processing(df=df)
        LINE_coords_in_genes = line_processing.get_positions_last_LINE_interval()

        df_with_exons_coverage = pd.read_csv(args.input, sep='\t')
        exons_cluster = ExonsClusterization(df=df_with_exons_coverage)
        coverage_before_LINE, coverage_after_LINE, undefined = exons_cluster.clusterization_exons(LINE_coords_in_genes)

        coverage_before_and_after = exons_cluster.clusters_to_common_dict(coverage_before_LINE, coverage_after_LINE)
        df = pd.DataFrame(coverage_before_and_after)
        df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
        delta = TranscriptsDelta()
        delta.delta()


