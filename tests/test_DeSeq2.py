import os
import pickle as pkl
import pandas as pd
from multiprocessing import Manager, Pool, Process
import time

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

def classical_deseq():

    SAVE = False  # whether to save the outputs of this notebook
    OUTPUT_PATH = "../output_files/synthetic_example/"
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    if SAVE:
        # Replace this with the path to directory where you would like results to be saved
        OUTPUT_PATH = "../output_files/synthetic_example"
        os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

    counts_df = load_example_data(
        modality="raw_counts",
        dataset="synthetic",
        debug=False,
    )

    metadata = load_example_data(
        modality="metadata",
        dataset="synthetic",
        debug=False,
    )

    samples_to_keep = ~metadata.condition.isna()
    counts_df = counts_df.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]

    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]


    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
        # n_cpus=8, # n_cpus can be specified here or in the inference object
    )

    dds.deseq2()


    if SAVE:
        with open(os.path.join(OUTPUT_PATH, "dds.pkl"), "wb") as f:
            pkl.dump(dds, f)

    print(dds)
    print(dds.varm["dispersions"])
    print(dds.varm["LFC"])

    stat_res = DeseqStats(dds, inference=inference)

    stat_res.summary()
    stat_res.results_df.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))

    if SAVE:
        with open(os.path.join(OUTPUT_PATH, "stat_results.pkl"), "wb") as f:
            pkl.dump(stat_res, f)

    stat_res.summary(lfc_null=0.1, alt_hypothesis="greaterAbs")
    stat_res.plot_MA(s=20)

def multi_cores_deseq():
    sample_dict = {
        'sample1': pd.DataFrame({'gene1': [10, 20], 'gene2': [5, 15]}),
        'sample2': pd.DataFrame({'gene1': [30, 40], 'gene2': [10, 20]}),
        'sample3': pd.DataFrame({'gene1': [50, 60], 'gene2': [15, 25]})
    } 
    start_time = time.time()
    with Manager() as manager:
        shared_dict = manager.dict()
        processes = [Process(target=process_sample, args=(key, sample_dict[key], shared_dict)) for key in sample_dict.keys()]

        for process in processes:
            process.start()
        
        for process in processes:
            process.join()

        end_time = time.time()
        total = end_time - start_time
        print(f'Итоговое время выполнения программы: {total:.2f} секунд')        

        print(dict(shared_dict))


def process_sample(sample_name, df: pd.DataFrame, shared_dict):
    df_normalized = df / df.sum()
    print(f'Обработали {sample_name}')
    shared_dict[sample_name] = df

# if __name__ == '__main__':
#     multi_cores_deseq()


import pandas as pd
from multiprocessing import Pool, cpu_count
import time

# Пример функции для обработки DataFrame
def process_sample(args):
    sample_name, df = args
    # Здесь выполняется обработка DataFrame
    # Например, нормализация данных
    df_normalized = df / df.sum()
    return sample_name, df_normalized

# Пример словаря с DataFrame
sample_dict = {
    'sample1': pd.DataFrame({'gene1': [10, 20], 'gene2': [5, 15]}),
    'sample2': pd.DataFrame({'gene1': [30, 40], 'gene2': [10, 20]}),
    'sample3': pd.DataFrame({'gene1': [50, 60], 'gene2': [15, 25]}),
    'sample4': pd.DataFrame({'gene1': [70, 80], 'gene2': [20, 30]}),
    'sample5': pd.DataFrame({'gene1': [90, 100], 'gene2': [25, 35]}),
    'sample6': pd.DataFrame({'gene1': [110, 120], 'gene2': [30, 40]})
}

def multi_cores_deseq():
    start_time = time.time()
    
    # Создание списка кортежей аргументов для map
    args = [(sample_name, df) for sample_name, df in sample_dict.items()]

    # Создание пула процессов
    num_processes = cpu_count()  # Количество доступных процессоров
    with Pool(num_processes) as pool:
        # Применение функции process_sample к каждому DataFrame параллельно
        results = pool.map(process_sample, args)

    end_time = time.time()
    total = end_time - start_time
    print(f'Итоговое время выполнения программы: {total:.2f} секунд')

    # Объединение результатов в словарь
    processed_dict = {sample_name: df_normalized for sample_name, df_normalized in results}
    print(processed_dict)

if __name__ == '__main__':
    multi_cores_deseq()