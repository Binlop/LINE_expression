import matplotlib.pyplot as plt
import pandas as pd

def get_transcript_indexes() -> pd.DataFrame:
    df = pd.read_csv('../clusters/combine_seed_in_transcripts.csv', sep='\t')
    return df

def make_genes_dict_with_transcript_indexes(df: pd.DataFrame) -> dict[str, dict]:
    genes_and_transcript_indexes = {}
    genes = df.groupby('gene')
    for gene, transcripts in genes:
        transcript_indexes = {trans['name']: {'control': trans.iloc[1], 'exp': list(trans.iloc[2:-2])} for _, trans in transcripts.iterrows()}
        genes_and_transcript_indexes[gene] = transcript_indexes
    return genes_and_transcript_indexes

def filtered_transcripts_by_delta(df_with_all_transcript_indexes: pd.DataFrame) -> pd.DataFrame:
    only_clustered_transcripts = pd.read_csv('../clusters/171024_placenta_transcripts.csv', sep='\t')
    filtered_df = df_with_all_transcript_indexes[df_with_all_transcript_indexes['name'].isin(only_clustered_transcripts['name'])]
    return filtered_df

def iterate_genes_to_box_plot(genes: dict[str, dict]):
    i = 0
    while i < 1:
        generate_box_plot('ICA1', genes['ICA1'])
        i += 1

def generate_box_plot(gene, transcripts):

    num_transcripts = len(transcripts)
    num_rows = int(num_transcripts**0.5) + 1
    num_cols = (num_transcripts + num_rows - 1) // num_rows
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 7))    
    fig.suptitle(f'Кол-во ридов у гена {gene}')
    red_patch = plt.Line2D([], [], color='red', marker='o', linestyle='None', markersize=10, label='Красные точки - кол-во ридов у контроля')
    fig.legend(handles=[red_patch], loc='lower right')

    if num_transcripts == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    i = 0
    for name, data in transcripts.items():
        ax = axes[i]
        ax.boxplot(data['exp'])
        ax.scatter([1], data['control'], color='red', label='Контроль')
        ax.set_title(f"Транскрипт {name}")
        ax.set_xticklabels([i+1])
        i += 1

    for j in range(num_transcripts, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    # plt.savefig(f'../images/indexes/{gene}.svg')
    plt.savefig(f'../images/seed/{gene}.svg')
    # plt.show()

def box_plots_to_gene_indexes():
    df = get_transcript_indexes()
    only_clustered_transcripts = filtered_transcripts_by_delta(df)
    genes_and_transcript_indexes = make_genes_dict_with_transcript_indexes(only_clustered_transcripts)
    iterate_genes_to_box_plot(genes_and_transcript_indexes)

def deseq_plot1():
    import matplotlib.pyplot as plt

    # Пример данных
    coverage_before = [10, 20, 30, 40, 50]
    coverage_after = [5, 15, 25, 35, 45]

    # Создание графика
    plt.figure(figsize=(10, 6))

    # График до нормализации
    plt.plot(coverage_before, label='До нормализации', marker='o')

    # График после нормализации
    plt.plot(coverage_after, label='После нормализации', marker='o')

    # Настройка графика
    plt.title('Изменение покрытия до и после нормализации')
    plt.xlabel('Индекс')
    plt.ylabel('Покрытие')
    plt.legend()
    plt.grid(True)

    # Отображение графика
    plt.show()

def plot2():
    import matplotlib.pyplot as plt

    # Пример данных
    coverage_before = [10, 20, 30, 40, 50]
    coverage_after = [5, 15, 25, 35, 45]

    # Создание графика
    plt.figure(figsize=(10, 6))

    # График рассеяния
    plt.scatter(coverage_before, coverage_after, color='blue', label='Покрытие')

    # Добавление линии равенства (y = x)
    plt.plot([min(coverage_before), max(coverage_before)], [min(coverage_before), max(coverage_before)], color='red', linestyle='--', label='y = x')

    # Настройка графика
    plt.title('Сравнение покрытия до и после нормализации')
    plt.xlabel('Покрытие до нормализации')
    plt.ylabel('Покрытие после нормализации')
    plt.legend()
    plt.grid(True)

    # Отображение графика
    plt.show()


def plot3():
    import matplotlib.pyplot as plt
    import numpy as np

    # Пример данных
    coverage_before = np.random.normal(30, 10, 100)
    coverage_after = np.random.normal(20, 5, 100)

    # Создание графика
    plt.figure(figsize=(10, 6))

    # Гистограмма до нормализации
    plt.hist(coverage_before, bins=30, alpha=0.5, label='До нормализации', color='blue')

    # Гистограмма после нормализации
    plt.hist(coverage_after, bins=30, alpha=0.5, label='После нормализации', color='green')

    # Настройка графика
    plt.title('Распределение покрытия до и после нормализации')
    plt.xlabel('Покрытие')
    plt.ylabel('Частота')
    plt.legend()
    plt.grid(True)

    # Отображение графика
    plt.show()

from coverage_seed import CoverageLINE_Seed

def filter_after_LINE_gt_zero(transcripts_coverage_df):
    """Фильтрация df, чтобы покрытие после вставки LINE было больше 0"""
    transcripts_coverage_df = transcripts_coverage_df[transcripts_coverage_df['coverage_after_LINE'] > 0]


def filter_after_LINE_greater_than_pre_LINE(transcripts_coverage_df):
    """Фильтрация df, чтобы покрытие после вставки LINE было больше, чем до LINE"""
    transcripts_coverage_df = transcripts_coverage_df[transcripts_coverage_df['coverage_after_LINE'] > transcripts_coverage_df['coverage_before_LINE']]
    return transcripts_coverage_df

def apply_filters_to_df_with_coverage(transcripts_coverage_df):
    transcripts_coverage_df = filter_after_LINE_gt_zero(transcripts_coverage_df)
    return filter_after_LINE_greater_than_pre_LINE(transcripts_coverage_df)

def plot4():
    import matplotlib.pyplot as plt

    df_norm = pd.read_csv('../gencode47/normalization_data/coverage/SRX_upd_cov.csv', sep='\t')
    df_raw = pd.read_csv('../gencode47/raw_data/coverage/SRX_coverage.csv', sep='\t')

    coverage_seed_norm = CoverageLINE_Seed(df=df_norm)
    coverage_seed_norm.apply_filters_to_df_with_coverage()
    df_norm = coverage_seed_norm.transcripts_coverage_df
    df_norm.drop(df_norm.columns[[1,3,4,5]], axis=1, inplace=True)

    coverage_seed_raw = CoverageLINE_Seed(df=df_raw)
    coverage_seed_raw.apply_filters_to_df_with_coverage()
    df_raw = coverage_seed_raw.transcripts_coverage_df
    df_raw.drop(df_raw.columns[[1,3,4,5]], axis=1, inplace=True)

    merge_df = df_norm.merge(df_raw, on='name', how='inner', suffixes=('', f'_raw'))
    print(merge_df.head())

    transcript_names = []
    coverage_norm = []
    coverage_raw = []
    # merge_df = merge_df[(merge_df['coverage_after_LINE'] < 200) & (merge_df['coverage_after_LINE_raw'] < 200)]    

    for i, row in merge_df.iterrows():
        transcript_names.append(row['name'])
        coverage_norm.append(row['coverage_after_LINE'])
        coverage_raw.append(row['coverage_after_LINE_raw'])

    # plt.figure(figsize=(10, 6))
    # plt.hist(merge_df['coverage_after_LINE'], bins=50, alpha=0.5, label='До нормализации', color='blue')
    # plt.hist(merge_df['coverage_after_LINE_raw'], bins=50, alpha=0.5, label='После нормализации', color='green')
    # plt.title('Распределение покрытия до и после нормализации')
    # plt.xlabel('Покрытие')
    # plt.ylabel('Частота')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # Создание графика
    plt.figure(figsize=(10, 6))

    # График покрытия до нормализации
    plt.plot(transcript_names, coverage_raw, label='Покрытие до нормализации', marker='o', linestyle='-', color='blue')

    # График покрытия после нормализации
    plt.plot(transcript_names, coverage_norm, label='Покрытие после нормализации', marker='o', linestyle='-', color='green')

    # Настройка графика
    plt.title('Сравнение покрытия до и после нормализации')
    plt.xlabel('Транскрипты')
    plt.ylabel('Покрытие')
    plt.legend()
    plt.grid(True)
    plt.xticks(rotation=45)  # Поворот названий транскриптов для лучшей читаемости

    # Отображение графика
    plt.tight_layout()  # Автоматическо
    plt.show()

if __name__ == "__main__":
    # box_plots_to_gene_indexes()
    plot4()
