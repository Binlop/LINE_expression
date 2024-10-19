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

if __name__ == "__main__":
    box_plots_to_gene_indexes()