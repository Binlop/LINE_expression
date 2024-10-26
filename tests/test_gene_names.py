from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi
from ncbi.datasets.openapi.models import V1GeneMatch

import pandas as pd
from math import ceil

def get_gene_information(gene_record: V1GeneMatch):

    if gene_record.warnings:
        print(gene_record.warnings)
    if gene_record.errors:
        print(gene_record.errors)

    # print gene metadata fields
    if gene_record.gene:

        gene_dictionary = gene_record.gene.to_dict()
        name = gene_dictionary['symbol']
        description = gene_dictionary['description']
        chr = gene_dictionary['chromosomes'][0]
        if 'genomic_ranges' in gene_dictionary:
            start = gene_dictionary['genomic_ranges'][0]['range'][0]['begin']
            end = gene_dictionary['genomic_ranges'][0]['range'][0]['end']
            strand = '-' if gene_dictionary['orientation'] == 'minus' else '+'
            return {'name': name, 'chr': chr, 'start': start, 'end': end, 'strand': strand, 'description': description}


def get_info():

    # Provide 1 taxon - may be an NCBI taxid, scientific name or common name
    taxon = "human"

    genes = list(pd.read_csv('../src/gene_borders.csv', sep='\t')['name'].unique())
    # Provide gene identifiers as a list of gene symbols
    gene_symbols: List[str] = genes

    with DatasetsApiClient() as api_client:
        gene_api = DatasetsGeneApi(api_client)
        try:
            genes_data = []
            step = 300
            counts = ceil(len(genes)/step)
            print(counts)
            for i in range(counts):
                part_genes = genes[i*step:i*step+step]
                gene_reply = gene_api.gene_metadata_by_tax_and_symbol(part_genes, taxon)
                for gene in gene_reply.genes:
                    data = get_gene_information(gene)
                    genes_data.append(data) if data is not None else None
                print(f'Осталось {len(genes) - (step + step*i)}')
        except DatasetsApiException as e: 
            print(f"Exception when calling GeneApi: {e}\n")


    df = pd.DataFrame(genes_data)
    df.to_csv('genes_data.csv', sep='\t', index=False)

if __name__ == "__main__":
    df = pd.read_csv('genes_data.csv', sep='\t')
    df = df.iloc[:, [1,2,3,0,4,5]]
    df['chr'] = df['chr'].apply(lambda x: 'chr' + str(x))
    df.to_csv('genes_data_1.csv', sep='\t', index=False)
