import pandas as pd

class LINE_Processing:
    def __init__(self, df: pd.DataFrame = None) -> None:
        self.df = df

    def filter_intersect_LINE_and_genes_by_strand(self):
        """Оставляем только строки, где направление гена противоположно направлению LINE"""
        columns = ['chrom', 'start', 'end', 'gene',  'strand', 'description', 'chrLINE', 'startLINE', 'endLINE', 'LINE_ID', 'strandLINE']
        self.df.columns = columns
        filter_df = self.df.loc[(self.df["strand"] != self.df["strandLINE"])]
        return filter_df
    
    def get_positions_last_LINE_interval(self) -> dict[str, dict]:
        """Получение для генов координат вставки LINE, наиболее ближайшей к концу гена.

        Функция группирует данные по генам и фильтрует вставки LINE, которые находятся внутри генов.
        Для каждого гена определяется вставка LINE, наиболее ближайшая к концу гена, в зависимости от
        направления гена (`strand`). Если `strand` равен '+', то выбирается вставка с максимальными
        координатами, если `strand` равен '-', то выбирается вставка с минимальными координатами.

        Parameters
        ----------
        self : object
            Объект, содержащий DataFrame с данными о генах и вставках LINE.

        Returns
        -------
        LINE_coords_in_genes : 
            Словарь, где ключи — это имена генов, а значения — словари с координатами вставки LINE,
            наиболее ближайшей к концу гена. Каждый словарь содержит два ключа:
            - 'start': int — начальная координата вставки LINE.
            - 'end': int — конечная координата вставки LINE.

        Raises
        -------
        ValueError
            Если значение `strand` не равно '+' или '-'.
        """
        LINE_coords_in_genes = {}
        grouped = self.df.groupby('gene')
        for gene, group in grouped:
            strand = group['strand'].iloc[0]
            startGene = group['start'].iloc[0]
            endGene = group['end'].iloc[0]
            filter_group = group.loc[(group['startLINE'] > startGene) & (group['endLINE'] < endGene)]
            if not filter_group.empty:
                if strand == '+':
                    max_start = filter_group['startLINE'].max()
                    max_end = filter_group['endLINE'].max()
                elif strand == '-':
                    max_start = filter_group['startLINE'].min()
                    max_end = filter_group['endLINE'].min()
                else:
                    raise ValueError(f'This {strand} was not possible')
                
                LINE_coords_in_genes[gene] = {'start': max_start, 'end': max_end}

        return LINE_coords_in_genes
    
    def get_positions_LINE_seed_in_genes(self) -> dict:
        LINE_coords_in_genes = self.get_positions_last_LINE_interval()

        LINE_seed_coords_in_genes = {}
        variation = 15
        for gene in LINE_coords_in_genes:
            strand_gene = self.df[self.df['gene'] == gene].iloc[0]['strand']
            chrom = self.df[self.df['gene'] == gene].iloc[0]['chrom']
            if strand_gene == '+':
                end_LINE = LINE_coords_in_genes[gene]['end']
                LINE_seed_coords_in_genes[gene] = {'chrom': chrom, 'start': end_LINE-variation, 'end': end_LINE + variation}
            elif strand_gene == '-':
                start_LINE = LINE_coords_in_genes[gene]['start']
                LINE_seed_coords_in_genes[gene] = {'chrom': chrom, 'start': start_LINE-variation, 'end': start_LINE + variation}
            else:
                print('нелогично')
                print(strand_gene)

        return LINE_seed_coords_in_genes

    def validate_min_distance_between_exon_and_LINE(self) -> list[str]:
        """Проверка, что между интервалом LINE и экзоном более 30 нуклеотидов"""
        gene_segments_df = pd.read_csv('../genes/exons_positions_in_genome.csv', sep='\t')
        # only_exons_df = gene_segments_df.loc[(gene_segments_df['name'].str.contains('exon'))]

        LINE_genes = self.filter_by_LINE_genes()
        LINE_exons_df = gene_segments_df.loc[(gene_segments_df['gene'].isin(LINE_genes))]

        LINE_coords = self.get_positions_last_LINE_interval()
        # print(LINE_coords)
        grouped = LINE_exons_df.groupby('transcript')

        filtered_transcripts = []
        for transcript, group in grouped:
            gene = group['gene'].iloc[0]
            strand = group['strand'].iloc[0]
            for index, exon in group.iterrows():                
                start_exon = exon['start']
                end_exon = exon['end']
                start_LINE = LINE_coords[gene]['start']
                end_LINE = LINE_coords[gene]['end']

                if strand == '+':
                    if self.check_exon_coord_and_LINE(start_exon, end_exon, start_LINE, end_LINE):
                        filtered_transcripts.append(transcript)

                elif strand == '-':
                    if self.check_exon_coord_and_LINE(start_exon, end_exon, start_LINE, end_LINE, forward_strand=False):
                        filtered_transcripts.append(transcript)

        return list(set(filtered_transcripts))

    def check_exon_coord_and_LINE(self, start_exon: int, end_exon: int, start_LINE: int, end_LINE: int, forward_strand: bool = True) -> bool:
        min_distance = 35
        if forward_strand:
            return True if start_exon - end_LINE > min_distance  or end_LINE - end_exon > 0 else False
        else:
            return True if start_LINE - end_exon > min_distance or start_LINE - start_exon > 0 else False

    def filter_by_LINE_genes(self):
        intersect_LINE_and_genes_df = pd.read_csv("../new_intersect_genes_and_LINE.bed", sep='\t')
        LINE_genes = list(set(intersect_LINE_and_genes_df['gene']))
        return LINE_genes

    def get_LINE_intervals_to_genome_browser(self):
        df_to_coords = self.df[['gene', 'chrLINE', 'startLINE', 'endLINE', 'LINE_ID', 'strandLINE']]
        
        df_to_coords = df_to_coords.assign(intensive=1000, color='255,0,0')
        df_to_coords = df_to_coords.assign(thickStart=df_to_coords['startLINE'], thickEnd=df_to_coords['endLINE'])

        df_to_coords = df_to_coords[['chrLINE', 'startLINE', 'endLINE', 'gene', 'intensive', 'strandLINE', 'thickStart', 'thickEnd', 'color']]
        df_to_coords.to_csv("LINE_coords_in_browser.bed", sep='\t', index=False)

class ManagementLINE:

    def __init__(self):
        pass

    def get_LINE_seed_coords():
        df = pd.read_csv('../new_intersect_genes_and_LINE.bed', sep='\t')
        line_processing = LINE_Processing(df=df)
        LINE_seed_coords_in_genes = line_processing.get_positions_LINE_seed_in_genes()
        LINE_seed_coords_list = [{'gene': gene, **values} for gene, values in LINE_seed_coords_in_genes.items()]
        df = pd.DataFrame(LINE_seed_coords_list)
        df = df[['chrom', 'start', 'end', 'gene']]
        df.to_csv('LINE_seed_coords_in_genes.bed', index=False, sep='\t')

def remove_duplicate():
    df = pd.read_csv('../new_LINE_coords_in_browser.bed', sep='\t')
    cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']
    df.columns = cols
    print(df.head())
    df = df.drop_duplicates(subset=['name'])
    df.to_csv('../LINE_coords_in_genes.bed', sep='\t', index=False)

if __name__ == "__main__":
    # df = pd.read_csv('../intersect_genes_and_LINE.bed', sep='\t')
    # processing = LINE_Processing(df)
    # filter_df = processing.filter_intersect_LINE_and_genes_by_strand()
    # print(len(set(filter_df['gene'])))
    # filter_df.to_csv('new_intersect_genes_and_LINE.bed', sep='\t', index=False)
    # res = line_processing.validate_min_distance_between_exon_and_LINE()
    # line_processing.get_LINE_intervals_to_genome_browser()
    ManagementLINE.get_LINE_seed_coords()
    # remove_duplicate()