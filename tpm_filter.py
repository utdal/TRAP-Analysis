import pandas as pd
df = pd.read_csv(r'C:\Users\NXI220005\Desktop\raw_data_to_merged_df_for_analysis\merged_tpm_files.csv')
ALL_GENE_PATH = r'C:\Users\NXI220005\Desktop\raw_data_to_merged_df_for_analysis\all_genes.xlsx'
all_gene_df = pd.read_excel(ALL_GENE_PATH)[['Gene ID', 'gene_biotype', 'chromosome_name']]
gene_biotype_map = dict(zip(all_gene_df['Gene ID'], all_gene_df['gene_biotype']))
chromosome_name_map = dict(zip(all_gene_df['Gene ID'], all_gene_df['chromosome_name']))
df['Gene Biotype'] = df['ensmbl_gene_id'].map(gene_biotype_map)
df['Chromosome Name'] = df['ensmbl_gene_id'].map(chromosome_name_map)
df = df.loc[(df['Gene Biotype']=='protein_coding') &
            (df['Chromosome Name']!='MT') &
            (~df['gene_name'].str.contains(r'^Gm\d'))]
df.to_csv(r'C:\Users\NXI220005\Desktop\kevin_tpmn_files1.csv', index=False)
