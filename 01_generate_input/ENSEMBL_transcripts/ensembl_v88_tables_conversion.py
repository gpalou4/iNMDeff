#!/usr/bin/python3

# You need to activate the following conda --> conda activate /home/gpalou/anaconda3_envs/general

from gtfparse import read_gtf
import pandas as pd

# 1) Obtain GTF GENCODE v26 (ENSEMBL v88)
ENSEMBLv88_gtf_path = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf"

# 2) Read GTF
# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
ENSEMBLv88_gtf_df = read_gtf(ENSEMBLv88_gtf_path)

# 3) Obtain Coding transcripts

# 3.1) Filter transcripts and coding 
# The nonsense_mediated_decay feature is for transcript_type only, so we use gene_type = protein_coding in the gene level not transcript level
ENSEMBLv88_coding_genes_transcripts = ENSEMBLv88_gtf_df[ENSEMBLv88_gtf_df["feature"] == "transcript"]
ENSEMBLv88_coding_genes_transcripts = ENSEMBLv88_coding_genes_transcripts[ENSEMBLv88_coding_genes_transcripts['gene_type'] == "protein_coding"]
ENSEMBLv88_coding_genes_transcripts = ENSEMBLv88_coding_genes_transcripts[['transcript_id','gene_id',]]
ENSEMBLv88_coding_genes_transcripts["transcript_id"] = pd.DataFrame(ENSEMBLv88_coding_genes_transcripts["transcript_id"].str.replace(r'\..*',''))
ENSEMBLv88_coding_genes_transcripts["gene_id"] = pd.DataFrame(ENSEMBLv88_coding_genes_transcripts["gene_id"].str.replace(r'\..*',''))

# 3.2) Save output
ENSEMBLv88_coding_transcripts_path = '/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_coding_transcripts.txt'
ENSEMBLv88_coding_genes_transcripts.to_csv(path_or_buf=ENSEMBLv88_coding_transcripts_path, sep='\t', header=False, index=False)

# 4) Obtain Gene symbols to ENSEMBL conversion table
# 4.1) Filter gene symbols, ensemble gene and transcript ID
ENSEMBLv88_gtf_df_filt = ENSEMBLv88_gtf_df[['gene_id','gene_name','transcript_id']]
ENSEMBLv88_gtf_df_filt.iloc[:,0] = ENSEMBLv88_gtf_df_filt.iloc[:,0].str.replace(r'\..*','')
ENSEMBLv88_gtf_df_filt.iloc[:,2] = ENSEMBLv88_gtf_df_filt.iloc[:,2].str.replace(r'\..*','')
ENSEMBLv88_gene_symbols_conversion = ENSEMBLv88_gtf_df_filt.drop_duplicates()

# 4.2) Save output
ENSEMBLv88_gene_symbols_conversion_path = '/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt'
ENSEMBLv88_gene_symbols_conversion.to_csv(path_or_buf=ENSEMBLv88_gene_symbols_conversion_path, sep='\t', header=True, index=False)
