#!/usr/bin/python3

import pandas as pd

# 1) Read one random RSEM file
RSEM_file_path = "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-ACC/TCGA-OR-A5JC/RSEM_quantification/TCGA-OR-A5JC.rsem.isoforms.results"
RSEM_file = pd.read_csv(filepath_or_buffer = RSEM_file_path, sep = "\t")

# 2) Keep Length and transcript ID columns
RSEM_file_filt = RSEM_file[['transcript_id','length']]

# 3) Remove version
RSEM_file_filt.iloc[:,0] = RSEM_file_filt.iloc[:,0].str.replace(r'\..*','')

# 4) Save output
transcript_length_df_path = '/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_transcripts_length.txt'
RSEM_file_filt.to_csv(path_or_buf=transcript_length_df_path, sep='\t', header=True, index=False)
