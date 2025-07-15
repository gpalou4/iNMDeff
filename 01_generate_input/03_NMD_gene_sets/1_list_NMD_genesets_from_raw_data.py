#!/opt/anaconda3/bin/python3

import pandas as pd
import xlrd

# Give the current working directory 
raw_data_path = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/raw/papers" 
output_path = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/NMD_targets/raw/"

# 1) Obtain the list of NMD target genes from each dataset

# 1.1) Tani (2012)
# All are UPF1 NMD targets

# Group A (248), indirect targets with > 2-fold upregulation but unaltered decay rates
Tani_UPF1_A = pd.read_excel(raw_data_path+"Tani_NMD_targets.xls", sheet_name="Table S2", header=0, skiprows=3)
print(Tani_UPF1_A.head())
Tani_UPF1_A["symbol"].to_csv(path_or_buf = output_path + "Tani_UPF1_ind_conf_gene_symbols.txt", sep='\t', header=False, index=False)

# Group B (709), direct targets with unchanged expression levels with > 2-fold prolonged decay rates
Tani_UPF1_B = pd.read_excel(raw_data_path+"Tani_NMD_targets.xls", sheet_name="Table S3", header=0, skiprows=3)
print(Tani_UPF1_B.head())
Tani_UPF1_B["symbol"].to_csv(path_or_buf = output_path + "Tani_UPF1_dir_unconf_gene_symbols.txt", sep='\t', header=False, index=False)

# Group C (76), direct UPF1 targets with > 2-fold upregulation with > 2-fold prolonged decay rates --> Bona-fide UPF1 targets
Tani_UPF1_C = pd.read_excel(raw_data_path+"Tani_NMD_targets.xls", sheet_name="Table S4", header=0, skiprows=3)
print(Tani_UPF1_C.head())
Tani_UPF1_C["symbol"].to_csv(path_or_buf = output_path + "Tani_UPF1_dir_conf_gene_symbols.txt", sep='\t', header=False, index=False)

# Match length of each subgroup with the paper
len(Tani_UPF1_A['symbol']) == 248
len(Tani_UPF1_B['symbol']) == 709
len(Tani_UPF1_C['symbol']) == 76 

## Combine the three spreedsheets and obtain a single column with the IDs of the three groups
Tani_UPF1 = pd.concat([Tani_UPF1_A, Tani_UPF1_B,Tani_UPF1_C])
print(Tani_UPF1["symbol"])

# 1033
Tani_UPF1["symbol"].to_csv(path_or_buf=output_path + "Tani_UPF1_all_gene_symbols.txt", sep='\t', header=False, index=False)

# 1.2) Colombo (2017)

Colombo = pd.read_excel(raw_data_path+"Colombo_NMD_targets.xlsx", header=0)
print(Colombo.head())

# 1.2.1) UPF1 NMD targets
Colombo_UPF1 = Colombo[Colombo["meta_meta"]<0.05]
# 2725
Colombo_UPF1[["gene","gene_name"]].to_csv(path_or_buf = output_path + "Colombo_UPF1_ensembl_gene_symbol.txt", sep='\t', header=True, index=False)

# 1.2.2) SMG6 NMD targets
# UPF1_FDR + meta_SMG6
Colombo_SMG6 = Colombo[(Colombo["UPF1_FDR"]<0.05) & (Colombo["meta_SMG6"]<0.05)]
# 4500 --> 1780 transcripts
Colombo_SMG6[["gene","gene_name"]].to_csv(path_or_buf = output_path + "Colombo_SMG6_ensembl_gene_symbol.txt", sep='\t', header=True, index=False)

# 1.2.3) SMG7 NMD targets
# UPF1_FDR + meta_SMG7
Colombo_SMG7 = Colombo[(Colombo["UPF1_FDR"]<0.05) & (Colombo["meta_SMG7"]<0.05)]
# 2251 --> 1047 SMG7 transcripts
Colombo_SMG7[["gene","gene_name"]].to_csv(path_or_buf = output_path + "Colombo_SMG7_ensembl_gene_symbol.txt", sep='\t', header=True, index=False)

# 1.3) Karousis (2021)
# 1911
Karousis = pd.read_excel(raw_data_path+"Karousis_NMD_targets.xlsx", sheet_name="Sup. table 1 NMD isoform pairs", header=0)
Karousis[["NMD_isoforms","non_NMD_isoforms","gene_id","gene_name"]].to_csv(path_or_buf = output_path + "Karousis_ensembl_gene_symbol.txt", sep='\t', header=True, index=False)

# 1.4) Courtney (2020)
# 2793
Courtney = pd.read_csv(raw_data_path+"Courtney_NMD_targets.txt", sep="\t")
Courtney[["Gene_Name","Ensembl_Gene_ID"]].to_csv(path_or_buf = output_path + "Courtney_ensembl_gene_symbol.txt", sep='\t', header=True, index=False)

# 1.5) Smichdt (2015)
#233
Schmidt = pd.read_excel(raw_data_path+"Schmidt_NMD_targets.xlsx", sheet_name="Sheet1", header=0, skiprows=2)
Schmidt_filt = Schmidt.loc[1:417,["Gene","NMD Trigger Feature17","Category18"]]
Schmidt_filt.columns = ["gene_symbol","NMD_feature","NMD_likely_triggering_feature"]
Schmidt_filt2 = Schmidt_filt[~Schmidt_filt["gene_symbol"].duplicated()]
Schmidt_filt2.to_csv(path_or_buf = output_path + "Schmidt_SMG6_gene_symbol.txt", sep='\t', header=True, index=False)
