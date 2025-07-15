#!/usr/bin/python3
# You need to activate the following conda --> conda activate /home/gpalou/anaconda3_envs/general

# Used
import pandas as pd
import time
import sys

# 1) Obtain the data

loading_time = time.time()

print("...1) LOADING THE DATA...")
wd = "/g/strcombio/fsupek_cancer1/gpalou/"
TCGA_cancer = sys.argv[1]
print("CANCER --> ",TCGA_cancer)

# 1.1) RNA-seq matrix

print("...RNA-SEQ...")
RNAseq_matrix_path = wd + "TCGA_RNAseq_quantification/primary_tumor/TCGA-"+TCGA_cancer+"/TCGA-"+TCGA_cancer+"_RNAseq_matrix_raw_transcript.txt"
RNAseq_matrix = pd.read_csv(filepath_or_buffer = RNAseq_matrix_path, sep = "\t", index_col=0)

print("Dimensions --> ", RNAseq_matrix.shape)
print("...COMPLETED...")

# 1.2) CNA file

print("...CNV...")
CNA_path = wd + "TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_"+TCGA_cancer+"-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt"
if (TCGA_cancer == "SKCM"):
	CNA_path = wd + "TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_"+TCGA_cancer+"-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt"
elif (TCGA_cancer == "LAML"):
    CNA_path = wd + "TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_"+TCGA_cancer+"-TB.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt"
CNA = pd.read_csv(filepath_or_buffer = CNA_path, sep = "\t", low_memory=False)
# Short the Tumor Sample Barcode to the first 16 characters, that will be used to compare both OMICs data
CNA.columns.values[3:] = CNA.columns[3:].map(lambda x: str(x)[:-12])
# Some Gene Symbols have a different format --> [Gene Symbol | ENSEMBL Gene ID] 
# Let's split them in two columns
CNA[['gene_symbol', 'ENSEMBL_gene_ID']] = CNA['Gene Symbol'].str.split(pat="|", n=-1, expand=True)
CNA["ENSEMBL_gene_ID"] = CNA["ENSEMBL_gene_ID"].str.replace(r'\..*','')
# Short the Tumor Sample Barcode to the first 12 characters
CNA.columns = CNA.columns.map(lambda x: str(x)[:12])
print("Dimensions --> ", CNA.shape)
print("...COMPLETED...")

# 1.3) WES mc3.maf file
print("...WES...")
WES = pd.read_csv(filepath_or_buffer = wd + "TCGA_WES/mc3.v0.2.8.PUBLIC.maf", sep = "\t", low_memory=False)
# Short the Tumor Sample Barcode to the first 12 characters, that will be used to compare both OMICs data
WES["Tumor_Sample_Barcode_short"] = WES["Tumor_Sample_Barcode"].map(lambda x: str(x)[:12])
print("Dimensions --> ", WES.shape)
print("...COMPLETED...")

# 1.4) ID's conversion tables
conversor_tables_path = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/"

# 1.4.1) ENSEMBL (v88) gene ID - transcript ID - Gene Symbol
ensembl_v88_gene_transcript_genesymbol = pd.read_csv(filepath_or_buffer = conversor_tables_path + "ensembl_v88_gene_transcript_genesymbol.txt", sep = "\t", low_memory=False)
print("Dimensions --> ", ensembl_v88_gene_transcript_genesymbol.shape)
# 257340

# 2) Create dictionaries to make the conversions

# 2.1) ENSEMBL gene ID - Gene Symbol (for CNA) 
# Create a dictionary with ENSEMBL transcript ID as key and Gene Symbol as value.
ensembl_v88_gene_transcript_genesymbol_filt = ensembl_v88_gene_transcript_genesymbol.dropna()
ensembl_v88_gene_transcript_genesymbol_T = ensembl_v88_gene_transcript_genesymbol_filt.set_index("transcript_id").T
ensembl_v88_gene_transcript_genesymbol_dict = ensembl_v88_gene_transcript_genesymbol_T.to_dict("list")
len(ensembl_v88_gene_transcript_genesymbol_dict)

elapsed_time_script = (time.time() - loading_time)
print("TOTAL TIME TAKEN FOR LOADING DATA: ", elapsed_time_script)

## 3) Look for transcripts that contains nonsense/frameshift/splice site mutations or have CNA ##

print("...2) SEARCHING FOR TRANSCRIPTS WITH MUTATIONS OR CNA...")
script_start = time.time()
transcripts_to_remove = {}

# For each individual of the RNAseq matrix
for sample in range(0,len(RNAseq_matrix.columns)):

    start = time.time()

    # 3.1) Filter WES dataframe by TCGA barcode and obtain all mutations/CNA from the specific individual

    # Obtain TCGA Barcode
    TCGA_barcode_sample = RNAseq_matrix.columns[sample]
    print("Inside sample:", TCGA_barcode_sample)
    mutations_to_filter = ["Nonsense_Mutation", "Nonstop_Mutation","In_Frame_Ins","Frame_Shift_Ins", "In_Frame_Del","Frame_Shift_Del", "Splice_Site"]

    # Filter WES dataframe
    WES_sample_check = False
    if ((WES["Tumor_Sample_Barcode_short"] == TCGA_barcode_sample).sum() != 0):
        WES_sample_check = True
        WES_sample = WES[WES["Tumor_Sample_Barcode_short"] == TCGA_barcode_sample]
        print("Number of total mutations in ENSEMBL transcripts in this sample:", WES_sample.shape[0])
        number_mutations = WES_sample["Variant_Classification"].isin(mutations_to_filter).sum()
        print("Number of nonsense/frameshift/splice site mutations is:",  number_mutations)
    else:
        print("This sample is not found in WES data")
    # Filter CNA dataframe
    CNA_sample_check = False
    if ((CNA.columns == TCGA_barcode_sample).sum() != 0):
        CNA_sample_check = True
        CNA_sample = CNA[["gene_symbol","ENSEMBL_gene",TCGA_barcode_sample]]
    else:
        print("This sample is not found in CNA data")

    # It's better to have separate list for each type of mutation, because the same transcript can have a nonsense + CNV, then the "key" from the dictionary will be re-updated with the last type of variation.
    # Also you have a better way to count each type of mutation for each individual
    transcripts_somatic_mutations, transcripts_CNA = [], []
    #transcripts_nonsense, transcripts_fs_ins, transcripts_fs_del, transcripts_splice_site, \
    #transcripts_CNA_amp, transcripts_CNA_del = [], [], [], [], [], []

    # 3.2) Check transcripts with nonsense/frameshift/splice site mutations if sample has mutations

    if (WES_sample_check == True):

        transcripts_somatic_mutations = list(WES_sample["Transcript_ID"])

    if (CNA_sample_check == True):

        # Check the Gene Symbol/gene ID in the remaining CNA from the filtered CNA dataframe
        #CNA_sample_del = CNA_sample[CNA_sample.iloc[:,2] < 0]
        #CNA_sample_amp = CNA_sample[CNA_sample.iloc[:,2] > 0]
        CNA_sample_amp_del = CNA_sample[CNA_sample.iloc[:,2] != 0]
        filt = ensembl_v88_gene_transcript_genesymbol["gene_name"].isin(CNA_sample_amp_del["gene_symbol"]) |  ensembl_v88_gene_transcript_genesymbol["gene_id"].isin(CNA_sample_amp_del["ENSEMBL_gene"])
        CNA_sample_filt = ensembl_v88_gene_transcript_genesymbol[filt]
        transcripts_CNA = list(CNA_sample_filt["transcript_id"].dropna())

    # Log prints for debugging
    #number_mut_to_remove = len(transcripts_nonsense) + len(transcripts_fs_ins) + len(transcripts_fs_del) + len(transcripts_splice_site)
    number_mut_to_remove = len(transcripts_somatic_mutations)
    print("Number of ENSEMBL_transcript_ID with mutations to remove in this sample")
    #print("----> Nonsense: ", len(transcripts_nonsense))
    #print("----> FS insertions: ", len(transcripts_fs_ins))
    #print("----> FS deletions: ", len(transcripts_fs_del))
    #print("----> Splice site: ", len(transcripts_splice_site))
    print("----> TOTAL: ", number_mut_to_remove)
    print("Number of ENSEMBL_transcript_ID with CNA to remove in this sample")
    #print("----> CNA amp: ", len(transcripts_CNA_amp))
    #print("----> CNA del: ", len(transcripts_CNA_del))
    print("----> CNA amp and del: ", len(transcripts_CNA))

    # 4) Create a unique list of transcripts to be removed for that specific sample

    #ENSEMBL_transcript_ID_to_remove = set(transcripts_nonsense + transcripts_fs_ins + transcripts_fs_del + transcripts_splice_site + transcripts_CNA_amp + transcripts_CNA_del)
    ENSEMBL_transcript_ID_to_remove = set(transcripts_somatic_mutations + transcripts_CNA)
    transcripts_to_remove[TCGA_barcode_sample] = list(ENSEMBL_transcript_ID_to_remove)

    elapsed_time_sample = (time.time() - start) 
    print("Time taken for this sample:", elapsed_time_sample)
    print("<--------------------NEXT SAMPLE-------------------->")

elapsed_time_script = (time.time() - script_start)
print("TOTAL TIME TAKEN FOR ALL SCRIPT: ", elapsed_time_script)

# Transform de dictionary of sample --> transcripts to remove, to a dataframe of 0 (not to remove) and 1 (to remove)

dataframe_transcripts = RNAseq_matrix

for sample in transcripts_to_remove.keys():
    print("Sample ",sample," has ",len(transcripts_to_remove[sample]), " transcripts to remove")
    dataframe_transcripts[sample] = dataframe_transcripts.index.isin(transcripts_to_remove[sample]).astype("int")

# Save the file
dataframe_transcripts.to_csv(path_or_buf= wd + "NMD_project/transcripts_to_remove/ENSEMBL/" + TCGA_cancer + "_transcripts_somatic_mut_CNV.txt", sep='\t', header=True, index=True)
#dataframe_transcripts.to_csv(path_or_buf= wd + "NMD_project/transcripts_to_remove/ENSEMBL/" + TCGA_cancer + "_transcripts_all_somatic_mut.txt", sep='\t', header=True, index=True)
