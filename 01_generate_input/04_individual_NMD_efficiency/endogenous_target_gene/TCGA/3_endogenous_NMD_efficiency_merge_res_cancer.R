
args <- commandArgs(trailingOnly=TRUE)
TCGA_cancer <- args[1]
NMD_method <- args[2]
print(TCGA_cancer)

# 1) Paths
TCGA_cancers_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/[X1]/"
TCGA_cancer_path <- gsub("\\[X1\\]",TCGA_cancer,TCGA_cancers_path)

# 2) NMD genesets

if (NMD_method == "endogenous") {
    NMD_gene_sets_names <- c("NMD_Colombo","NMD_Karousis","NMD_Tani","NMD_Courtney","NMD_ensembl","NMD_global",
                        "NMD_global_all_shared","NMD_global_2_shared","SMG6","SMG7","non_NMD_neg_control",
                        "non_NMD_neg_control_with_NMD_features")
    NMD_eff <- "NMD_global_2_shared"
}

# 2) Merge NB results from samples into one dataframe for the cancer
TCGA_samples <- list.files(TCGA_cancer_path)
common_cols <- c("sample","subtype","purity","LF")
NB_reg_cancer <- c()
for (TCGA_sample in TCGA_samples) {
    NB_reg_sample <- c()
    for (NMD_gene_set_name in NMD_gene_sets_names) {
        # Check if file exists
        NMD_geneset_NB_reg_file_path <- paste0(TCGA_cancer_path,TCGA_sample,"/NMD_",NMD_method,"_transcripts/",NMD_gene_set_name,"_NB_regression_results.txt")
        if (file.exists(NMD_geneset_NB_reg_file_path)) {
            # Open PTCs file
            NMD_geneset_NB_reg_file <- read.table(file = NMD_geneset_NB_reg_file_path, header = TRUE, sep = "\t")
            if (NMD_gene_set_name != NMD_eff) {
                NMD_geneset_NB_reg_file[,c("error","num_NMD_targets")] <- NULL
            }
            if (length(NB_reg_sample) == 0) {
                NB_reg_sample <- NMD_geneset_NB_reg_file
            } else {
                NB_reg_sample <- merge(NB_reg_sample,NMD_geneset_NB_reg_file, by.x = common_cols, by.y = common_cols, all.x = TRUE)
            }
        }
    }
# Join Sample to the cancer data frame 
NB_reg_cancer <- rbind(NB_reg_cancer,NB_reg_sample)
}

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA_cancer,"/NB_results/",NMD_method,"/",TCGA_cancer,"_",NMD_method,"_NMD_efficiencies.txt")
print(paste0("Writting results to --> ",output_path))
write.table(NB_reg_cancer, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
