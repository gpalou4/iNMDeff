
args <- commandArgs(trailingOnly=TRUE)
GTEx_tissue <- args[1]
NMD_method <- args[2]
print(GTEx_tissue)

# 1) Paths

GTEx_tissues_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/NMDeff/[X1]/"
GTEx_tissue_path <- gsub("\\[X1\\]",GTEx_tissue,GTEx_tissues_path)

# 2) NMD genesets

if (NMD_method == "endogenous") {
    NMD_gene_sets_names <- c("NMD_Colombo","NMD_Karousis","NMD_Tani","NMD_Courtney","NMD_ensembl","NMD_global",
                        "NMD_global_all_shared","NMD_global_2_shared","SMG6","SMG7","non_NMD_neg_control",
                        "non_NMD_neg_control_with_NMD_features")
    NMD_eff <- "NMD_global_2_shared"
    NMD_gene_sets_names_cols <- gsub("\\_","\\.",NMD_gene_sets_names)
    NMD_method_char <- "NMD_endogenous"
} else if (NMD_method == "PTC") {
    NMD_gene_sets_names <- c("stopgain_NMD_triggering","stopgain_NMD_evading","synonymous_NMD_triggering")
    NMD_eff <- "stopgain_NMD_triggering"
    NMD_method_char <- "PTC"
}

other_columns <- c("sample","num_NMD_targets","error")
df_cols <- length(other_columns)+(length(NMD_gene_sets_names_cols)*5)

# 2) Merge NB results from samples into one dataframe for the cancer
GTEx_samples <- list.files(GTEx_tissue_path)
GTEx_samples <- GTEx_samples[grep("GTEX",GTEx_samples)]
common_cols <- c("sample","error")
NB_reg_tissue <- c()
for (GTEx_sample in GTEx_samples) {
    print(GTEx_sample)
    NB_reg_sample <- data.frame(matrix(ncol=df_cols,nrow=0, dimnames=list(NULL, c(other_columns,NMD_gene_sets_names_cols,paste0(NMD_gene_sets_names_cols,".pval"),
                                    paste0(NMD_gene_sets_names_cols,".sd"),paste0(NMD_gene_sets_names_cols,".CI_2.5"),paste0(NMD_gene_sets_names_cols,".CI_97.5")))))
    NB_reg_sample[1,"sample"] <- GTEx_sample
    for (NMD_gene_set_name in NMD_gene_sets_names) {
        # Check if file exists
        NMD_geneset_NB_reg_file_path <- paste0(GTEx_tissue_path,GTEx_sample,"/",NMD_method_char,"_transcripts/",NMD_gene_set_name,"_NB_regression_results.txt")
        if (file.exists(NMD_geneset_NB_reg_file_path)) {
            # Open NMD endogenous file
            NMD_geneset_NB_reg_file <- read.table(file = NMD_geneset_NB_reg_file_path, header = TRUE, sep = "\t")
            if (NMD_gene_set_name != NMD_eff) {
                NMD_geneset_NB_reg_file[,c("error","num_NMD_targets")] <- NULL
            }
            NB_reg_sample[,colnames(NMD_geneset_NB_reg_file)] <- NMD_geneset_NB_reg_file
        }
    }
# Join Sample to the cancer data frame 
NB_reg_tissue <- rbind(NB_reg_tissue,NB_reg_sample)
}
# Remove NA samples
NB_reg_tissue <- NB_reg_tissue[!rowSums(is.na(NB_reg_tissue)) %in% (ncol(NB_reg_tissue)-1),]

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/tissues/",GTEx_tissue,"/NB_results/",NMD_method,"/GTEx-",GTEx_tissue,"_",NMD_method,"_NMD_efficiencies.txt")
print(paste0("Writting results to --> ",output_path))
write.table(NB_reg_tissue, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)