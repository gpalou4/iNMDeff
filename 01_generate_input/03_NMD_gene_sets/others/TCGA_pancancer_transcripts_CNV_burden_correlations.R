library("dplyr")

# 1) TCGA TPM pancancer
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt", 
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print("Dimensions -->")
print(dim(TCGA_RNAseq_TPM_pancancer))
colnames(TCGA_RNAseq_TPM_pancancer) <- gsub("\\.","-",colnames(TCGA_RNAseq_TPM_pancancer))

# 2) Sample CNV burden
input_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
TCGA_CNV_burden <- read.table(file = input_path, 
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TCGA_CNV_burden <- TCGA_CNV_burden[,c("sample","CNV_burden")]
# 3) Correlations gene exp with CNV_burden
TCGA_RNAseq_TPM_pancancer_filt <- t(TCGA_RNAseq_TPM_pancancer)
# Add CNV_burden
TCGA_RNAseq_TPM_pancancer_filt <- merge(TCGA_RNAseq_TPM_pancancer_filt,TCGA_CNV_burden, by.x = "row.names", by.y = "sample", all.x = TRUE)
rownames(TCGA_RNAseq_TPM_pancancer_filt) <- TCGA_RNAseq_TPM_pancancer_filt$Row.names
TCGA_RNAseq_TPM_pancancer_filt$Row.names <- NULL
# Correlations with CNV_burden
CNV_burden_corr <- cor(TCGA_RNAseq_TPM_pancancer_filt[ , colnames(TCGA_RNAseq_TPM_pancancer_filt) != "CNV_burden"],  # Calculate correlations
                TCGA_RNAseq_TPM_pancancer_filt$CNV_burden, method = "pearson", use = "pairwise.complete.obs")
CNV_burden_corr <- data.frame(corr = CNV_burden_corr)  
CNV_burden_corr <- purity_corr[order(CNV_burden_corr$corr),, drop = FALSE]
head(CNV_burden_corr)

# Save
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_pancancer_transcripts_CNV_burden_correlation.txt")
write.table(CNV_burden_corr, file = output_path, quote=FALSE, sep='\t',row.names = TRUE, col.names = TRUE)
