library("dplyr")

# 1) TCGA TPM pancancer
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt", 
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print("Dimensions -->")
print(dim(TCGA_RNAseq_TPM_pancancer))
colnames(TCGA_RNAseq_TPM_pancancer) <- gsub("\\.","-",colnames(TCGA_RNAseq_TPM_pancancer))

# 2) Sample Purity
TCGA_tumor_purity <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_purity.txt",
                           header = TRUE, sep = "\t")
TCGA_tumor_purity$TCGA_sample <- gsub("(.*)-01","\\1",TCGA_tumor_purity$array)
# 3) Correlations gene exp with purity
#TCGA_RNAseq_TPM_pancancer_filt <- TCGA_RNAseq_TPM_pancancer %>% filter(rownames(TCGA_RNAseq_TPM_pancancer) %in% NMD_geneset$ensembl_transcript_id)
TCGA_RNAseq_TPM_pancancer_filt <- t(TCGA_RNAseq_TPM_pancancer)
# Add purity
TCGA_RNAseq_TPM_pancancer_filt <- merge(TCGA_RNAseq_TPM_pancancer_filt,TCGA_tumor_purity[,c("TCGA_sample","purity")], by.x = "row.names", by.y = "TCGA_sample", all.x = TRUE)
rownames(TCGA_RNAseq_TPM_pancancer_filt) <- TCGA_RNAseq_TPM_pancancer_filt$Row.names
TCGA_RNAseq_TPM_pancancer_filt$Row.names <- NULL
# Correlations with purity

purity_corr <- cor(TCGA_RNAseq_TPM_pancancer_filt[ , colnames(TCGA_RNAseq_TPM_pancancer_filt) != "purity"],  # Calculate correlations
                TCGA_RNAseq_TPM_pancancer_filt$purity, method = "pearson", use = "pairwise.complete.obs")
purity_corr <- data.frame(corr = purity_corr)  
purity_corr <- purity_corr[order(purity_corr$corr),, drop = FALSE]
head(purity_corr)

# Save
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_pancancer_transcripts_purity_correlation.txt")
write.table(purity_corr, file = output_path, quote=FALSE, sep='\t',row.names = TRUE, col.names = TRUE)