library("dplyr")

glm_PCs_CNV_gene_exp <- function(RNAseq_TPM_pancancer, TCGA_CNV_PCA_ind, NMD_geneset, CNV_PCs = 100) {
    # Num of PCs
    PCs <- c()
    for (i in 1:as.numeric(CNV_PCs)) {
        PC <- paste0("Dim.",i)
        PCs <- c(PCs,PC)
    }
    CNV_PCs_glm_char <- paste(PCs,collapse=" + ")
    # Add NMD target and NMD control to our RNAseq matrix
    RNAseq_TPM_pancancer <- filter(RNAseq_TPM_pancancer, rownames(RNAseq_TPM_pancancer) %in% c(NMD_geneset$ensembl_transcript_id))
    RNAseq_TPM_pancancer <- t(RNAseq_TPM_pancancer)
    # Merge
    RNAseq_CNV_PCs <- merge(RNAseq_TPM_pancancer,TCGA_CNV_PCA_ind[,1:CNV_PCs], by = "row.names", all.x = TRUE)
    #NMD_targets <- NMD_geneset[NMD_geneset$NMD_type_final == "NMD_target","ensembl_transcript_id"]
    NMD_transcripts <- NMD_geneset$ensembl_transcript_id
    NMD_transcripts_CNV_PCs_correlation <- data.frame(NMD_transcript = NMD_transcripts, CNV_PCs_adj_r_squared = NA)
    for (transcript in NMD_transcripts) {
        # Transcript expression
        index <- grep(paste0(transcript,"|Dim"),colnames(RNAseq_CNV_PCs))
        RNAseq_CNV_PCs_filt <- RNAseq_CNV_PCs[,index]
        colnames(RNAseq_CNV_PCs_filt)[1] <- "transcript_expression" 
        # GLM regression
        glm_char <- ""
        glm_model <- paste0("lm( transcript_expression ~ ",CNV_PCs_glm_char,", data = RNAseq_CNV_PCs_filt, na.action = na.exclude)")
        glm_res <- eval(parse(text=glm_model))
        glm_res <- summary(glm_res)
        adj_r_squared <- glm_res$adj.r.squared
        NMD_transcripts_CNV_PCs_correlation[NMD_transcripts_CNV_PCs_correlation$NMD_transcript %in% transcript,"CNV_PCs_adj_r_squared"] <- adj_r_squared
    }
    return(NMD_transcripts_CNV_PCs_correlation[order(NMD_transcripts_CNV_PCs_correlation$CNV_PCs_adj_r_squared),])
}

# 1.1) NMD targets genesets
NMD_genesets <- c("NMD_Courtney","NMD_Karousis","NMD_Colombo","NMD_Tani","NMD_ensembl","NMD_global","NMD_global_2_shared","non_NMD_genes","SMG6","SMG7")
NMD_targets_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts"

for (NMD_geneset_name in NMD_genesets) {
    if (NMD_geneset_name == "non_NMD_genes") {
        char <- "ensembl.txt"
    } else {
        char <- "ensembl_filt.txt"
    }
    NMD_geneset <- read.table(file = paste0(NMD_targets_path,"/",NMD_geneset_name,"_",char), 
                                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    NMD_geneset <- NMD_geneset[-which(NMD_geneset$tag == "PAR"),] 
    # Create variable
    assign(NMD_geneset_name, NMD_geneset, envir = parent.frame())
}
NMD_geneset <- NMD_global_2_shared

# 1.2) TCGA TPM pancancer
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt", 
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print("Dimensions -->")
print(dim(TCGA_RNAseq_TPM_pancancer))
colnames(TCGA_RNAseq_TPM_pancancer) <- gsub("\\.","-",colnames(TCGA_RNAseq_TPM_pancancer))

# 1.3) sparse-PCA from TCGA somatic CNV

TCGA_CNV_PCA_ind_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_ind.txt"
TCGA_CNV_PCA_ind <- read.table(file = TCGA_CNV_PCA_ind_path, header = TRUE, sep = "\t")
rownames(TCGA_CNV_PCA_ind) <- gsub("\\.","-",rownames(TCGA_CNV_PCA_ind))

# 2) Correlations with somatic CNV
# Regression of gene expression of NMD targets by sparse-PCs (1-20) from somatic CNV 
TCGA_NMD_transcripts_CNV_PCs_correlation <- glm_PCs_CNV_gene_exp(RNAseq_TPM_pancancer = TCGA_RNAseq_TPM_pancancer, TCGA_CNV_PCA_ind, NMD_geneset, CNV_PCs = 20)
TCGA_NMD_transcripts_CNV_PCs_correlation <- TCGA_NMD_transcripts_CNV_PCs_correlation[order(TCGA_NMD_transcripts_CNV_PCs_correlation$CNV_PCs_adj_r_squared),]

# Save
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/TCGA_NMD_transcripts_CNV_PCs_correlation.txt")
write.table(TCGA_NMD_transcripts_CNV_PCs_correlation, file = output_path, quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
