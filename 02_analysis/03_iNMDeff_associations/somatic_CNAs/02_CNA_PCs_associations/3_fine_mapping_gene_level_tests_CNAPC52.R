genes_update <- function(TCGA_CNV_genes, ensembl_v88_gene_transcript_genesymbol) {

  # Add gene symbols / ensembl gene id to TCGA_CNV_genes
  genes_symbols <- unlist(lapply(strsplit(rownames(TCGA_CNV_genes),"\\|"), function(x){
    x[1]
  }))
  ensembl_gene_id <- unlist(lapply(strsplit(rownames(TCGA_CNV_genes),"\\|"), function(x){
    x[2]
  }))
  ensembl_gene_id <- gsub("(.*)\\..*","\\1",ensembl_gene_id)
  TCGA_CNV_genes$gene_symbols <- genes_symbols
  TCGA_CNV_genes$ensembl_gene_id <- ensembl_gene_id

  # 1.8.1) Update TCGA_CNV_genes ENSEMBLE GENE IDs (90%)
  #### Genes with no ENSEMBL gene ID are UNIQUE, so let's add it from our table
  missing_rows <- which(is.na(TCGA_CNV_genes$ensembl_gene_id))
  matching_genes <- ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name %in% TCGA_CNV_genes[,"gene_symbols"],c("gene_id","gene_name")]
  matching_genes <- matching_genes[!duplicated(matching_genes),]
  TCGA_CNV_genes_1 <- TCGA_CNV_genes[missing_rows,]
  TCGA_CNV_genes_1 <- merge(TCGA_CNV_genes_1,matching_genes, by.x = "gene_symbols", by.y = "gene_name", all.x = TRUE)
  TCGA_CNV_genes_1$ensembl_gene_id <- TCGA_CNV_genes_1$gene_id
  TCGA_CNV_genes_1$gene_id <- NULL
  TCGA_CNV_genes_1 <- TCGA_CNV_genes_1[-which(duplicated(TCGA_CNV_genes_1$gene_symbols)),]
  TCGA_CNV_genes_2 <- TCGA_CNV_genes[-missing_rows,]
  TCGA_CNV_genes_2 <- merge(TCGA_CNV_genes_2,matching_genes, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_2$gene_name <- NULL
  TCGA_CNV_genes_updated <- rbind(TCGA_CNV_genes_1,TCGA_CNV_genes_2)
  # Add Chromosome and genome location
  TCGA_CNV_genes_updated <- merge(TCGA_CNV_genes_updated,ensembl_v88_gtf_filt, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_updated <- TCGA_CNV_genes_updated[!is.na(TCGA_CNV_genes_updated$genome_location),]
  TCGA_CNV_genes_updated$genome_location_str <- substr(TCGA_CNV_genes_updated$genome_location,1,8)
  TCGA_CNV_genes_updated$chr <- gsub("(.*):.*","\\1",TCGA_CNV_genes_updated$genome_location)
  # Cytogenetic locations of genes
  TCGA_CNV_genes_updated <- merge(TCGA_CNV_genes_updated, anno_gene[,c("ensembl_gene_id","band")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  TCGA_CNV_genes_updated$chr_arm <- factor(substr(TCGA_CNV_genes_updated$band,1,1))
  # Colnames
  colnames(TCGA_CNV_genes_updated) <- gsub("\\.","-",colnames(TCGA_CNV_genes_updated))
  return(TCGA_CNV_genes_updated)

}

modify_NMDeff_dataframe <- function(sample_NMDeff, dataset, scale = FALSE) {
  # Convert some columns to factors
  if (dataset == "TCGA") {
    factor_cols <- c("cancer_type","cancer_type_MSI","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                    "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
    # Remove "TCGA" from the cancer type
    sample_NMDeff$cancer_type <- gsub("TCGA-","",sample_NMDeff$cancer_type)
    sample_NMDeff$cancer_type_strat <- gsub("TCGA-","",sample_NMDeff$cancer_type_strat)
    sample_NMDeff$cancer_type_MSI <- gsub("TCGA-","",sample_NMDeff$cancer_type_MSI)
  } else if (dataset == "GTEx") {
    factor_cols <- c("tissue","sample")
  }
  sample_NMDeff[factor_cols] <- lapply(sample_NMDeff[factor_cols], factor) 
  # Rename NMD genesets for the 3 methods
  all_NMD_genesets <- c(endogenous_NMD_genesets,ASE_NMD_genesets,"NMDeff_mean")
  # Endogenous
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global_2_shared","endogenous_NMD_global_2_shared_randomized")] <- c("endogenous_NMD_Consensus","endogenous_NMD_Consensus_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global","endogenous_NMD_global_randomized")] <- c("endogenous_NMD_all","endogenous_NMD_all_randomized")
  # ASE
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.2","ASE_stopgain_0.2_randomized")] <- c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_triggering_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.01","ASE_stopgain_0.01_randomized")] <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_triggering_0.01_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.2","ASE_stopgain_NMD_evading_0.2_randomized")] <- c("ASE_PTC_NMD_evading_0.2","ASE_PTC_NMD_evading_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.01","ASE_stopgain_NMD_evading_0.01_randomized")] <- c("ASE_PTC_NMD_evading_0.01","ASE_PTC_NMD_evading_0.01_randomized")

  # Scale NMD genesets for the three methods
  # Change the sign (coefficients are reversed), so higher values means high NMDeff
  sample_NMDeff[,all_NMD_genesets] <- -sample_NMDeff[,all_NMD_genesets]
  # Scale
  if (isTRUE(scale)) {
      sample_NMDeff[,all_NMD_genesets] <- scale(sample_NMDeff[,all_NMD_genesets])
  }
  # Filter samples with low PTC number in ASE
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.2 < 3),c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","NMDeff_mean")] <- NA
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.01 < 3),c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01")] <- NA

  return(sample_NMDeff)
}

order_genome_location_by_chr <- function(genes_chr) {
  genes_chr$genome_location_chr <- as.numeric(gsub(".*:(.*)-.*","\\1",genes_chr$genome_location))
  genes_chr$chr <- as.numeric(genes_chr$chr)
  genes_chr <- genes_chr[order(genes_chr$chr, genes_chr$genome_location_chr, decreasing = FALSE),]
  return(genes_chr)
}

# # facet_nested
library(ggh4x)
library(scales)
library(biomaRt)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(corrplot)
library(ggrepel)
library(reshape2)
library(tibble)
library(ggcorrplot)

# 1) Data
# 1.1) ENSEMBL gene id and Gene Symbol conversion
ensembl_v88_gene_transcript_genesymbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_v88_gene_transcript_genesymbol <- ensembl_v88_gene_transcript_genesymbol[!duplicated(ensembl_v88_gene_transcript_genesymbol[,1:2]),1:2]

# 1.2) TCGA TPM RNA-seq data
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt"
RNAseq_TCGA_TPM_all <- read.table(file = output_path, header = TRUE, sep = "\t", row.names = 1)
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))
# Check Samples with NAs
na_counts <- colSums(is.na(RNAseq_TCGA_TPM_all))
cols_na <- names(na_counts[na_counts > 0])
RNAseq_TCGA_TPM_all <- RNAseq_TCGA_TPM_all[, !colnames(RNAseq_TCGA_TPM_all) %in% cols_na]
# Check Genes with NAs
na_counts <- rowSums(is.na(RNAseq_TCGA_TPM_filt))
rows_na <- names(na_counts[na_counts > 0])
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))

# Ensembl gene id to Gene Symbol
RNAseq_TCGA_TPM_df <- merge(RNAseq_TCGA_TPM_all,ensembl_v88_gene_transcript_genesymbol, by.x = "row.names", by.y = "gene_id", all.x = TRUE)

# 1.3 ) Cytogenetic locations of genes
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# anno_gene <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","band", "gene_biotype"),mart=ensembl, useCache = FALSE)
file_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/biomaRt_ensembl_table.txt"
anno_gene <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1.4) Co-dependency scores in DepMap (siRNA and CRISPR KO)
# Pairs of correlations between all genes and UPF1, SMG1
# For every Chr 1q Cell Line classification separately
file_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/DepMap/NMD_codependeny_scores.txt"
NMD_codependency_scores <- read.table(file = file_path, header = TRUE, sep = "\t")

# 1.5) DepMap similarity network (Co-Dep scores Normalized with Onion and RPCA) 
file_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/snf_run_rpca_7_5_5.Rdata"
load(file_path)
DepMap_snf_rpca <- w
w <- NULL

# 1.5) Sample NMD efficiencies TCGA
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

# 1.6) ENSEMBL gene annotation
# ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$type == "gene",]
ensembl_v88_gtf_filt$genome_location <- gsub("chr","",paste0(ensembl_v88_gtf_filt$seqnames,":",ensembl_v88_gtf_filt$start,"-",ensembl_v88_gtf_filt$end))
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[,c("gene_name","gene_id","genome_location")]
ensembl_v88_gtf_filt$gene_id <- gsub("(.*)\\..*","\\1",ensembl_v88_gtf_filt$gene_id)

# 1.7) TCGA CNV genes data - focal
TCGA_CNV_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/TCGA_pancancer_CNV_focal.txt")
TCGA_CNV_genes_focal <- read.table(file = TCGA_CNV_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(dim(TCGA_CNV_genes_focal))

# 1.8) NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1.9) Classification of CL based on arm-level CNV by --> Cohen Sharir Nat 2021
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/Cohen_Sharir_Nat_2021_DepMap_CNV_CL.csv"
DepMap_CNV_arm_level_type <- read.csv(file = output_path, header = TRUE, sep = ";")

# 2) Adding info
# Add chromosome band and location info to CNA data
TCGA_CNV_genes_focal_updated <- genes_update(TCGA_CNV_genes_focal,ensembl_v88_gene_transcript_genesymbol)
# Add NMD_CoDep_scores
NMD_codependency_scores_filt <- NMD_codependency_scores[NMD_codependency_scores$CL_chr1q == 0,]
TCGA_CNV_genes_focal_updated <- merge(TCGA_CNV_genes_focal_updated, NMD_codependency_scores_filt, by.x = "gene_symbols", by.y = "gene_symbol", all.x = TRUE)
# Candidate genes CNA-PC52 --> Chr 2 amplification
chr_char <- 2
TCGA_CNV_genes_focal_filt <- TCGA_CNV_genes_focal_updated %>%
                                filter(chr == chr_char) %>%
                                filter(chr_arm == "q")
# Create RNA-seq & CNA dataframe for tests
CNA_tmp <- t(TCGA_CNV_genes_focal_filt)
colnames(CNA_tmp) <- CNA_tmp[1,]
CNA_tmp <- CNA_tmp[grep("TCGA",rownames(CNA_tmp)),]
CNA_tmp <- data.frame(CNA_tmp)
# Common genes between CNA and RNA-seq
dup_genes <- RNAseq_TCGA_TPM_df$gene_name[duplicated(RNAseq_TCGA_TPM_df$gene_name)]
RNAseq_tmp <- RNAseq_TCGA_TPM_df[!RNAseq_TCGA_TPM_df$gene_name %in% dup_genes,]
common_genes <- intersect(RNAseq_tmp$gene_name, colnames(CNA_tmp))
RNAseq_tmp <- RNAseq_tmp[RNAseq_tmp$gene_name %in% common_genes,]
CNA_tmp <- CNA_tmp[,colnames(CNA_tmp) %in% common_genes]
colnames(CNA_tmp) <- paste0(colnames(CNA_tmp),"_CNA")
RNAseq_tmp <- t(RNAseq_tmp)
colnames(RNAseq_tmp) <- RNAseq_tmp["gene_name",]
RNAseq_tmp <- RNAseq_tmp[grep("TCGA",rownames(RNAseq_tmp)),]
RNAseq_tmp <- data.frame(RNAseq_tmp)
rownames(RNAseq_tmp) <- gsub("\\.","-",rownames(RNAseq_tmp))
colnames(RNAseq_tmp) <- paste0(colnames(RNAseq_tmp),"_RNAseq")
# Merge
CNA_RNAseq_df <- merge(RNAseq_tmp, CNA_tmp, by = "row.names", all.x = TRUE)
rownames(CNA_RNAseq_df) <- CNA_RNAseq_df$Row.names
CNA_RNAseq_df$Row.names <- NULL
CNA_RNAseq_df <- CNA_RNAseq_df[,order(colnames(CNA_RNAseq_df))]
# Convert to numeric
CNA_RNAseq_df <- data.frame(CNA_RNAseq_df) %>%
  mutate_all(~as.numeric(trimws(.)))

# 3) Analysis

CNA_PC <- "52"

# 3.1) Pair-wise correlations between CNA and RNA-seq

# All cancers
df <- CNA_RNAseq_df
corr_CNA_RNAseq <- seq(1, ncol(df), by = 2) %>%
  map_df(~{
    col1 <- df[[.x]]
    col2 <- df[[.x + 1]]
    cor_value <- cor(col1, col2, use = "pairwise.complete.obs") # you can adjust the "use" argument as needed
    tibble(
      col1 = names(df)[.x],
      col2 = names(df)[.x + 1],
      correlation = cor_value
    )
  }) %>% arrange(desc(correlation))
# Add Info
corr_CNA_RNAseq$gene_symbol <- gsub("(.*)\\_.*","\\1",corr_CNA_RNAseq$col1)
cols <- c("gene_symbols","chr","band","UPF1_KD_codependency_score","UPF1_KO_codependency_score","SMG1_KO_codependency_score","SMG1_KD_codependency_score")
corr_CNA_RNAseq <- merge(corr_CNA_RNAseq, TCGA_CNV_genes_focal_filt[,cols] , by.x = "gene_symbol", by.y = "gene_symbols")
# Order
corr_CNA_RNAseq <- corr_CNA_RNAseq %>%
  arrange( desc(correlation), desc(SMG1_KO_codependency_score), desc(UPF1_KD_codependency_score))
data.frame(corr_CNA_RNAseq[1:25,])
corr_CNA_RNAseq_final <- corr_CNA_RNAseq[which(abs(corr_CNA_RNAseq$correlation) >= 0.35),]

# Plot
# p <- ggplot(data = corr_CNA_RNAseq_final, aes(x = band, y = gene_symbol, fill = correlation))+
#     geom_tile(color = "white") + xlab("") + ylab("") +
#     geom_text(aes(label= round(correlation,1)))+
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#     midpoint = 0, limit = c(0,1), space = "Lab", 
#     name="Spearman\nCorrelation") +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(size = 15),
#             axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))
# output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/RNAseq_CNA_corr_top_genes_CNA_PC_",CNA_PC,".png")
# png(output_path, width = 3000, height = 3000, res = 300)
# print(p)
# dev.off()

# Convert to numeric
CNA_RNAseq_df <- data.frame(CNA_RNAseq_df) %>%
  mutate_all(~as.numeric(trimws(.)))
# Add NMDeff
CNA_RNAseq_iNMDeff <- merge(CNA_RNAseq_df,sample_NMD_efficiencies_TCGA[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2","sample")], by.x = "row.names", by.y = "sample", all.x = TRUE)
rownames(CNA_RNAseq_iNMDeff) <- CNA_RNAseq_iNMDeff$Row.names
CNA_RNAseq_iNMDeff$Row.names <- NULL

# Save
write.table(CNA_RNAseq_iNMDeff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15C_2.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(CNA_RNAseq_iNMDeff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15C_2.RData")

# 3.2) Pair-wise correlations between iNMDeff and RNA-seq
corr_RNAseq_iNMDeff <- merge(RNAseq_tmp,sample_NMD_efficiencies_TCGA[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2","sample")], by.x = "row.names", by.y = "sample", all.x = TRUE)
rownames(corr_RNAseq_iNMDeff) <- corr_RNAseq_iNMDeff$Row.names
corr_RNAseq_iNMDeff$Row.names <- NULL
# All cancers
df <- corr_RNAseq_iNMDeff
# Convert to numeric
df <- data.frame(df) %>%
  mutate_all(~as.numeric(trimws(.)))
# Correlations (take only the ones from the iNMDeff)
df <- cor(df,use = "pairwise.complete.obs", method = "pearson")
df <- data.frame(df)
corr_RNAseq_iNMDeff_final <- df[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")]
corr_df <- data.frame(t(round(corr_RNAseq_iNMDeff_final,1)))
corr_melt <- melt(as.matrix(corr_df), na.rm = TRUE)
genes_to_keep <- corr_melt[which(abs(corr_melt$value) >= 0.2),"Var2"]
corr_melt <- corr_melt[corr_melt$Var2 %in% genes_to_keep,]
corr_melt <- corr_melt[order(corr_melt$value),]
corr_melt$gene_symbol <- gsub("(.*)\\_.*","\\1",corr_melt$Var2)
corr_melt <- corr_melt[!corr_melt$Var2 %in% c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2"),]
corr_RNAseq_iNMDeff_top_genes <- unique(corr_melt$gene_symbol)

# Plot
# p <- ggplot(data = corr_melt, aes(x = Var1, y = gene_symbol, fill = value))+
#     geom_tile(color = "white") + xlab("") + ylab("") +
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#     midpoint = 0, limit = c(-1,1), space = "Lab", 
#     name="Spearman\nCorrelation") +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(size = 15),
#             axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))
# output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/RNAseq_iNMDeff_corr_top_genes_CNA_PC_",CNA_PC,".png")
# png(output_path, width = 3000, height = 3000, res = 300)
# print(p)
# dev.off()

# Save
corr_RNAseq_iNMDeff <- data.frame(corr_RNAseq_iNMDeff) %>%
  mutate_all(~as.numeric(trimws(.)))

# 3.3) Scatterplot of corr_RNAseq_CNA vs corr_RNAseq_iNMDeff
rownames(corr_RNAseq_iNMDeff_final) <- unique(gsub("(.*)\\_.*","\\1",rownames(corr_RNAseq_iNMDeff_final)))
# Merge
corr_CNA_RNAseq_iNMDeff <- merge(corr_CNA_RNAseq,corr_RNAseq_iNMDeff_final, by.x = "gene_symbol", by.y = "row.names", all.x = TRUE)
# Bands & top genes
bands <- c("q31.1","q31.2","q31.3","q32.1","q32.2","q32.3","q33.1","q33.2","q33.3","q34","q35","q36.1","q36.2","q36.3")
corr_CNA_RNAseq_iNMDeff <- corr_CNA_RNAseq_iNMDeff %>%
        mutate(selected_bands = if_else(band %in% bands, band,"other"))
colnames(corr_CNA_RNAseq_iNMDeff) <- c("gene_symbol","gene_CNA","gene_RNAseq","CNA_RNA_corr","chr","band","UPF1_KD_codependency_score",
                    "UPF1_KO_codependency_score","SMG1_KO_codependency_score","SMG1_KD_codependency_score",
                    "ETG_iNMDeff_RNA_corr","ASE_iNMDeff_RNA_corr","selected_bands")
corr_CNA_RNAseq_iNMDeff <- corr_CNA_RNAseq_iNMDeff %>%
        mutate(top_hits = if_else(CNA_RNA_corr >= 0.20 & ETG_iNMDeff_RNA_corr <= -0.10 & selected_bands != "other", gene_symbol,"NA"))
# Scatterplot
p <- ggplot(data = corr_CNA_RNAseq_iNMDeff, aes(x = CNA_RNA_corr, y = ETG_iNMDeff_RNA_corr, 
                fill = factor(selected_bands), color = SMG1_KO_codependency_score))+
    labs(fill = "Chr1 bands") +
    geom_point(size = 4, alpha = 1) + xlab("mRNA vs CNA correlation") + ylab("mRNA vs iNMDeff correlation") +
    geom_smooth(fill = "black", method = lm, se=TRUE) + 
    geom_label_repel(data = subset(corr_CNA_RNAseq_iNMDeff, gene_symbol %in% top_hits), #%in% c(NMD_genes$gene_symbol)), 
                    aes(label = gene_symbol), color = "black", size = 3, arrow = arrow(length = unit(0.01, 'npc')),
                        point.padding = 0.1, nudge_x = .02, nudge_y = .02) +
    #geom_text(aes(label= round(correlation,1)))+
    scale_color_gradient2(low = "white", mid = "grey", high = "black", 
        midpoint = 0, limit = c(-0.3,0.3), space = "Lab", 
        name="SMG1 KO CoDep score") +
    theme_classic()+ 
      theme(plot.title = element_text(hjust = 0.5, size = 45),
            axis.title.x = element_text(color="black", size=30),
            axis.text.x = element_text(size = 25),
            axis.title.y = element_text(color="black", size=30),
            axis.text.y = element_text(color="black", size=25),
            legend.title=element_text(size=20),
            legend.text=element_text(size=17),
            legend.position = "right")
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/RNAseqCNA_vs_RNAseqiNMDeff_scatterplot_CNA_PC_",CNA_PC,".png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()

# Save
write.table(corr_CNA_RNAseq_iNMDeff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(corr_CNA_RNAseq_iNMDeff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15C.RData")
# corr_CNA_RNAseq_iNMDeff <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15C.RData")

# 3.4) Matrix of correlations of CoDep scores for Candidate Hits + NMD factors or scatterplot
# Candidate Hits
candidate_genes <- as.character(na.omit(corr_CNA_RNAseq_iNMDeff[corr_CNA_RNAseq_iNMDeff$top_hits != "NA","top_hits"]))
candidate_genes <- c(candidate_genes,"FARSB") # PC52
# Control Hits (random genes from the other bands with no sig correlations, remove bands that are closer to 1q21-23.1)
set.seed(401)
tmp_df <- corr_CNA_RNAseq_iNMDeff %>%
        filter(! band %in% bands) %>%
        filter(CNA_RNA_corr < 0.2 & ETG_iNMDeff_RNA_corr > -0.10 ) %>%
        filter(selected_bands == "other")
control_genes <- sample(tmp_df$gene_symbol)[1:25]
sort(table(corr_CNA_RNAseq_iNMDeff[corr_CNA_RNAseq_iNMDeff$gene_symbol %in%control_genes,"band"]))

# NMD factors
# NMD_factors <- NMD_genes[NMD_genes$NMD_type %in% c("EJC","NMD_ER","EJC_associated","NMD_factor","DECID_complex"),"gene_symbol"]
NMD_factors <- NMD_genes[,"gene_symbol"]
final_candidate_genes <- c(candidate_genes,NMD_factors,control_genes)

# Pair-wise correlations of CoDep scores between all these set of genes
# Do this for every Chr 1q CL classification separately

# 3.4.3) Corrected and normalized CRISPR KO by RPCA-Onion

final_candidate_genes[!final_candidate_genes %in% colnames(DepMap_snf_rpca)]
cols <- colnames(DepMap_snf_rpca) %in% final_candidate_genes
rows <- rownames(DepMap_snf_rpca) %in% final_candidate_genes
DepMap_snf_rpca_filt <- DepMap_snf_rpca[rows,cols]

# Limit
DepMap_snf_rpca_filt[DepMap_snf_rpca_filt >= 0.1] <- 0.1

# Get the reordered matrix from corrplot
reordered_matrix <- corrplot(as.matrix(DepMap_snf_rpca_filt), order = "hclust")
reorder_vec <- as.character(unique(reordered_matrix$corrPos$xName))
neutral_reorder_vec <- reorder_vec

# Color candidate genes
label_colors <- rep("black", length(neutral_reorder_vec))
# label_colors[neutral_reorder_vec %in% c(candidate_genes,"PLEKHO1","PSMD4","ISG20L2","MSTO1","RPRD2","ADAR")] <- "orange"
label_colors[neutral_reorder_vec %in% c(candidate_genes)] <- "orange"
label_colors[neutral_reorder_vec %in% control_genes] <- "blue"

df <- data.frame(reordered_matrix$corr)
df$label_colors <- label_colors

# df[,!colnames(df) %in% "label_colors"] <- log2(df[,!colnames(df) %in% "label_colors"])

# Save
write.table(df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15D.txt", 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15D.RData")

plot_final <- ggcorrplot(as.matrix(df[,!colnames(df) %in% "label_colors"]), method = "square", 
                        hc.order = TRUE, type = "full", lab = FALSE,
                        title = "Similarity Matrix of CRISPR KO CoDep scores",
                        #ggtheme = ggplot2::theme_classic,
                        tl.cex = 8,
                        tl.srt= 55,
                        tl.col = label_colors,
                        show.diag = TRUE)
# Extract the data from the ggplot object
plot_data <- levels(plot_final$data$Var1)

# Match the order of label_colors to the order of labels in the plot
ordered_colors <- label_colors[match(levels(plot_final$data$Var1), rownames(df))]
plot_final <- plot_final +
        scale_fill_gradient2(limits = c(0,0.1), low = "#2A0BEE", mid = "#EDFEFE", high = "#EE0B0B") +
        labs(fill = "Corr") +
        theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(color = ordered_colors),
                axis.text.y = element_text(color = ordered_colors))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/DepMap_Onion_RPCO_norm_data_CNA_PC_",CNA_PC,".png")
png(output_path, width = 4500, height = 4500, res = 300)
print(plot_final)
dev.off()

# 3.5) Manhattan plot of individual tests from all 1q arm using CoDep scores - Onion RPCA (RPCO).
# Color the candidate hits and NMD-related genes

# Obtain a p-value for each 1q gene 
# Filt RPCO CoDep scores for genes in NMD factors and 1q genes
# NMD factors
NMD_non_factors <- NMD_genes[NMD_genes$NMD_type %in% c("EJC","NMD_ER","EJC_associated","DECID_complex"),"gene_symbol"]
NMD_factors <- NMD_genes[NMD_genes$NMD_type %in% c("NMD_factor"),"gene_symbol"]
# 1q genes
genes_chr2q <- TCGA_CNV_genes_focal_filt$gene_name
genes_chr2q_candidate_region <- TCGA_CNV_genes_focal_filt[TCGA_CNV_genes_focal_filt$band %in% bands,"gene_name"]
# Filt
cols <- colnames(DepMap_snf_rpca) %in% c(NMD_factors,genes_chr2q)
rows <- rownames(DepMap_snf_rpca) %in% c(NMD_factors,genes_chr2q)
DepMap_snf_rpca_filt <- DepMap_snf_rpca[rows,cols]
DepMap_snf_rpca_filt$gene_type <- NA
# Label genes
DepMap_snf_rpca_filt$gene_type <- ifelse(rownames(DepMap_snf_rpca_filt) %in% NMD_factors, "NMD_gene", rownames(DepMap_snf_rpca_filt))
DepMap_snf_rpca_filt$gene_type <- ifelse(DepMap_snf_rpca_filt$gene_type %in% genes_chr2q_candidate_region, "2q31.1-36.3", DepMap_snf_rpca_filt$gene_type)
DepMap_snf_rpca_filt$gene_type <- ifelse(!DepMap_snf_rpca_filt$gene_type %in% c("NMD_gene","2q31.1-36.3"), "1q_other", DepMap_snf_rpca_filt$gene_type)
sort(table(DepMap_snf_rpca_filt$gene_type))
# Create dataframe with candidate genes as columns and NMD factors/control genes as rows
# candidate_genes <- row.names(DepMap_snf_rpca_filt[DepMap_snf_rpca_filt$gene_type == "2q31.1-36.3" | rownames(DepMap_snf_rpca_filt) %in% c("INTS3","SMG5","SMG7","RBM8A"),])
candidate_genes <- row.names(DepMap_snf_rpca_filt[DepMap_snf_rpca_filt$gene_type == "2q31.1-36.3",])
cols <- colnames(DepMap_snf_rpca_filt) %in% c(candidate_genes,"gene_type")
rows <- DepMap_snf_rpca_filt$gene_type != "2q31.1-36.3"
DepMap_snf_rpca_filt <- DepMap_snf_rpca_filt[rows,cols]
# Stack
DepMap_snf_rpca_filt_stacked <- stack(DepMap_snf_rpca_filt[,!colnames(DepMap_snf_rpca_filt) %in% "gene_type"])
DepMap_snf_rpca_filt_stacked$gene_type <- DepMap_snf_rpca_filt$gene_type
DepMap_snf_rpca_filt_stacked$gene_name <- rownames(DepMap_snf_rpca_filt)
colnames(DepMap_snf_rpca_filt_stacked) <- c("correlations","candidate_gene","gene_type","gene_name")
# Order
tmp_df <- aggregate(correlations ~ candidate_gene + gene_type, data = DepMap_snf_rpca_filt_stacked, median)
tmp_df <- tmp_df[tmp_df$gene_type == "NMD_gene",]
tmp_df <- tmp_df[order(tmp_df$correlations),]
order_genes <- as.character(tmp_df$candidate_gene)
DepMap_snf_rpca_filt_stacked$candidate_gene <- factor(DepMap_snf_rpca_filt_stacked$candidate_gene, levels = order_genes)
# Obtain a p-value from a wilcox-test between mean of RPCO scores of NMD_gene vs 1q_other for each 1q gene
# Test is not two-sided, we only want to test that the means between 1q_other and NMD_gene is negative ("less")
all_final_df <- c()
final_df <- data.frame(gene = NA, p_value = NA, diff = NA)
for (gene in na.omit(unique(DepMap_snf_rpca_filt_stacked$candidate_gene))) {
        final_df$gene <- gene
        df <- DepMap_snf_rpca_filt_stacked[DepMap_snf_rpca_filt_stacked$candidate_gene == gene,]
        X <- df[df$gene_type == "1q_other","correlations"]
        Y <- df[df$gene_type == "NMD_gene","correlations"]
        final_df$mean_diff <- mean(Y) - mean(X)
        final_df$median_diff <- median(Y) - median(X)
        # Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
        # Test means that X (1q_other) group are sistematically producing lower values than Y (NMD_gene)
        # measured in terms of rank sums not means/medians (although it is more similar to medians).
        # So both groups could potentially have similar or equal means/medians
        # But still have a difference in "location shift" based on mu parameter.
        # https://stats.stackexchange.com/questions/573544/which-alternative-option-should-i-choose-for-a-wilcoxon-rank-sum-test
        res <- wilcox.test(X,Y, data = df, alternative = "less", paired = FALSE)
        final_df$p_value <- res$p.value
        if (length(all_final_df) == 0) {
                all_final_df <- final_df   
        } else {
                all_final_df <- rbind(all_final_df,final_df)
        }
}

all_final_df$p_value_FDR_adjust <- p.adjust(all_final_df$p_value, method = "fdr")
# Add genome location and gene type
all_final_df <- merge(all_final_df,TCGA_CNV_genes_focal_filt[,c("gene_symbols","chr","band","chr_arm","genome_location")], by.x = "gene", by.y = "gene_symbols", all.x = TRUE)
all_final_df <- order_genome_location_by_chr(all_final_df)
all_final_df$gene_type <- NA
all_final_df$gene_type <- ifelse(all_final_df$gene %in% c(NMD_factors,"SF3B1","CWC22","FARSB","NOP58"), "NMD_related", "chr2q")
all_final_df <- all_final_df[order(all_final_df$p_value_FDR_adjust),]
table(all_final_df$gene_type)

plot_final <- ggplot(data = all_final_df, aes(x = as.numeric(as.character(genome_location_chr)), y = -log10(p_value_FDR_adjust), #fill = as.numeric(diff),
                        colour = factor(gene_type))) +
        geom_point(size = 3) +
        geom_text_repel(data = subset(all_final_df, p_value_FDR_adjust < 0.05 | gene %in% c("SMG5","SMG7","SF3B4","INTS3","RBM8A")),
                aes(label = gene), size = 2.5, arrow = arrow(length = unit(0.005, 'npc'), type = "closed"),
                        point.padding = 0.005, nudge_x = .01, nudge_y = .003, max.overlaps = 200) +
        theme_classic()+
        labs(title = "", x = "Chr 2q31.1-36.3", y = "FDR adjusted p-value") +
        geom_hline(yintercept= -log10(0.05), size = 0.75, color = "black")
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/Manhattan_plot_chr2q_RPCO_CoDep_p_values_CNA_PC_",CNA_PC,".png")
png(output_path, width = 3000, height = 3000, res = 300)
print(plot_final)
dev.off()

# Save
write.table(all_final_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(all_final_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/SuppFig15E.RData")

# 3.6) Manhattan plot of individual tests from all 1q arm using focal CNA-NMDeff associations
# ETG
file_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/pancancer_endogenous_CGC_somatic_mut_CNV_PCs_1.txt"
ETG_iNMDeff_somatic_associations <- read.table(file = file_path, header = TRUE, sep = "\t")
ETG_iNMDeff_somatic_associations <- ETG_iNMDeff_somatic_associations[ETG_iNMDeff_somatic_associations$Gene_symbol %in% genes_chr2q,]
ETG_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust <- p.adjust(ETG_iNMDeff_somatic_associations$som_CNV_amp_pval, method = "fdr")
ETG_iNMDeff_somatic_associations <- ETG_iNMDeff_somatic_associations[order(ETG_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust),]
ETG_iNMDeff_somatic_associations[which(ETG_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust < 0.05),c("Gene_symbol","som_CNV_amp_coeff")]
# ASE
# file_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/pancancer_ASE_0.2_CGC_somatic_mut_CNV_PCs_1.txt"
# ASE_iNMDeff_somatic_associations <- read.table(file = file_path, header = TRUE, sep = "\t")
# ASE_iNMDeff_somatic_associations <- ASE_iNMDeff_somatic_associations[ASE_iNMDeff_somatic_associations$Gene_symbol %in% all_final_df$gene,]
# ASE_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust <- p.adjust(ASE_iNMDeff_somatic_associations$som_CNV_amp_pval, method = "fdr")
# ASE_iNMDeff_somatic_associations <- ASE_iNMDeff_somatic_associations[order(ASE_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust),]
# ASE_iNMDeff_somatic_associations[which(ASE_iNMDeff_somatic_associations$som_CNV_amp_pval_FDR_adjust < 0.1),c("Gene_symbol","som_CNV_amp_coeff")]
# Add genome location

ETG_iNMDeff_somatic_associations <- merge(ETG_iNMDeff_somatic_associations,TCGA_CNV_genes_focal_filt[,c("gene_symbols","chr","band","chr_arm","genome_location")], 
                                            by.x = "Gene_symbol", by.y = "gene_symbols", all.x = TRUE)
ETG_iNMDeff_somatic_associations <- order_genome_location_by_chr(ETG_iNMDeff_somatic_associations)
ETG_iNMDeff_somatic_associations <- ETG_iNMDeff_somatic_associations[!is.na(ETG_iNMDeff_somatic_associations$genome_location_chr),]

plot_final <- ggplot(data = ETG_iNMDeff_somatic_associations, aes(x = as.numeric(as.character(genome_location_chr)), y = -log10(som_CNV_amp_pval_FDR_adjust), #fill = as.numeric(diff),
                        colour = factor(band))) +
        geom_point(size = 3) +
        geom_text_repel(data = subset(ETG_iNMDeff_somatic_associations, Gene_symbol %in% c("SMG5","SMG7","SF3B4","INTS3","RBM8A","PLEKHO1","MSTO1","ADAR","PMF1","ISG20L2","GON4L","PRPF3","VPS72","PSMD4","PSMB4","RPRD2")),
                aes(label = Gene_symbol), color = "black", size = 3, arrow = arrow(length = unit(0.05, 'npc'), type = "closed"),
                        point.padding = 0.05, nudge_x = .05, nudge_y = .01, max.overlaps = 350) +
        theme_classic()+
        labs(title = "CNA vs iNMDeff gene-level associations", x = "Chr 1q", y = "raw p-value") +
        geom_hline(yintercept= -log10(0.05), size = 0.75, color = "black")
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/Manhattan_plot_chr2q_CNA_vs_iNMDeff_p_values.png"
png(output_path, width = 3000, height = 3000, res = 300)
print(plot_final)
dev.off()
