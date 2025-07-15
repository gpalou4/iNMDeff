rm(list=ls())

# # Libraries
# library("matrixStats")
library("readxl")
# # plots
# library("ggplot2")
# library("plyr")
# # Quantile Normalization
library("preprocessCore")
# # ComBat
library("sva")
# # PCA
# library("FactoMineR")
# library("factoextra")

# 1) Load data
# Arguments and paths

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/list_NMD_transcripts/NMD_targets_controls_outliers_PATHS_cluster.txt",
                    header = TRUE, sep = ",", stringsAsFactors = FALSE)
#TCGA.cancer <- "TCGA-LIHC"

NMD_targets_path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
CL_UPF1KD_path <- paths[paths$folder_or_object=="CL_UPF1KD","path_or_filename"]
CL_metadata_path <- paths[paths$folder_or_object=="CL_metadata_path","path_or_filename"]
RNAseq_path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor_tables_path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
TCGA_names_path <- paths[paths$folder_or_object=="TCGA_names_path","path_or_filename"]

# 1.1) Cell line RNA-seq data

CL_UPF1KD <- read.table(file = paste(CL_UPF1KD_path,paths[paths$folder_or_object=="CL_UPF1KD_TPM","path_or_filename"], sep = ""), 
                                         header = TRUE, sep = "\t", row.names = 1)
CL_UPF1KD <- CL_UPF1KD[,-1]
CL_UPF1KD <- CL_UPF1KD[grep(".*PAR.*",rownames(CL_UPF1KD), invert = TRUE),]
rownames(CL_UPF1KD) <- sub("\\..*","", rownames(CL_UPF1KD))

# Merge replicates together
CL_UPF1KD_merged <- as.data.frame(matrix(nrow = nrow(CL_UPF1KD)), row.names = rownames(CL_UPF1KD) )
CL_UPF1KD_merged[,"GSE152435_C"] <- rowMeans(CL_UPF1KD[,c(1:5)])
CL_UPF1KD_merged[,"GSE152435"] <- rowMeans(CL_UPF1KD[,c(6:10)])
CL_UPF1KD_merged[,"GSE86148_C"] <- rowMeans(CL_UPF1KD[,c(11:13)])
CL_UPF1KD_merged[,"GSE86148"] <- rowMeans(CL_UPF1KD[,c(14:16)])
CL_UPF1KD_merged[,"GSE88140"] <- rowMeans(CL_UPF1KD[,c(17:18)])
CL_UPF1KD_merged[,"GSE88148"] <- rowMeans(CL_UPF1KD[,c(19:20)])
CL_UPF1KD_merged[,"GSE88266"] <- rowMeans(CL_UPF1KD[,c(21:22)])
CL_UPF1KD_merged[,"GSE88466"] <- rowMeans(CL_UPF1KD[,c(23:24)])
CL_UPF1KD_merged <- CL_UPF1KD_merged[,-1]

# Merge HeLa cell lines?
# CL_UPF1KD_merged[,"GSE152435_GSE86148"] <- rowMeans(CL_UPF1KD_merged[,c("GSE152435","GSE86148")])
# CL_UPF1KD_merged <- CL_UPF1KD_merged[,!colnames(CL_UPF1KD_merged)%in%c("GSE152435","GSE86148")]
# head(CL_UPF1KD_merged)

# Split controls from KD
CL_controls <- c("GSE152435_C","GSE86148_C","GSE88148","GSE88266")

CL_UPF1KD_controls <- CL_UPF1KD_merged[,colnames(CL_UPF1KD_merged)%in%CL_controls]
CL_UPF1KD_merged <- CL_UPF1KD_merged[,!colnames(CL_UPF1KD_merged)%in%CL_controls]

# 1.2) TCGA RNAseq TPM and raw matrices

# 1.3) TCGA RNAseq TPM for each tissue (I still don't have this for ENSEMBL)

TCGA.tissues <- as.character(read.table(file = paste0(TCGA_names_path,paths[paths$folder_or_object=="TCGA_names","path_or_filename"]),sep = "\t")$V1)
#TCGA.tissues <- c("TCGA-LIHC")

print("1) RNAseq raw")
RNAseq_sample_names <- list()
RNAseq_TCGA_raw_all <- c()
for (i in seq(1:length(TCGA.tissues))) {
    TCGA.cancer <- as.character(TCGA.tissues[i])
    print(paste0(i," --> ",TCGA.cancer))
    RNAseq.error <- FALSE
    tryCatch( {   
    RNAseq_TCGA_raw <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq_path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])),
                                                header = TRUE, sep = "\t", row.names = 1)
    RNAseq_TCGA_raw <- RNAseq_TCGA_raw[order(rownames(RNAseq_TCGA_raw)),]
    # Filter for PRAD
    RNAseq_TCGA_raw <- t(na.omit(t(RNAseq_TCGA_raw)))
    # Round values
    RNAseq_TCGA_raw <- round(RNAseq_TCGA_raw)
    colnames(RNAseq_TCGA_raw) <- gsub("\\.","-",substr(colnames(RNAseq_TCGA_raw),1,12))
    # keep track of sample names for that cancer
    RNAseq_sample_names[[TCGA.cancer]] <- colnames(RNAseq_TCGA_raw)
    # Join matrices in one
    if (length(RNAseq_TCGA_raw_all)==0) {RNAseq_TCGA_raw_all <- RNAseq_TCGA_raw}
    else{RNAseq_TCGA_raw_all <- cbind(RNAseq_TCGA_raw_all,RNAseq_TCGA_raw)}
    # Check rownames are the same
    print(paste0("ENSEMBL transcripts IDs are the same? ",table(rownames(RNAseq_TCGA_raw)==rownames(RNAseq_TCGA_raw_all))))
    }
    ,error = function(e) {
    RNAseq.error <<- TRUE
    print(e)
    })
    if (isTRUE(RNAseq.error)) {next}
}
print("Dimensions -->")
print(dim(RNAseq_TCGA_raw_all))


# 1.1.2) RNAseq TPM
# RNAseq TPM
# RNAseq_TCGA_TPM <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq_path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
#                               header = TRUE, sep = "\t", row.names = 1)

# 1.3) Conversor tables
# 1.3.1) Protein coding transcripts
ensembl_v88_coding_transcripts <- read.table(file = paste0(conversor_tables_path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t",colClasses = "vector")

# 1.4) NMD targets genesets

# 1.4.1) NMD Colombo
NMD_Colombo <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_Colombo","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.2) NMD Karousis
NMD_Karousis <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_Karousis","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.3) NMD Tani
NMD_Tani <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_Tani","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.4) NMD Courtney
NMD_Courtney <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_Courtney","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.5) NMD ENSEMBL
NMD_Ensembl <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_ensembl","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.6) NMD global
NMD_global <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_global","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.7) NMD global all shared
NMD_global_all_shared <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_global_all_shared","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.8) NMD global 2 shared
NMD_global_2_shared <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="NMD_global_2_shared","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.9) SMG6
SMG6 <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="SMG6_ensembl","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.10) SMG7
SMG7 <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="SMG7_ensembl","path_or_filename"])), 
                                        header = TRUE, sep = "\t")
# 1.4.11) Negative control (non NMD)
non_NMD_neg_control <- read.table(file = gsub("TCGA-X",TCGA.cancer, paste0(NMD_targets_path,paths[paths$folder_or_object=="non_NMD_ensembl","path_or_filename"])), 
                                        header = TRUE, sep = "\t")

# 2) Clean data and match shared transcripts between TCGA RNAseq and rest of cell lines

# 2.1) Filter low-expressed genes and non-coding genes
# Filter TCGA RNAseq TPM matrix by the raw data matrix used in the NB regression

RNAseq_TCGA_raw_all_filt <- RNAseq_TCGA_raw_all[rowSums(log2(RNAseq_TCGA_raw_all) >= 1) >= round(length(colnames(RNAseq_TCGA_raw_all)) * 0.50),]
#RNAseq_TCGA_TPM_filt <- RNAseq_TCGA_TPM[rownames(RNAseq_TCGA_TPM)%in%rownames(RNAseq_TCGA_raw_filt),]
RNAseq_TCGA_raw_all_filt <- RNAseq_TCGA_raw_all_filt[rownames(RNAseq_TCGA_raw_all_filt)%in%ensembl_v88_coding_transcripts$V1,]

# Cell lines
CL_UPF1KD_filt <- CL_UPF1KD_merged[(rowSums(CL_UPF1KD_merged >= 0.1)) >= round(length(colnames(CL_UPF1KD_merged)) * 0.50),]
CL_UPF1KD_controls_filt <- CL_UPF1KD_controls[(rowSums(CL_UPF1KD_controls >= 0.1)) >= round(length(colnames(CL_UPF1KD_controls)) * 0.50),]

# 2.2) Share transcripts between TCGA RNAseq and Cell Lines
CL_UPF1KD_filt <- CL_UPF1KD_filt[rownames(CL_UPF1KD_filt)%in%rownames(RNAseq_TCGA_raw_all_filt),]
CL_UPF1KD_controls_filt <- CL_UPF1KD_controls_filt[rownames(CL_UPF1KD_controls_filt)%in%rownames(RNAseq_TCGA_raw_all_filt),]

# 2.3) Remove Cell Line vs tissue differences by ComBat

# 2.3.1) Merge Cell Lines and Tissues

# CL_UPF1KD_merged_all <- merge(CL_UPF1KD_filt,CL_UPF1KD_controls_filt, by = "row.names")
# rownames(CL_UPF1KD_merged_all) <- CL_UPF1KD_merged_all$Row.names
# CL_UPF1KD_merged_all$Row.names <- NULL
# TCGA_CL_merged <- merge(RNAseq_TCGA_filt,CL_UPF1KD_merged_all, by = "row.names")
# rownames(TCGA_CL_merged) <- TCGA_CL_merged$Row.names
# TCGA_CL_merged$Row.names <- NULL
# colnames(TCGA_CL_merged) <- substr(colnames(TCGA_CL_merged),1,12)

# batch <- 1:ncol(TCGA_CL_merged)
# batch[grep("TCGA",colnames(TCGA_CL_merged))] <- "TCGA"
# batch[grep("TCGA",batch, invert = TRUE)] <- "CL"
# batch <- factor(batch)

# # 2.3.2) Quantile normalization
# target <- normalize.quantiles.determine.target(x = as.matrix(TCGA_CL_merged[which(batch == "TCGA"),]))
# TCGA_CL_merged_qnorm <- normalize.quantiles.use.target(x = as.matrix(TCGA_CL_merged), target = target, copy=FALSE)

# gene_exp_dist <- function(gene_exp_medians_df, log2 = "NO") {

# if (log2 == "YES") {
    
#     p <- ggplot(log2(gene_exp_medians_df+0.001), aes(x=gene_exp_medians)) +
#     geom_histogram(aes(y=..density..), colour="black", show.legend = F, binwidth=0.5, alpha=0.75)+
#     theme(legend.position="right", plot.title = element_text(hjust = 0.5, size = 15),
#             axis.title.x = element_text(color="black", size=13, face="bold"),
#             axis.title.y = element_text(color="black", size=13, face="bold"),
#             legend.title = element_text(colour="black", size=13, face="bold")) +
#     ggtitle(paste("Gene expression distribution", sep = "")) +
#     labs(x = "log2(TPM)")
#     print(p)
# } else
# p <- ggplot(gene_exp_medians_df, aes(x=gene_exp_medians)) +
#     geom_histogram(aes(y=..density..), colour="black", show.legend = F, binwidth=0.5, alpha=0.75)+
#     theme(legend.position="right", plot.title = element_text(hjust = 0.5, size = 15),
#         axis.title.x = element_text(color="black", size=13, face="bold"),
#         axis.title.y = element_text(color="black", size=13, face="bold"),
#         legend.title = element_text(colour="black", size=13, face="bold")) +
#     ggtitle(paste("Gene expression distribution", sep = "")) +
#     labs(x = "TPM") + xlim(0,20)
# print(p)
# }

# 2.5.3) ComBat and gene expression distributions
# png(gsub("TCGA-X",TCGA.cancer,paste0(results.check.path,paths[paths$folder_or_object=="median_exp_TCGA_CL_distr","path_or_filename"])), width = 4500, height = 3500, res = 300)
# gene.exp.medians.df <- data.frame(gene_exp_medians = rowMedians(as.matrix(TCGA_CL_merged)))
# gene_exp_dist(gene_exp_medians_df = gene.exp.medians.df, log2 = "YES")
# dev.off()
# png(gsub("TCGA-X",TCGA.cancer,paste0(results.check.path,paths[paths$folder_or_object=="median_exp_qnorm_TCGA_CL_distr","path_or_filename"])), width = 4500, height = 3500, res = 300)
# gene.exp.medians.df <- data.frame(gene_exp_medians = rowMedians(as.matrix(TCGA_CL_merged.qnorm)))
# gene_exp_dist(gene_exp_medians_df = gene.exp.medians.df, log2 = "YES")
# dev.off()

# 2.3.3) ComBat
# TCGA_CL_merged_adjusted <- ComBat(dat = TCGA_CL_merged_qnorm, batch = batch, 
#                     mod = NULL, ref.batch = 'TCGA', 
#                     par.prior = TRUE)

# After combat, gene exp distribution is the same
# png(gsub("TCGA-X",TCGA.cancer,paste0(results.check.path,paths[paths$folder_or_object=="median_exp_qnorm_combat_TCGA_CL_distr","path_or_filename"])), width = 4500, height = 3500, res = 300)
# gene.exp.medians.df <- data.frame(gene_exp_medians = rowMedians(as.matrix(TCGA_CL_merged_adjusted)))
# gene_exp_dist(gene_exp_medians_df = gene.exp.medians.df, log2 = "YES")
# dev.off()

# PCA before and after combat

# 2.5.4) PCA to check if comBat fixes cell type vs tissue differences

# png(gsub("TCGA-X",TCGA.cancer,paste0(results.check.path,paths[paths$folder_or_object=="qnorm_TCGA_CL_PCA","path_or_filename"])), width = 4500, height = 3500, res = 300)
# res.pca <- PCA(t(TCGA_CL_merged.qnorm), scale.unit = TRUE, ncp = 1000, graph = FALSE)
# fviz_pca_ind(res.pca,repel = TRUE, habillage = batch, geom = "point")
# dev.off()

# png(gsub("TCGA-X",TCGA.cancer,paste0(results.check.path,paths[paths$folder_or_object=="qnorm_combat_TCGA_CL_PCA","path_or_filename"])), width = 4500, height = 3500, res = 300)
# res.pca <- PCA(t(TCGA_CL_merged.adjusted), scale.unit = TRUE, ncp = 1000, graph = FALSE)
# fviz_pca_ind(res.pca,repel = TRUE, habillage = batch, geom = "point")
# dev.off()

# Split adjusted matrix 

# CL_UPF1KD_filt <- as.data.frame(TCGA_CL_merged_adjusted[,colnames(CL_UPF1KD_filt)])
# CL_UPF1KD_controls_filt <- as.data.frame(TCGA_CL_merged_adjusted[,colnames(CL_UPF1KD_controls_filt)])
# RNAseq_TCGA_filt <- as.data.frame(TCGA_CL_merged_adjusted[,grep("TCGA",colnames(TCGA_CL_merged_adjusted))])

# For each NMD geneset...

NMD_gene_sets_names <- c("NMD_Colombo","NMD_Karousis","NMD_Tani","NMD_Courtney","NMD_Ensembl","NMD_global",
                        "NMD_global_all_shared","NMD_global_2_shared","SMG6","SMG7","non_NMD_neg_control")
                        
for (NMD_gene_set_name in NMD_gene_sets_names) {

    NMD_gene_set <- eval(parse(text = NMD_gene_set_name))
  
    print(paste("NMD geneset ----------------->",NMD_gene_set_name))

    # For each NMD target
    outlier_genes <- c()
    for (ensembl_gene in NMD_gene_set$ensembl_gene_id) {

        NMD_gene_set_gene <- NMD_gene_set[NMD_gene_set$ensembl_gene_id%in%ensembl_gene,]
        NMD_target <- as.character(NMD_gene_set_gene[NMD_gene_set_gene$final_consensus=="NMD_target","ensembl_transcript_id"])
        control <- as.character(NMD_gene_set_gene[NMD_gene_set_gene$final_consensus=="control","ensembl_transcript_id"])

        # TCGA-cell Adjusted
        #CL_UPF1KD_gene <- CL_UPF1KD_filt[rownames(CL_UPF1KD_filt)%in%NMD_gene_set_gene$ensembl_transcript_id,]
        #CL_UPF1KD_controls_gene <- CL_UPF1KD_controls_filt[rownames(CL_UPF1KD_controls_filt)%in%NMD_gene_set_gene$ensembl_transcript_id,]

        # Non-adjusted
        CL_UPF1KD_gene <- CL_UPF1KD_merged[rownames(CL_UPF1KD_merged)%in%c(NMD_target,control),]
        CL_UPF1KD_controls_gene <- CL_UPF1KD_controls[rownames(CL_UPF1KD_controls)%in%c(NMD_target,control),]

        # NMD_ratio
        NMD_KD_ratios <- CL_UPF1KD_gene[rownames(CL_UPF1KD_gene)%in%NMD_target,] / CL_UPF1KD_gene[rownames(CL_UPF1KD_gene)%in%control,]
        NMD_KD_ratios_mean <- mean(NMD_KD_ratios[(NMD_KD_ratios!=Inf & NMD_KD_ratios!= 0)],na.rm = TRUE)
        NMD_WT_ratios <- CL_UPF1KD_controls_gene[rownames(CL_UPF1KD_controls_gene)%in%NMD_target,] / CL_UPF1KD_controls_gene[rownames(CL_UPF1KD_controls_gene)%in%control,]
        NMD_WT_ratios_mean <- mean(NMD_WT_ratios[(NMD_WT_ratios!=Inf & NMD_WT_ratios!= 0)],na.rm = TRUE)

        if (is.na(NMD_KD_ratios_mean) | is.na(NMD_WT_ratios_mean)) {
            print("cannot make ratio")
            next
        }
     
        # Two requirements:
        # 1) Is NMD target 
        # (optional, bc we already know this by NMD features, this could be as extra confirmation)
        # Look if WT ratio is lower than KD ratio
        # Bc NMD works better in WT, there is less NMD target, so less ratio
        
        # 2) NMD target expression in WT cells is lower or equal than control
        # Look if WT ratio is lower or equal than 1

        if (NMD_WT_ratios_mean >= 1) {

            print(ensembl_gene)
            outlier_genes <- c(outlier_genes,ensembl_gene)



        }

        write.table(unique(outlier_genes), file = paste0("/g/strcombio/fsupek_home/gpalou/data/NMD_project/NMD_targets/ENSEMBL_transcripts/",NMD_gene_set_name,"_outliers.txt"), sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)


    }

}









