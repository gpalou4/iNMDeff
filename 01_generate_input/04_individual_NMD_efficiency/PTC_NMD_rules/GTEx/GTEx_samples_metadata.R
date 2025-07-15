rm(list=ls())


################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################


PCA_type_RNAseq <- function(RNAseq.matrix, PCs = 4, scale = FALSE, center = FALSE) {

    PCs.names.var <- c()
    PCs.names.ind <- c()

    print("... PCA for whole tissue ...")
    # VST transformation
    RNAseq.GTEx.raw.vst <- vst(object = as.matrix(RNAseq.matrix), nsub = 10000)
    RNAseq.GTEx.raw.vst.tp <- t(RNAseq.GTEx.raw.vst)
    # Scale and center
    #RNAseq.GTEx.raw.vst.tp <- scale(RNAseq.GTEx.raw.vst.tp, center = center, scale = scale)
    # PCA
    res.pca <- PCA(RNAseq.GTEx.raw.vst.tp, scale.unit = TRUE, ncp = 1000, graph = FALSE)
    # Plot
    PCA.path <- paste0(GTEx.RNAseq.PCA.path,GTEx.tissue,"_tissue_PCA_ind.png")
    png(PCA.path, width = 4500, height = 3000, res = 300)
    p <- fviz_pca_ind(res.pca, col.ind = "cos2",
                label = "none" # hide individual labels,
                #habillage = as.factor(RNAseq.GTEx.raw.vst.tp.group$cluster)
                #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                )
    print(p)
    dev.off()
    PCs.ob <- ncol(res.pca$ind$coord)
    if (PCs.ob == 1) {
        PCs.to.keep <- 1
    } else if (PCs < PCs.ob) {
        PCs.to.keep <- PCs
    } else {
        PCs.to.keep <- PCs.ob
    }
    # Extract the results for variables (genes)
    res.pca.var <- get_pca_var(res.pca)
    # Extract the results for individuals
    res.pca.ind <- get_pca_ind(res.pca)
    print(paste0("PCs --> ", PCs.to.keep))
    # Variables
    PC.var <- "PCA.RNAseq.var"
    assign(PC.var, data.frame(res.pca.var$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
    PCs.names.var <- c(PCs.names.var, PC.var)
    # Individuals
    PC.ind <- "PCA.RNAseq.ind"
    assign(PC.ind, data.frame(res.pca.ind$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
    PCs.names.ind <- c(PCs.names.ind, PC.ind)
  return(c(PCs.names.var,PCs.names.ind))
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################
.libPaths( rev( .libPaths() ) )
library("DESeq2")
# PCA
library("factoextra")
library("FactoMineR")

##############################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.tissue <- args[2]

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
GTEx.RNAseq.PCA.path <- paths[paths$folder_or_object=="GTEx_RNAseq_PCA_path","path_or_filename"]
GTEx.metadata.path <- paths[paths$folder_or_object=="GTEx_metadata_path","path_or_filename"]

print(paste0("GTEx cancer --> ",GTEx.tissue))

# 1.1) RNAseq matrix raw counts
RNAseq.GTEx.raw <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])),
                                              header = TRUE, row.names = 1)
#colnames(RNAseq.GTEx.raw) <- gsub("\\.","-",gsub("(GTEX\\.\\w{4,5})\\..*","\\1",colnames(RNAseq.GTEx.raw)))
RNAseq.GTEx.raw$gene_id <- NULL
# Remove X-Y genes
RNAseq.GTEx.raw <- RNAseq.GTEx.raw[-grep("PAR",rownames(RNAseq.GTEx.raw)),]
rownames(RNAseq.GTEx.raw) <- gsub("(ENST.*)\\..*","\\1",rownames(RNAseq.GTEx.raw))
tissue.samples <- colnames(RNAseq.GTEx.raw)
# Round values
RNAseq.GTEx.raw <- round(RNAseq.GTEx.raw)
# 1.2) RNAseq TPM
RNAseq.GTEx.TPM <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
                              header = TRUE, row.names = 1)
RNAseq.GTEx.TPM$gene_id <- NULL
# Remove X-Y genes
RNAseq.GTEx.TPM <- RNAseq.GTEx.TPM[-grep("PAR",rownames(RNAseq.GTEx.TPM)),]
rownames(RNAseq.GTEx.TPM) <- gsub("(ENST.*)\\..*","\\1",rownames(RNAseq.GTEx.TPM))
tissue.samples <- colnames(RNAseq.GTEx.TPM)
# Round values
RNAseq.GTEx.TPM <- round(RNAseq.GTEx.TPM)
# 1.3) Protein coding transcripts
ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t")

# 2) Add metadata information for all samples from that cancer
GTEx.samples.metadata <- data.frame(sample = tissue.samples)
# 2.1) Sample library size
# Estimate normalization factor for each sample
library.size.norm.factors <- estimateSizeFactorsForMatrix(RNAseq.GTEx.raw)
library.size.norm.factors <- data.frame(sample_lib_size = library.size.norm.factors)
GTEx.samples.metadata <- merge(GTEx.samples.metadata,library.size.norm.factors, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 2.2) PCA of raw VST-normalized RNAseq matrix to take PC weigths for transcripts
# 2.2.1) Filter low-expressed genes
RNAseq.GTEx.raw.filt <- RNAseq.GTEx.raw[rowSums(log2(RNAseq.GTEx.raw) >= 1) >= round(length(colnames(RNAseq.GTEx.raw)) * 0.50),]
# 2.2.2) Filter out non-coding transcripts
RNAseq.GTEx.raw.filt <- RNAseq.GTEx.raw.filt[rownames(RNAseq.GTEx.raw.filt)%in%ensembl.v88.coding.transcripts$V1,]
# 2.2.3) PCA
PC.RNAseq.names.tissue <- PCA_type_RNAseq(RNAseq.matrix = RNAseq.GTEx.raw.filt, PCs = 4, scale = TRUE, center = TRUE)

# 2.3) Coefficient of variation (SD/mean) per transcript
coeff.var.transcripts <- apply(RNAseq.GTEx.TPM,1, function(transcript) {
    mean.transcript <- mean(na.omit(as.numeric(transcript)))
    SD.transcript <- sd(transcript,na.rm = TRUE)
    coeff.var <- SD.transcript/mean.transcript
    coeff.var 
})
coeff.var.transcripts.df <- data.frame(coeff_var = coeff.var.transcripts)

# 2.4) Median transcript expression
# 2.4.1) TPM
median.TPM.expression.transcripts <- apply(RNAseq.GTEx.TPM,1, function(transcript) {
    median.exp.transcript <- median(na.omit(as.numeric(transcript)))
    median.exp.transcript
})
median.TPM.expression.transcripts.df <- data.frame(median_TPM_exp_transcript = median.TPM.expression.transcripts)
# 2.4.2) raw
median.raw.expression.transcripts <- apply(RNAseq.GTEx.raw,1, function(transcript) {
    median.exp.transcript <- median(na.omit(as.numeric(transcript)))
    median.exp.transcript
})
median.raw.expression.transcripts.df <- data.frame(median_raw_exp_transcript = median.raw.expression.transcripts)
# Merge
median.expression.transcripts <- merge(median.TPM.expression.transcripts.df,median.raw.expression.transcripts.df, by = "row.names", all.x = TRUE)
transcripts.df <- merge(median.expression.transcripts,coeff.var.transcripts.df, by.x = "Row.names", by.y ="row.names", all.x = TRUE)
rownames(transcripts.df) <- transcripts.df$Row.names
transcripts.df$Row.names <- NULL

# 3) Save results
GTEx.metadata.path <- gsub("\\[X1\\]",GTEx.tissue, GTEx.metadata.path)

# 3.1) Samples Metadata
write.table(GTEx.samples.metadata, file = gsub("\\[X2\\]",GTEx.tissue, paste0(GTEx.metadata.path,paths[paths$folder_or_object=="GTEx_samples_metadata","path_or_filename"])), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# 3.2) Transcripts coefficient variation 
write.table(transcripts.df, file = gsub("\\[X2\\]",GTEx.tissue, paste0(GTEx.metadata.path,paths[paths$folder_or_object=="GTEx_transcripts_exp","path_or_filename"])), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# 3.3) PCA results
# Tissue
for (view in c("var","ind")) {
    # RNAseq
    PCA.path <- paste0(GTEx.RNAseq.PCA.path,GTEx.tissue,"_tissue_PCA_",view,".txt")
    PC.RNAseq.names.tissue.filt <- PC.RNAseq.names.tissue[grep(view, PC.RNAseq.names.tissue)]
    write.table(eval(parse(text = PC.RNAseq.names.tissue.filt)), file = PCA.path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
}









