rm(list=ls())


################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################


PCA_type_RNAseq <- function(RNAseq.matrix, type, cancer.subtypes = NULL, PCs = 4, scale = FALSE, center = FALSE) {

  PCs.names.var <- c()
  PCs.names.ind <- c()
  if (type == "tissue") {
    print("... PCA for whole tissue ...")
    # VST transformation
    RNAseq.matrix.raw.vst <- vst(object = as.matrix(RNAseq.matrix), nsub = 10000)
    RNAseq.matrix.raw.vst.tp <- t(RNAseq.matrix.raw.vst)
    # Scale and center
    #RNAseq.matrix.raw.vst.tp <- scale(RNAseq.matrix.raw.vst.tp, center = center, scale = scale)
    # Add group
    RNAseq.matrix.raw.vst.tp.group <- merge(RNAseq.matrix.raw.vst.tp,cancer.subtypes[,"cluster", drop = FALSE], by.x = "row.names", by.y = "row.names", all.x = TRUE)
    rownames(RNAseq.matrix.raw.vst.tp.group) <- RNAseq.matrix.raw.vst.tp.group$Row.names
    RNAseq.matrix.raw.vst.tp.group$Row.names <- NULL
    # PCA
    res.pca <- PCA(RNAseq.matrix.raw.vst.tp.group[,-length(colnames(RNAseq.matrix.raw.vst.tp.group))], scale.unit = TRUE, ncp = 1000, graph = FALSE)
    # Plot
    PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_tissue_PCA_ind.png")
    png(PCA.path, width = 4500, height = 3000, res = 300)
    p <- fviz_pca_ind(res.pca, col.ind = "cos2",
                label = "none", # hide individual labels,
                habillage = as.factor(RNAseq.matrix.raw.vst.tp.group$cluster)
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
  } else if (type == "subtissue") {
    print("... PCA for each subtype ...")
    for (i in 1:max(cancer.subtypes$cluster)) {
        # Filter RNAseq matrix for each cluster or subtype
        cluster <- cancer.subtypes[cancer.subtypes$cluster == i,]
        samples <- substr(colnames(RNAseq.matrix),1,12)%in%rownames(cluster)
        print(paste0("subtype: ",i, " with ", nrow(cluster), " samples"))
        if ( (nrow(cluster) < 5) || (sum(samples) < 5) ) {
            print("Less than 5 samples, skipping...")
            next
        }
        RNAseq.matrix.raw.cluster <- RNAseq.matrix[,samples]
        # VST transformation
        RNAseq.matrix.raw.cluster.vst <- vst(object = as.matrix(RNAseq.matrix.raw.cluster), nsub = 10000)
        RNAseq.matrix.raw.cluster.vst.tp <- t(RNAseq.matrix.raw.cluster.vst)
        # Scale and center
        #RNAseq.matrix.raw.cluster.vst.tp <- scale(RNAseq.matrix.raw.cluster.vst.tp, center = center, scale = scale)
        # PCA
        res.pca <- PCA(RNAseq.matrix.raw.cluster.vst.tp, scale.unit = TRUE, ncp = 1000, graph = FALSE)
        PCs.ob <- ncol(res.pca$ind$coord)
        if (PCs.ob == 1) {
            PCs.to.keep <- 1
        } else if (PCs < PCs.ob) {
            PCs.to.keep <- PCs
        } else {
            PCs.to.keep <- PCs.ob
        }
        # Plot
        PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_subtissue_",i,"_PCA_ind.png")
        png(PCA.path, width = 4500, height = 3000, res = 300)
        p <- fviz_pca_ind(res.pca, col.ind = "cos2",
            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
            label = "none" # hide individual labels
            )
        print(p)
        dev.off()
        # Extract the results for variables (genes)
        res.pca.var <- get_pca_var(res.pca)
        # Extract the results for individuals
        res.pca.ind <- get_pca_ind(res.pca)
        # Create a new variable for each PC to be used later on the regression
        print(paste0("PCs --> ", PCs.to.keep))
        # Variables
        PC.var <- paste0("PCA.RNAseq.var.",i)
        assign(PC.var, data.frame(res.pca.var$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
        PCs.names.var <- c(PCs.names.var, PC.var)
        # Individuals
        PC.ind <- paste0("PCA.RNAseq.ind.",i)
        assign(PC.ind, data.frame(res.pca.ind$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
        PCs.names.ind <- c(PCs.names.ind, PC.ind)
    }
  }
  return(c(PCs.names.var,PCs.names.ind))
}

PCA_type_CNV <- function(CNV.file, type, cancer.subtypes = NULL, PCs = 4, scale = FALSE, center = FALSE) {

  PCs.names.var <- c()
  PCs.names.ind <- c()
  if (type == "tissue") {
    print("... PCA somatic CNV for whole tissue ...")
    # PCA
    CNV.file.t <- t(CNV.file)
    CNV.file.t <- scale(CNV.file.t, center = center, scale = scale)
    res.pca <- PCA(CNV.file.t, scale.unit = FALSE, ncp = 1000, graph = FALSE)
    # Plot
    PCA.path <- paste0(TCGA.CNV.PCA.path,TCGA.cancer,"_tissue_PCA_ind.png")
    png(PCA.path, width = 4500, height = 3000, res = 300)
    p <- fviz_pca_ind(res.pca, col.ind = "cos2",
                axes = c(1, 2),
                label = "none")
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
    PC.var <- "PCA.CNV.var"
    assign(PC.var, data.frame(scale(res.pca.var$contrib[,1:PCs.to.keep])), envir = .GlobalEnv)
    PCs.names.var <- c(PCs.names.var, PC.var)
    # Individuals
    PC.ind <- "PCA.CNV.ind"
    assign(PC.ind, data.frame(scale(res.pca.ind$contrib[,1:PCs.to.keep])), envir = .GlobalEnv)
    PCs.names.ind <- c(PCs.names.ind, PC.ind)
  } else if (type == "subtissue") {
    print("... PCA CNV for each subtype ...")
    for (i in 1:max(cancer.subtypes$cluster)) {
        # Filter RNAseq matrix for each cluster or subtype
        cluster <- cancer.subtypes[cancer.subtypes$cluster == i,]
        samples <- substr(colnames(CNV.file),1,12)%in%rownames(cluster)
        print(paste0("subtype: ",i, " with ", nrow(cluster), " samples"))
        if ( (nrow(cluster) < 5) || (sum(samples) < 5) ) {
            print("Less than 5 samples, skipping...")
            next
        }
        CNV.file.subtype <- CNV.file[,samples]
        # PCA
        CNV.file.subtype <- t(CNV.file.subtype)
        CNV.file.subtype <- scale(CNV.file.subtype, center = center, scale = scale)
        res.pca <- PCA(CNV.file.subtype, scale.unit = FALSE, ncp = 1000, graph = FALSE)
        PCs.ob <- ncol(res.pca$ind$coord)
        if (PCs.ob == 1) {
            PCs.to.keep <- 1
        } else if (PCs < PCs.ob) {
            PCs.to.keep <- PCs
        } else {
            PCs.to.keep <- PCs.ob
        }
        # Plot
        PCA.path <- paste0(TCGA.CNV.PCA.path,TCGA.cancer,"_subtissue_",i,"_PCA_ind.png")
        png(PCA.path, width = 4500, height = 3000, res = 300)
        p <- fviz_pca_ind(res.pca, col.ind = "cos2",
            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
            label = "none" # hide individual labels
            )
        print(p)
        dev.off()
        # Extract the results for variables (genes)
        res.pca.var <- get_pca_var(res.pca)
        # Extract the results for individuals
        res.pca.ind <- get_pca_ind(res.pca)
        # Create a new variable for each PC to be used later on the regression
        print(paste0("PCs --> ", PCs.to.keep))
        # Variables
        PC.var <- paste0("PCA.CNV.var.",i)
        assign(PC.var, data.frame(res.pca.var$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
        PCs.names.var <- c(PCs.names.var, PC.var)
        # Individuals
        PC.ind <- paste0("PCA.CNV.ind.",i)
        assign(PC.ind, data.frame(res.pca.ind$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
        PCs.names.ind <- c(PCs.names.ind, PC.ind)
    }
  }
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
TCGA.cancer <- args[2]

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/collabs/Ignasi/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
CNV.path <- paths[paths$folder_or_object=="CNV_file_path","path_or_filename"]
firehose.subtypes.path <- paths[paths$folder_or_object=="firehose_subtypes_path","path_or_filename"]
TCGA.RNAseq.PCA.path <- paths[paths$folder_or_object=="TCGA_RNAseq_PCA_path","path_or_filename"]
TCGA.CNV.PCA.path <- paths[paths$folder_or_object=="TCGA_CNV_PCA_path","path_or_filename"]
TCGA.TMB.path <- paths[paths$folder_or_object=="TCGA_TMB_path","path_or_filename"]
TCGA.metadata.path <- paths[paths$folder_or_object=="TCGA_metadata_path","path_or_filename"]

print(paste0("TCGA cancer --> ",TCGA.cancer))

# 1.1) RNAseq matrix raw counts
RNAseq.matrix.raw <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])),
                                              header = TRUE, sep = "\t", row.names = 1)
RNAseq.matrix.raw <- RNAseq.matrix.raw[order(rownames(RNAseq.matrix.raw)),]
cancer.samples <- substr(colnames(RNAseq.matrix.raw),1,12)
# Filter for PRAD
RNAseq.matrix.raw <- t(na.omit(t(RNAseq.matrix.raw)))
# Round values
RNAseq.matrix.raw <- round(RNAseq.matrix.raw)

# 1.2) RNAseq TPM
RNAseq.TCGA.TPM <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
                              header = TRUE, sep = "\t", row.names = 1)
cancer.samples <- substr(colnames(RNAseq.TCGA.TPM),1,12)
# 1.3) Firehose subtypes (mRNA NMF) and CNV
if (TCGA.cancer == "TCGA-SKCM") {
    cancer.subtypes.path <- gsub("-TP","-TM",paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"]))
    CNV.file.path <- gsub("-TP","-TM",paste0(CNV.path,paths[paths$folder_or_object=="CNV_file","path_or_filename"]))
} else if (TCGA.cancer == "TCGA-LAML")  {
    cancer.subtypes.path <- gsub("-TP","-TB",paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"]))
    CNV.file.path <- gsub("-TP","-TB",paste0(CNV.path,paths[paths$folder_or_object=="CNV_file","path_or_filename"]))
} else {
    cancer.subtypes.path <- paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"])
    CNV.file.path <- paste0(CNV.path,paths[paths$folder_or_object=="CNV_file","path_or_filename"])
}
# Subtypes
cancer.subtypes <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), cancer.subtypes.path), 
                                           header = TRUE, sep = "\t", skip = 1)
cancer.subtypes$SampleName<- gsub("-",".", cancer.subtypes$SampleName, fixed = TRUE)
rownames(cancer.subtypes) <- substr(cancer.subtypes$SampleName,1,12)
# CNV
CNV.file <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), CNV.file.path ), 
                                           header = TRUE, sep = "\t")
rownames(CNV.file) <- CNV.file$Gene.Symbol
CNV.file <- CNV.file[,grep("TCGA*",colnames(CNV.file))]
colnames(CNV.file) <- substr(colnames(CNV.file),1,12)
# 1.4) Leukocyte Fraction
leukocyte.fraction <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="leuk_fraction","path_or_filename"]),
                                 header = FALSE, sep = "\t")
colnames(leukocyte.fraction) <- c("cancer","sample","LF")
leukocyte.fraction$sample <- gsub("-",".",substr(leukocyte.fraction$sample,1,12))
# 1.5) Tumor Purity
tumor.purity <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="tumor_purity","path_or_filename"]),
                           header = TRUE, sep = "\t")
tumor.purity$sample <- gsub("-",".",substr(tumor.purity$sample,1,12))
# Merge and filter by cancer samples
LF.purity <- merge(leukocyte.fraction, tumor.purity, by.x = "sample", by.y = "sample", all.x = TRUE)
LF.purity <- LF.purity[LF.purity$sample%in%cancer.samples,c("sample","cancer","LF","purity")]
# 1.6) Protein coding transcripts
ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t")
# 1.7) Tumor Mutation Burden
TMB_TCGA <- read.csv(file = paste0(TCGA.TMB.path,paths[paths$folder_or_object=="TCGA_TMB","path_or_filename"]), row.names = 1)
# Some samples are duplicated
TMB_TCGA <- TMB_TCGA[-which(duplicated(substr(rownames(TMB_TCGA),1,12))),]
rownames(TMB_TCGA) <- substr(rownames(TMB_TCGA),1,12)

# 2) Add metadata information for all samples from that cancer
TCGA.samples.metadata <- data.frame(sample = cancer.samples)
# 2.1) Leukocyte Fraction and purity
# Values
TCGA.samples.metadata <- merge(TCGA.samples.metadata,LF.purity[,c("sample","LF","purity")], by.x = "sample", by.y = "sample", all.x = TRUE)
# Lowest (LF) or first (purity) decile to remove?
LF.deciles <- quantile(na.omit(LF.purity$LF), prob=seq(0, 1, length = 11))
purity.deciles <- quantile(na.omit(LF.purity$purity), prob=seq(0, 1, length = 11))
LF.samples.to.remove <- LF.purity[which(LF.purity$LF >=LF.deciles[10]),"sample"]
purity.samples.to.remove <- LF.purity[which(LF.purity$purity <= purity.deciles[2]),"sample"]
TCGA.samples.metadata$LF_remove <- "no"
TCGA.samples.metadata[TCGA.samples.metadata$sample %in% LF.samples.to.remove,"LF_remove"] <- "yes"
TCGA.samples.metadata$purity_remove <- "no"
TCGA.samples.metadata[TCGA.samples.metadata$sample %in% purity.samples.to.remove,"purity_remove"] <- "yes"
# 2.2) CNA file available and counts
# samples not found in CNV file
CNV.missing.samples <- colnames(RNAseq.TCGA.TPM)[!colnames(RNAseq.TCGA.TPM) %in% colnames(CNV.file)]
TCGA.samples.metadata[,"CNV_available"] <- "yes"
TCGA.samples.metadata[TCGA.samples.metadata$sample %in% CNV.missing.samples,"CNV_available"] <- "no"
# CNV burden
CNV.burden <- data.frame(CNV_burden = colSums(abs(CNV.file)))
TCGA.samples.metadata <- merge(TCGA.samples.metadata,CNV.burden, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 2.3) PCA of somatic CNA to take PC weights for samples and genes
# Split CNV variable in two: amp and del for each gene
TCGA_pancancer_CNV_split <- apply(CNV.file,2,function(gene) {
  amp <- ifelse(gene < 0, 0, as.numeric(gene))
  #rownames(amp) <- paste0(rownames(amp),"_amp")
  del <- ifelse(gene > 0, 0, as.numeric(gene))
  #rownames(del) <- paste0(rownames(del),"_amp")
  gene_amp_del <- rbind(amp,del)
  gene_amp_del
})
TCGA_pancancer_CNV_split <- as.data.frame(TCGA_pancancer_CNV_split)
amp_index <- seq(from = 1, to = nrow(TCGA_pancancer_CNV_split),by = 2)
del_index <- seq(from = 2, to = nrow(TCGA_pancancer_CNV_split),by = 2)
rownames(TCGA_pancancer_CNV_split)[amp_index] <- paste0(rownames(CNV.file),"_amp")
rownames(TCGA_pancancer_CNV_split)[del_index] <- paste0(rownames(CNV.file),"_del")
PC.CNV.names.tissue <- PCA_type_CNV(CNV.file = TCGA_pancancer_CNV_split, type = "tissue", cancer.subtypes = cancer.subtypes, PCs = 1000, scale = FALSE, center = FALSE)
PC.CNV.names.subtissue <- PCA_type_CNV(CNV.file = TCGA_pancancer_CNV_split, type = "subtissue", cancer.subtypes = cancer.subtypes, PCs = 1000, scale = FALSE, center = FALSE)
# 2.4) RNAseq NMF subtype
TCGA.samples.metadata <- merge(TCGA.samples.metadata,cancer.subtypes[,"cluster", drop = FALSE], by.x = "sample", by.y = "row.names", all.x = TRUE)
colnames(TCGA.samples.metadata)[colnames(TCGA.samples.metadata) %in% "cluster"] <- "RNAseq_NMF_cluster"
# 2.5) Sample library size
# Estimate normalization factor for each sample
library.size.norm.factors <- estimateSizeFactorsForMatrix(RNAseq.matrix.raw)
library.size.norm.factors <- data.frame(sample_lib_size = library.size.norm.factors)
TCGA.samples.metadata <- merge(TCGA.samples.metadata,library.size.norm.factors, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 2.6) PCA of raw VST-normalized RNAseq matrix to take PC weigths for transcripts
# 2.6.1) Filter low-expressed genes
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw[rowSums(log2(RNAseq.matrix.raw) >= 1) >= round(length(colnames(RNAseq.matrix.raw)) * 0.50),]
# 2.6.2) Filter out non-coding transcripts
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw.filt[rownames(RNAseq.matrix.raw.filt)%in%ensembl.v88.coding.transcripts$V1,]
# 2.6.3) PCA
PC.RNAseq.names.tissue <- PCA_type_RNAseq(RNAseq.matrix = RNAseq.matrix.raw.filt, type = "tissue", cancer.subtypes = cancer.subtypes, PCs = 4, scale = TRUE, center = TRUE)
PC.RNAseq.names.subtissue <- PCA_type_RNAseq(RNAseq.matrix = RNAseq.matrix.raw.filt, type = "subtissue", cancer.subtypes = cancer.subtypes, PCs = 4, scale = TRUE, center = TRUE)
# 2.7) Coefficient of variation (SD/mean) per transcript
coeff.var.transcripts <- apply(RNAseq.TCGA.TPM,1, function(transcript) {
    mean.transcript <- mean(na.omit(as.numeric(transcript)))
    SD.transcript <- sd(transcript,na.rm = TRUE)
    coeff.var <- SD.transcript/mean.transcript
    coeff.var 
})
coeff.var.transcripts.df <- data.frame(coeff_var = coeff.var.transcripts)
# 2.8) Tumor Mutation Burden
TCGA.samples.metadata <- merge(TCGA.samples.metadata,TMB_TCGA, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 2.9) Median transcript expression
# 2.9.1) TPM
median.TPM.expression.transcripts <- apply(RNAseq.TCGA.TPM,1, function(transcript) {
    median.exp.transcript <- median(na.omit(as.numeric(transcript)))
    median.exp.transcript
})
median.TPM.expression.transcripts.df <- data.frame(median_TPM_exp_transcript = median.TPM.expression.transcripts)
# 2.9.2) raw
median.raw.expression.transcripts <- apply(RNAseq.matrix.raw,1, function(transcript) {
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
# 3.1) Samples Metadata
write.table(TCGA.samples.metadata, file = gsub("\\[X\\]",TCGA.cancer, paste0(TCGA.metadata.path,paths[paths$folder_or_object=="TCGA_samples_metadata","path_or_filename"])), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# 3.2) Transcripts expression + coefficient variation 
write.table(transcripts.df, file = gsub("\\[X\\]",TCGA.cancer, paste0(TCGA.metadata.path,paths[paths$folder_or_object=="transcripts_exp","path_or_filename"])), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# 3.3) PCA results
# Tissue
for (view in c("var","ind")) {
    # RNAseq
    PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_tissue_PCA_",view,".txt")
    PC.RNAseq.names.tissue.filt <- PC.RNAseq.names.tissue[grep(view, PC.RNAseq.names.tissue)]
    write.table(eval(parse(text = PC.RNAseq.names.tissue.filt)), file = PCA.path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
    # CNV
    PCA.path <- paste0(TCGA.CNV.PCA.path,TCGA.cancer,"_tissue_PCA_",view,".txt")
    PC.CNV.names.tissue.filt <- PC.CNV.names.tissue[grep(view, PC.CNV.names.tissue)]
    write.table(eval(parse(text = PC.CNV.names.tissue.filt)), file = PCA.path , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
}
# Subtissues
for (i in 1:max(cancer.subtypes$cluster)) {
    for (view in c("var","ind")) {
        # RNAseq
        PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_subtissue_",i,"_PCA_",view,".txt")
        PC.RNAseq.names.tissue.filt <- PC.RNAseq.names.subtissue[grep(paste0(view,".",i), PC.RNAseq.names.subtissue)]
        write.table(eval(parse(text = PC.RNAseq.names.tissue.filt)), file = PCA.path , sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = TRUE)
        # CNV
        PCA.path <- paste0(TCGA.CNV.PCA.path,TCGA.cancer,"_subtissue_",i,"_PCA_",view,".txt")
        PC.CNV.names.tissue.filt <- PC.CNV.names.subtissue[grep(paste0(view,".",i), PC.CNV.names.subtissue)]
        write.table(eval(parse(text = PC.CNV.names.tissue.filt)), file = PCA.path , sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = TRUE)
    }
}









