PCA_CNV <- function(CNV.file, PCs = 4, scale = FALSE, center = FALSE) {
  PCs.names.var <- c()
  PCs.names.ind <- c()
  print("... PCA somatic CNV ...")
  # PCA
  CNV.file.t <- t(CNV.file)
  CNV.file.t <- scale(CNV.file.t, center = center, scale = scale)
  res.pca <- PCA(CNV.file.t, scale.unit = FALSE, ncp = PCs, graph = FALSE)
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
  assign(PC.var, data.frame(res.pca.var$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
  PCs.names.var <- c(PCs.names.var, PC.var)
  # Individuals
  PC.ind <- "PCA.CNV.ind"
  assign(PC.ind, data.frame(res.pca.ind$contrib[,1:PCs.to.keep]), envir = .GlobalEnv)
  PCs.names.ind <- c(PCs.names.ind, PC.ind)
  # Output
  return(c(PCs.names.var,PCs.names.ind))
}

sparse_PCA_CNV <- function(CNV.file, PCs = 10, alpha = 0.0001, scale = FALSE, center = FALSE, max.iter, robust) {
  PCs.names.var <- c()
  PCs.names.ind <- c()
  print("... PCA somatic CNV ...")
  # sparse PCA
  CNV.file.t <- t(CNV.file)
  CNV.file.t <- scale(CNV.file.t, center = center, scale = scale)
  if (robust == "yes") {
    # This is not robust
    sparsePCA <- rspca(X = CNV.file.t, k = PCs,  alpha = alpha, center = center, scale = scale, max_iter = max.iter)
  } else {
    sparsePCA <- spca(X = CNV.file.t, k = PCs,  alpha = alpha, center = center, scale = scale, max_iter = max.iter)
  }
  PCs.ob <- ncol(sparsePCA$scores)
  if (PCs.ob == 1) {
      PCs.to.keep <- 1
  } else if (PCs < PCs.ob) {
      PCs.to.keep <- PCs
  } else {
      PCs.to.keep <- PCs.ob
  }
  # Extract the results for variables (genes)
  res.pca.var <- data.frame(sparsePCA$loadings)
  rownames(res.pca.var) <- colnames(CNV.file.t)
  colnames(res.pca.var) <- gsub("X","Dim\\.",colnames(res.pca.var))
  # Extract the results for individuals
  res.pca.ind <- data.frame(sparsePCA$scores)
  colnames(res.pca.ind) <- gsub("X","Dim\\.",colnames(res.pca.ind))
  print(paste0("PCs --> ", PCs.to.keep))
  # Variables
  PC.var <- "PCA.CNV.var"
  assign(PC.var, res.pca.var[,1:PCs.to.keep], envir = .GlobalEnv)
  PCs.names.var <- c(PCs.names.var, PC.var)
  # Individuals
  PC.ind <- "PCA.CNV.ind"
  assign(PC.ind, res.pca.ind[,1:PCs.to.keep], envir = .GlobalEnv)
  PCs.names.ind <- c(PCs.names.ind, PC.ind)
  # Output
  return(c(PCs.names.var,PCs.names.ind))
}

PCA_RNAseq <- function(RNAseq.matrix) {

  PCs.names.var <- c()
  PCs.names.ind <- c()

    print("... PCA RNAseq for pancancer ...")
    # VST transformation
    RNAseq.matrix.raw.vst <- vst(object = as.matrix(RNAseq.matrix), nsub = 10000)
    RNAseq.matrix.raw.vst.tp <- t(RNAseq.matrix.raw.vst)
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
    if (ncol(res.pca$ind$coord) == 1) {
        PCs.to.keep <- 1
    } else {PCs.to.keep <- 4}
    # Extract the results for variables (genes)
    res.pca.var <- get_pca_var(res.pca)
    # Extract the results for individuals
    res.pca.ind <- get_pca_ind(res.pca)
    print(paste0("PCs --> ", PCs.to.keep))
    # Variables
    PC.var <- "PCA.RNAseq.var"
    assign(PC.var, data.frame(res.pca.var$contrib[,1:4]), envir = .GlobalEnv)
    PCs.names.var <- c(PCs.names.var, PC.var)
    # Individuals
    PC.ind <- "PCA.RNAseq.ind"
    assign(PC.ind, data.frame(res.pca.ind$contrib[,1:4]), envir = .GlobalEnv)
    PCs.names.ind <- c(PCs.names.ind, PC.ind)

  return(c(PCs.names.var,PCs.names.ind))
}

# PCA
.libPaths( rev( .libPaths() ) )
library("factoextra")
library("FactoMineR")
library("sparsepca")

args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
alpha <- as.numeric(args[2])
robust_SPCA <- args[3]
num_PCs <- as.numeric(args[4])

print(args)

# 1) Open files

# 1.1) CNV files
TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1

TCGA_CNV <- c()
print(type)
if (type == "pancancer") {
  for (TCGA_cancer in TCGA_cancers) {
    print(TCGA_cancer)
    tryCatch({
      CNV_file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_",gsub("TCGA-","",TCGA_cancer),"-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt")
      if (TCGA_cancer == "TCGA-SKCM") {
        CNV_file_path <- gsub("-TP","-TM",paste0(CNV_file_path))
      }
      # CNV
      CNV_file <- read.table(file = CNV_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      colnames(CNV_file) <- substr(colnames(CNV_file),1,12)
      rownames(CNV_file) <- CNV_file$Gene.Symbol
      CNV_file <- CNV_file[,grep("TCGA*",colnames(CNV_file))]
      if (length(TCGA_CNV) == 0) {
        TCGA_CNV <- CNV_file
      } else {
        TCGA_CNV <- merge(TCGA_CNV, CNV_file, by.x= "row.names", by.y = "row.names", all.x = TRUE)
        rownames(TCGA_CNV) <- TCGA_CNV$Row.names
        TCGA_CNV$Row.names <- NULL
      }
    }, error = function(e){
        print(e)
        }
    )
  }
} else {
  TCGA_cancer <- type
  tryCatch({
    CNV_file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_",gsub("TCGA-","",TCGA_cancer),"-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt")
    if (TCGA_cancer == "TCGA-SKCM") {
      CNV_file_path <- gsub("-TP","-TM",paste0(CNV_file_path))
    } else if (TCGA_cancer == "TCGA-LAML") {
      CNV_file_path <- gsub("-TP","-TB",paste0(CNV_file_path))
    }
    # CNV
    CNV_file <- read.table(file = CNV_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(CNV_file) <- substr(colnames(CNV_file),1,12)
    rownames(CNV_file) <- CNV_file$Gene.Symbol
    CNV_file <- CNV_file[,grep("TCGA*",colnames(CNV_file))]
    if (length(TCGA_CNV) == 0) {
      TCGA_CNV <- CNV_file
    } else {
      TCGA_CNV <- merge(TCGA_CNV, CNV_file, by.x= "row.names", by.y = "row.names", all.x = TRUE)
      rownames(TCGA_CNV) <- TCGA_CNV$Row.names
      TCGA_CNV$Row.names <- NULL
    }
  }, error = function(e){
      print(e)
      }
  )
}
print(dim(TCGA_CNV))

# 1.2) RNAseq files

# 2) PCAs

# 2.1) PCA CNV

# Split CNV variable in two: amp and del for each gene
TCGA_CNV_split <- apply(TCGA_CNV,2,function(gene) {
  amp <- ifelse(gene < 0, 0, as.numeric(gene))
  #rownames(amp) <- paste0(rownames(amp),"_amp")
  del <- ifelse(gene > 0, 0, as.numeric(gene))
  #rownames(del) <- paste0(rownames(del),"_amp")
  gene_amp_del <- rbind(amp,del)
  gene_amp_del
})
TCGA_CNV_split <- as.data.frame(TCGA_CNV_split)
amp_index <- seq(from = 1, to = nrow(TCGA_CNV_split),by = 2)
del_index <- seq(from = 2, to = nrow(TCGA_CNV_split),by = 2)
rownames(TCGA_CNV_split)[amp_index] <- paste0(rownames(TCGA_CNV),"_amp")
rownames(TCGA_CNV_split)[del_index] <- paste0(rownames(TCGA_CNV),"_del")

# Perform PCA
if ( type == "pancancer") {
  nPCs <- num_PCs
} else {
  nPCs <- ncol(TCGA_CNV_split)
}
print(nPCs)
#PC_CNV_names <- PCA_CNV(CNV.file = TCGA_CNV_split, PCs = nPCs, scale = FALSE, center = FALSE)

# 2.2) Sparse PCA CNV
sparse_PC_CNV_names <- sparse_PCA_CNV(CNV.file = TCGA_CNV_split, alpha = alpha, PCs = nPCs, scale = FALSE, center = FALSE, max.iter = 200, robust = robust_SPCA)

#for (type in c("PCA","sparse_PCA")) {
for (PCA_type in c("sparse_PCA")) {
  for (view in c("var","ind")) {
    if (PCA_type == "PCA") {
      PC_CNV_names_filt <- PC_CNV_names[grep(view, PC_CNV_names)]
    } else if (PCA_type == "sparse_PCA") {
      PC_CNV_names_filt <- sparse_PC_CNV_names[grep(view, sparse_PC_CNV_names)]
    }
    PCA_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",type,"_",PCA_type,"_",view,"_",alpha,"_robust_",robust_SPCA,"_num_PCs_",num_PCs,".txt")
    print(paste0("Writting results to --> ",PCA_path))
    write.table(eval(parse(text = PC_CNV_names_filt)), file = PCA_path , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
  }
}


