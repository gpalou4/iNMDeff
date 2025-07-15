rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

somatic_transcripts_filtering <- function(RNAseq.matrix, transcripts.filtering.matrix, TCGA.barcode) {
  sample.transcripts.to.remove <- rownames(transcripts.filtering.matrix[transcripts.filtering.matrix[,TCGA.barcode]==1,])
  RNAseq.matrix <- RNAseq.matrix[!rownames(RNAseq.matrix)%in%sample.transcripts.to.remove,, drop = FALSE]
  if(is.vector(RNAseq.matrix)){return(NULL)}
  # Remove genes with 0 expression (sample-specific...)
  #RNAseq.matrix <- RNAseq.matrix[RNAseq.matrix[,sample] != 0,]
  return(RNAseq.matrix)
}

germline_transcripts_filtering <- function(RNAseq.matrix, transcripts.filtering.matrix, TCGA.barcode) {
  germ.mut <- "stopgain|[^non]frameshift deletion|[^non]frameshift insertion|splicing|nonframeshift deletion|nonframeshift insertion|startloss|stoploss"
  transcripts.germ.mut.sample <- transcripts.filtering.matrix[,TCGA.barcode, drop = FALSE]
  sample.transcripts.to.remove <- rownames(transcripts.germ.mut.sample[grep(germ.mut,transcripts.germ.mut.sample[,1]),, drop = FALSE])
  RNAseq.matrix <- RNAseq.matrix[!rownames(RNAseq.matrix)%in%sample.transcripts.to.remove,, drop = FALSE]
  if(is.vector(RNAseq.matrix)){return(NULL)}
  return(RNAseq.matrix)
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

# Library home

#install.packages("",lib="/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6")
#.libPaths("/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6")
.libPaths( rev( .libPaths() ) )
# Bayesian for NB regression
library("V8")
# library("shinystan")
library("rstan")
library("rstanarm")
#.libPaths( c( "/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6" , .libPaths() ) )
# Library Size and VST
library("DBI")
library("DESeq2")
# PCA
library("factoextra")
library("FactoMineR")
#library("RColorBrewer")
# Broken Stick PC
library("PCDimension")

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA/endogenous_NMD_efficiency_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#TCGA.cancer <- "TCGA-CHOL"
#nextflow.params <- c("somatic_transcripts_filt","germline_transcripts_filt")

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
TCGA.cancer <- args[2]
nextflow.params <- args[3:length(args)]

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_raw_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
firehose.subtypes.path <- paths[paths$folder_or_object=="firehose_subtypes_path","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
transcripts.to.remove.path <- paths[paths$folder_or_object=="transcripts_to_remove_path","path_or_filename"]
nb.coeff.res.path <- paths[paths$folder_or_object=="NB_results","path_or_filename"]
CNV.path <- paths[paths$folder_or_object=="CNV_file_path","path_or_filename"]
transcripts.germ.mut.path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
NMD.endogenous.transcripts.path <- paths[paths$folder_or_object=="NMD_endogenous_transcripts_path","path_or_filename"]
scripts.WD  <- "/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA/"

print(paste0("TCGA cancer --> ",TCGA.cancer))
print("nextflow parameters:")
print(nextflow.params)

# 1.1) RNAseq matrix raw counts
RNAseq.matrix.raw <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])),
                                              header = TRUE, sep = "\t", row.names = 1)
RNAseq.matrix.raw <- RNAseq.matrix.raw[order(rownames(RNAseq.matrix.raw)),]

cancer.samples <- substr(colnames(RNAseq.matrix.raw),1,12)
# Filter for PRAD
RNAseq.matrix.raw <- t(na.omit(t(RNAseq.matrix.raw)))
# Round values
RNAseq.matrix.raw <- round(RNAseq.matrix.raw)

# 1.1.2) RNAseq TPM
# RNAseq TPM
# RNAseq.TCGA.TPM <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
#                               header = TRUE, sep = "\t", row.names = 1)

# 1.2) Transcripts length
ensembl.v88.transcripts.length <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_transcripts_length","path_or_filename"]), 
                                header = TRUE, sep = "\t")
# 1.3) Protein coding transcripts
ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t",colClasses = "vector")             
# 1.4) transcripts with somatic mutations
transcripts.to.remove <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(transcripts.to.remove.path,paths[paths$folder_or_object=="transcripts_to_remove","path_or_filename"])),
                                    header = TRUE, sep = "\t", row.names = 1)
# 1.5) transcripts with germline mutations
transcripts.germ.mut <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(transcripts.germ.mut.path,paths[paths$folder_or_object=="transcripts_germ_mut","path_or_filename"])),
                                    header = TRUE, sep = "\t", row.names = 1)
transcripts.germ.mut[] <- lapply(transcripts.germ.mut , as.character)

# 1.5) Firehose subtypes (mRNA NMF) and CNV
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
CNV.file <- CNV.file[,grep("TCGA*",colnames(CNV.file))]
colnames(CNV.file) <- substr(colnames(CNV.file),1,12)

# 1.6) Leukocyte Fraction
leukocyte.fraction <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="leuk_fraction","path_or_filename"]),
                                 header = FALSE, sep = "\t")
colnames(leukocyte.fraction) <- c("cancer","sample","LF")
leukocyte.fraction$sample <- gsub("-",".",substr(leukocyte.fraction$sample,1,12))
# 1.7) Tumor Purity
tumor.purity <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="tumor_purity","path_or_filename"]),
                           header = TRUE, sep = "\t")
tumor.purity$sample <- gsub("-",".",substr(tumor.purity$sample,1,12))

# Merge and filter by cancer samples
LF.purity <- merge(leukocyte.fraction, tumor.purity, by.x = "sample", by.y = "sample", all.x = TRUE)
LF.purity <- LF.purity[LF.purity$sample%in%cancer.samples,c("sample","cancer","LF","purity")]

# 1.8) NMD targets genesets

# 1.8.1) NMD Colombo
NMD.Colombo <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_Colombo","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.2) NMD Karousis
NMD.Karousis <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_Karousis","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.3) NMD Tani
NMD.Tani <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_Tani","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.4) NMD Courtney
NMD.Courtney <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_Courtney","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.5) NMD ENSEMBL
NMD.ensembl <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_ensembl","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.6) NMD global
NMD.global <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.7) NMD global all shared
NMD.global.all.shared <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_all_shared","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.8) NMD global 2 shared
NMD.global.2.shared <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_2_shared","path_or_filename"]), 
                                         header = TRUE, sep = "\t")
# 1.8.9) SMG6
SMG6 <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="SMG6","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.10) SMG7
SMG7 <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="SMG7","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 1.8.11) Negative control (non NMD)
non.NMD.neg.control <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="non_NMD_genes","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
non.NMD.neg.control.with.NMD.features <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="non_NMD_genes_with_NMD_features","path_or_filename"]), 
                                        header = TRUE, sep = "\t")
# 2) Controls, offsets, PCA and filterings

# 2.1) Filter samples with highest Leukocyte Fraction (last decile)
LF.deciles <- quantile(na.omit(LF.purity$LF), prob=seq(0, 1, length = 11)) 
samples.to.remove <- LF.purity[which(LF.purity$LF >=LF.deciles[10]),"sample"]
# 2.2) Filter samples with lowest purity (first decile)
purity.deciles <- quantile(na.omit(LF.purity$purity), prob=seq(0, 1, length = 11)) 
samples.to.remove <- c(samples.to.remove,LF.purity[which(LF.purity$purity <= purity.deciles[2]),"sample"])
# 2.3) Filter samples with no CNA file
#colnames(transcripts.to.remove)[colSums(transcripts.to.remove) < 1000]
# samples not found in CNV file
CNV.samples.to.remove <- colnames(transcripts.to.remove)[!colnames(transcripts.to.remove)%in%colnames(CNV.file)]
# samples found in CNV file but with <10 CNV (weird)
#CNV.samples.to.remove <- c(CNV.samples.to.remove,colnames(CNV.file)[abs(colSums(CNV.file))<=10])
samples.to.remove <- c(samples.to.remove,CNV.samples.to.remove)
# Join
samples.to.remove <- unique(samples.to.remove)
#RNAseq.matrix.raw <- RNAseq.matrix.raw[,!substr(colnames(RNAseq.matrix.raw),1,12)%in%samples.to.remove]

# 2.2) Filter low-expressed genes
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw[rowSums(log2(RNAseq.matrix.raw) >= 1) >= round(length(colnames(RNAseq.matrix.raw)) * 0.50),]

# 2.3) Filter out non-coding transcripts
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw.filt[rownames(RNAseq.matrix.raw.filt)%in%ensembl.v88.coding.transcripts$V1,]

# 3) Negative Binomial regression on raw count data 
# Dataframe for coefficients for different NMD gene sets
NMD.gene.sets.names <- c("NMD.Colombo","NMD.Karousis","NMD.Tani","NMD.Courtney","NMD.ensembl","NMD.global",
                        "NMD.global.all.shared","NMD.global.2.shared","SMG6","SMG7",
                        "non.NMD.neg.control","non.NMD.neg.control.with.NMD.features")
NMD.gene.sets.names <- c("non.NMD.neg.control")
other.columns <- c("sample","subtype","purity","LF","num_NMD_targets","error")
df.cols <- length(other.columns)+(length(NMD.gene.sets.names)*5)

sbatch.commands <- data.frame(sbatch_command = NA)

total.num.NMDtargets.sample <- c()
for (sample in seq(1,ncol(RNAseq.matrix.raw))) {

  nb.coeff.res <- data.frame(matrix(ncol=df.cols,nrow=0, dimnames=list(NULL, c(other.columns,NMD.gene.sets.names,paste0(NMD.gene.sets.names,".pval"),
                                  paste0(NMD.gene.sets.names,".sd"),paste0(NMD.gene.sets.names,".CI_2.5"),paste0(NMD.gene.sets.names,".CI_97.5")))))
  start_time_total <- Sys.time()
  TCGA.barcode <- substr(colnames(RNAseq.matrix.raw)[sample],1,12)
  nb.coeff.res[1,"sample"] <- TCGA.barcode
  print(paste("SAMPLE", sample, " ----> ",TCGA.barcode),sep = "")
  sample.cluster <- cancer.subtypes[TCGA.barcode,"cluster"]
  purity <- LF.purity[LF.purity$sample%in%TCGA.barcode,"purity"]
  LF <- LF.purity[LF.purity$sample%in%TCGA.barcode,"LF"]
  
  # 4.1) Check LF and purity
  if (TCGA.barcode%in%samples.to.remove) {
    print("High Leukocyte Fraction / Low tumor purity sample, skipping...")
    nb.coeff.res[1,"error"] <- "High Leukocyte Fraction / Low tumor purity"
    #next
  } else if (is.na(sample.cluster)) { # Check if it sample contains mRNAseq cluster
    print("No sample mRNAseq cluster found, skipping...")
    nb.coeff.res[1,"error"] <- "No sample mRNAseq cluster"
    #next
  } else {
    # Save sample cluster
    nb.coeff.res[1,"subtype"] <- sample.cluster
  }
  if (TCGA.cancer == "TCGA-THYM") {
    nb.coeff.res[1,"purity"] <- NA
    purity <- NULL
    LF <- NA
  } else {
    nb.coeff.res[1,"purity"] <- mean(purity)
    nb.coeff.res[1,"LF"] <- mean(LF)
  }
  # 4.3) Check sample has germline VCF file
  if (!TCGA.barcode%in%colnames(transcripts.germ.mut)) {
    nb.coeff.res[1,"error"] <- "VCF germline file not found"
    print("VCF germline file not found, skipping...")
    next
  }
  
  # 4.4) For each NMD geneset
  for (NMD.gene.set.name in NMD.gene.sets.names) {
    
    start_time <- Sys.time()
    nb.res <- NA
    print(paste0("------------------------",toupper(NMD.gene.set.name),"------------------------------"))
    if (exists(NMD.gene.set.name)) {
      NMD.gene.set <- eval(parse(text = paste0(NMD.gene.set.name)))
      # RNAseq matrix for the sample
      RNAseq.filt.sample <- RNAseq.matrix.raw.filt[,colnames(RNAseq.matrix.raw.filt)%in%TCGA.barcode, drop = FALSE]
    } else {next}
    # 4.5) Filter sample-specific somatic CNA/mutations
    if ("somatic_transcripts_filt"%in%nextflow.params) {
      if (sample == 1) {print("... Transcripts filtering (somatic CNV/mutations for each sample) ...")}
      RNAseq.filt <- somatic_transcripts_filtering(RNAseq.matrix = RNAseq.filt.sample, transcripts.filtering.matrix = transcripts.to.remove, TCGA.barcode = TCGA.barcode)
    }
    # 4.6) Filter sample-specific germline mutations
    if ("germline_transcripts_filt"%in%nextflow.params) {
      if (sample == 1) {print("... Transcripts filtering (germline mutations for each sample) ...")}
      RNAseq.filt <- germline_transcripts_filtering(RNAseq.matrix = RNAseq.filt, transcripts.filtering.matrix = transcripts.germ.mut)
    }    
    if (nrow(RNAseq.filt) <= 1 || is.null(RNAseq.filt) & ( NMD.gene.set.name == "NMD.global.2.shared" )) {
      nb.coeff.res[1,"error"] <- "<=1 total genes"
      next
    }
    # 4.7) Merge Gene expression to the NMD targets-controls gene set
    NMD.gene.set.gene.exp <- merge(NMD.gene.set,RNAseq.filt, by.x= "ensembl_transcript_id", by.y = "row.names",all.x = TRUE)
    colnames(NMD.gene.set.gene.exp)[colnames(NMD.gene.set.gene.exp) %in% TCGA.barcode] <- "raw_gene_exp"
    # 4.8) Add transcript length
    NMD.gene.set.gene.exp <- merge(NMD.gene.set.gene.exp,ensembl.v88.transcripts.length, by.x= "ensembl_transcript_id", by.y = "transcript_id",all.x = TRUE)
    # 4.9) Filter transcripts
    # Remove transcripts without gene expression (filtered ones)
    NMD.gene.set.gene.exp.filt <- NMD.gene.set.gene.exp[!is.na(NMD.gene.set.gene.exp[,"raw_gene_exp"]),]
    # Remove genes where only 1 transcript has expression
    genes.ids.to.keep <- NMD.gene.set.gene.exp.filt[duplicated(NMD.gene.set.gene.exp.filt$ensembl_gene_id),"ensembl_gene_id"]
    NMD.final.df <- NMD.gene.set.gene.exp.filt[NMD.gene.set.gene.exp.filt$ensembl_gene_id %in% genes.ids.to.keep,]
    # Number of NMD.global.2.shared" NMD targets
    if ( NMD.gene.set.name == "NMD.global.2.shared" ) {
      num.NMDtargets <- length(unique(NMD.final.df$ensembl_gene_id))
      genes.size <- num.NMDtargets
      nb.coeff.res[1,"num_NMD_targets"] <- num.NMDtargets
      print(paste0("Final number of NMD.global.2.shared NMD targets --> ", num.NMDtargets))
      total.num.NMDtargets.sample <- c(total.num.NMDtargets.sample,num.NMDtargets)
      if (num.NMDtargets == 0) {
        nb.coeff.res[1,"error"] <- "0 NMD targets"
        #next
      }
    } else {
      genes.size <- 50
      num.NMDtargets <- length(unique(NMD.final.df$ensembl_gene_id))
      if (num.NMDtargets > genes.size ) {
        num.NMDtargets <- genes.size
      }
      nb.coeff.res[1,"num_NMD_targets"] <- num.NMDtargets
    }
    # Save input data.frame
    res.path <- gsub("\\[X1\\]",TCGA.cancer, paste0(NMD.endogenous.transcripts.path))
    res.path <- gsub("\\[X2\\]",gsub("\\.","-",TCGA.barcode), res.path)
    dir.create(res.path,showWarnings = FALSE)
    NMD.transcripts.res.path <- paste0(res.path,paths[paths$folder_or_object=="NMD_endogenous_transcripts","path_or_filename"])
    NMD.transcripts.res.path <- gsub("\\[X\\]",gsub("\\.","_",NMD.gene.set.name), NMD.transcripts.res.path)
    write.table(NMD.final.df, file = NMD.transcripts.res.path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
    # Save nb.coeff.res empty table
    NB.reg.res.path <- paste0(res.path,paths[paths$folder_or_object=="NB_regression_results","path_or_filename"])
    NB.reg.res.path <- gsub("\\[X\\]",gsub("\\.","_",NMD.gene.set.name), NB.reg.res.path)
    cols <- c(other.columns,c(NMD.gene.set.name,paste0(NMD.gene.set.name,c(".pval",".sd",".CI_2.5",".CI_97.5"))))
    nb.coeff.res.final <- nb.coeff.res[,cols]
    write.table(nb.coeff.res.final, file = NB.reg.res.path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

    # Send job to perform NB regression
    command <- paste0(scripts.WD,"endogenous_NMD_efficiency_regression.sh endogenous_NMD_efficiency_PATHS_cluster.txt ",TCGA.cancer," ",gsub("\\.","-",TCGA.barcode)," ",NMD.gene.set.name)
    # Save the command in a new file
    sbatch.commands[(nrow(sbatch.commands)+sample),"sbatch_command"] <- command

    # 4.11) Save output
    end_time <- Sys.time()
    print(paste0(NMD.gene.set.name," time taken --> ",end_time - start_time))
  }
  end_time_total <- Sys.time()
  print(paste0("TOTAL time take for the sample -->",end_time_total - start_time_total))
}

output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA.cancer,"/",TCGA.cancer,"_endogenous_sbatch_commands_samples.txt")
print(paste0("output here -->",output_path))
write.table(na.omit(sbatch.commands), file = output_path, col.names = FALSE, row.names = FALSE, quote = FALSE)



