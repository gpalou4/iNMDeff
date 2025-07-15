rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

germline_transcripts_filtering <- function(RNAseq.matrix, transcripts.filtering.matrix) {
  germ.mut <- "stopgain|[^non]frameshift deletion|[^non]frameshift insertion|splicing|nonframeshift deletion|nonframeshift insertion|startloss|stoploss"
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

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/GTEx/endogenous_NMD_efficiency_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#GTEx.tissue <- "Liver"
#nextflow.params <- c("germline_transcripts_filt")

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.tissue <- args[2]
nextflow.params <- args[3:length(args)]

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_raw_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
nb.coeff.res.path <- paths[paths$folder_or_object=="NB_results","path_or_filename"]
transcripts.germ.mut.path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
NMD.endogenous.transcripts.path <- paths[paths$folder_or_object=="NMD_endogenous_transcripts_path","path_or_filename"]
scripts.WD  <- "/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/GTEx/"

print(paste0("GTEx tissue --> ",GTEx.tissue))
print("nextflow parameters:")
print(nextflow.params)

# 1.1) RNAseq matrix raw counts
RNAseq.matrix.raw <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])),
                                              header = TRUE, row.names = 1)
# Delete first column (ENSEMBL Gene ID)
RNAseq.matrix.raw <- RNAseq.matrix.raw[,-1]
# Delete genes in PAR Y region to avoid having same ENSEMBL transcript ID twice
RNAseq.matrix.raw <- RNAseq.matrix.raw[-grep("PAR",rownames(RNAseq.matrix.raw)),]
# Remove ID version
rownames(RNAseq.matrix.raw) <- gsub("(.*)\\..*","\\1",rownames(RNAseq.matrix.raw))
# Round values
RNAseq.matrix.raw <- round(RNAseq.matrix.raw)

# 1.2) Transcripts length
ensembl.v88.transcripts.length <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_transcripts_length","path_or_filename"]), 
                                header = TRUE, sep = "\t")
# 1.3) Protein coding transcripts
ensembl_v88_coding_transcripts <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_coding_transcripts.txt",
                                  header = FALSE, sep = "\t",colClasses = "vector")   
GTEx_RNAseq_raw_pantissue_filt <- GTEx_RNAseq_raw_pantissue_filt[ rownames(GTEx_RNAseq_raw_pantissue_filt) %in% 
                  unique(ensembl_v88_coding_transcripts$V2),]                                   
dim(GTEx_RNAseq_raw_pantissue_filt)        

# 1.5) transcripts with germline mutations
transcripts.germ.mut <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(transcripts.germ.mut.path,paths[paths$folder_or_object=="transcripts_germ_mut","path_or_filename"])),
                                    header = TRUE, sep = "\t", row.names = 1)
transcripts.germ.mut[] <- lapply(transcripts.germ.mut , as.character)

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

# 2.1) Filter low-expressed genes
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw[rowSums(log2(RNAseq.matrix.raw) >= 1) >= round(length(colnames(RNAseq.matrix.raw)) * 0.50),]

# 2.2) Filter out non-coding transcripts
RNAseq.matrix.raw.filt <- RNAseq.matrix.raw.filt[rownames(RNAseq.matrix.raw.filt)%in%ensembl.v88.coding.transcripts$V1,]

# 2.3) Offset 1 (transcript length)
# Already done

# 2.4) Filter NMD target genes outliers using cell line data
NMD.gene.sets.names <- c("NMD.Colombo","NMD.Karousis","NMD.Tani","NMD.Courtney","NMD.ensembl","NMD.global",
                        "NMD.global.all.shared","NMD.global.2.shared","SMG6","SMG7","non.NMD.neg.control",
                        "non.NMD.neg.control.with.NMD.features")
NMD.gene.sets.names <- c("non.NMD.neg.control")

# 4) Negative Binomial regression on raw count data 
# Dataframe for coefficients for different NMD gene sets
other.columns <- c("sample","num_NMD_targets","error")
df.cols <- length(other.columns)+(length(NMD.gene.sets.names)*5)
# Sbatch command table
sbatch.commands <- data.frame(sbatch_command = NA)

for (sample in seq(1,ncol(RNAseq.matrix.raw))) {
    
  nb.coeff.res <- data.frame(matrix(ncol=df.cols,nrow=0, dimnames=list(NULL, c(other.columns,NMD.gene.sets.names,paste0(NMD.gene.sets.names,".pval"),
                                  paste0(NMD.gene.sets.names,".sd"),paste0(NMD.gene.sets.names,".CI_2.5"),paste0(NMD.gene.sets.names,".CI_97.5")))))
  start_time_total <- Sys.time()
  GTEx.full.barcode <- colnames(RNAseq.matrix.raw)[sample]
  GTEx.barcode <- gsub("\\.","-",gsub("(GTEX\\.\\w{4,5})\\..*","\\1",GTEx.full.barcode))
  nb.coeff.res[1,"sample"] <- GTEx.barcode
  print(paste("SAMPLE", sample, " ----> ",GTEx.barcode),sep = "")
  # Sample transcripts with germline variants
  if (!GTEx.full.barcode%in%colnames(transcripts.germ.mut)) {
    nb.coeff.res[1,"error"] <- "VCF germline file not found"
    print("VCF germline file not found, skipping...")
    next
  } else {
    transcripts.germ.mut.sample <-  transcripts.germ.mut[,GTEx.full.barcode, drop = FALSE]
  }
  # 4.2) For each NMD geneset
  for (NMD.gene.set.name in NMD.gene.sets.names) {
    
    start_time <- Sys.time()
    nb.res <- NA
    print(paste0("------------------------",toupper(NMD.gene.set.name),"------------------------------"))
    if (exists(NMD.gene.set.name)) {
      NMD.gene.set <- eval(parse(text = paste0(NMD.gene.set.name)))
      # RNAseq matrix for the sample
      RNAseq.filt.sample <- RNAseq.matrix.raw.filt[,GTEx.full.barcode, drop = FALSE]
    } else {next}
    # 4.3) Filter sample-specific germline mutations
    if ("germline_transcripts_filt"%in%nextflow.params) {
      if (sample == 1) {print("... Transcripts filtering (germline mutations for the sample) ...")}
      RNAseq.filt <- germline_transcripts_filtering(RNAseq.matrix = RNAseq.filt.sample, transcripts.filtering.matrix = transcripts.germ.mut.sample)
    }    
    if (nrow(RNAseq.filt) <= 1 || is.null(RNAseq.filt)) {
      nb.coeff.res[1,"error"] <- "<=1 total genes"
      next
    }
    # 4.4) Merge Gene expression to the NMD targets-controls gene set
    NMD.gene.set.gene.exp <- merge(NMD.gene.set,RNAseq.filt, by.x= "ensembl_transcript_id", by.y = "row.names",all.x = TRUE)
    colnames(NMD.gene.set.gene.exp)[colnames(NMD.gene.set.gene.exp) %in% GTEx.full.barcode] <- "raw_gene_exp"
    # 4.5) Add transcript length
    NMD.gene.set.gene.exp <- merge(NMD.gene.set.gene.exp,ensembl.v88.transcripts.length, by.x= "ensembl_transcript_id", by.y = "transcript_id",all.x = TRUE)
    # 4.6) Filter transcripts
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
      #total.num.NMDtargets.sample <- c(total.num.NMDtargets.sample,num.NMDtargets)
      if (num.NMDtargets == 0) {
        nb.coeff.res[1,"error"] <- "0 NMD targets"
        next
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
    res.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(NMD.endogenous.transcripts.path))
    res.path <- gsub("\\[X2\\]",gsub("\\.","-",GTEx.barcode), res.path)
    dir.create(res.path,showWarnings = FALSE, recursive = TRUE)
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
    command <- paste0(scripts.WD,"endogenous_NMD_efficiency_regression.sh endogenous_NMD_efficiency_PATHS_cluster.txt ",GTEx.tissue," ",gsub("\\.","-",GTEx.barcode)," ",NMD.gene.set.name)
    # Save the command in a new file
    sbatch.commands[(nrow(sbatch.commands)+sample),"sbatch_command"] <- command
    # print(paste0("SENDING SAMPLE JOB --> ",command))
    # system(command)

    # 4.7) Save output
    end_time <- Sys.time()
    print(paste0(NMD.gene.set.name," time taken --> ",end_time - start_time))
  }
  end_time_total <- Sys.time()
  print(paste0("TOTAL time take for the sample -->",end_time_total - start_time_total))
}

output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/NMDeff/",GTEx.tissue,"/",GTEx.tissue,"_endogenous_sbatch_commands_samples.txt")
write.table(na.omit(sbatch.commands), file = output_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
