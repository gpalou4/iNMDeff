rm(list=ls())
library("biomaRt")
library("matrixStats")
library("ggplot2")
# Read GTF
library("rtracklayer")
# Ven Diagrams
library("VennDiagram")
library("UpSetR")

# 1) Data

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/list_NMD_transcripts/list_ensembl_transcripts_NMD_features_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
# NMD targets
NMD.targets.gene.symbols.updated <- paths[paths$folder_or_object=="NMD_targets_gene_symbols_updated","path_or_filename"]
NMD.targets.raw <- paths[paths$folder_or_object=="NMD_targets_raw","path_or_filename"]
#NMD.targets.ensembl.shared.path <- paths[paths$folder_or_object=="NMD_targets_ensembl_shared","path_or_filename"]
NMD.targets.ensembl.path <- paths[paths$folder_or_object=="NMD_targets_ensembl","path_or_filename"]
# Results
NMD.targets.res <- paths[paths$folder_or_object=="NMD_genesets_res_path","path_or_filename"]
TCGA.names.path <- paths[paths$folder_or_object=="TCGA_names_path","path_or_filename"]
RNAseq.path <- paths[paths$folder_or_object=="RNAseq_TPM_path","path_or_filename"]

# 1.1) ENSEMBL with NMD features
ensembl.NMD.features <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_NMD_features","path_or_filename"]), 
                                    header = TRUE, sep = "\t", colClasses = "character")


####3
ensembl.NMD.features <- read.table(file = paste0(NMD_targets_ensembl_path,"ensembl_NMD_features_all.txt"),
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# ensembl_NMD_features_filt <- ensembl_NMD_features[,c("ensembl_transcript_id","stop_codon_to_3UTR_splice_site")]
# colnames(ensembl_NMD_features_filt)[2] <- "stop_codon_to_3UTR_splice_site_new"


# ensembl.NMD.features <- merge(ensembl.NMD.features,ensembl_NMD_features_filt, by = "ensembl_transcript_id", all.x = TRUE)

# ensembl.NMD.features$stop_codon_to_3UTR_splice_site <- ensembl.NMD.features$stop_codon_to_3UTR_splice_site_new

######


# 1.2) ENSEMBL transcripts IDs hg19 GTF
ensembl.v88.gtf <- rtracklayer::import(paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
ensembl.v88.gtf <- as.data.frame(ensembl.v88.gtf)
# Remove versions from IDs
ensembl.v88.gtf[,c("gene_id","transcript_id")] <- sapply(ensembl.v88.gtf[,c("gene_id","transcript_id")], function(col) { 
  gsub("(.*)\\..*","\\1", col)
}) 

# Conversor Table for Gene Symbols - ENSEMBL gene ID
ensembl.v88.genesymbol <- ensembl.v88.gtf[,c("gene_name","gene_id")]
ensembl.v88.genesymbol <- ensembl.v88.genesymbol[!duplicated(ensembl.v88.genesymbol),]

# X) Conversor Affymetrix human genome GeneChip U133 plus array - ensembl table
# ensembl.affyU133plus.conversion <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensemblToU133Plus2","path_or_filename"]), 
#                                            header = FALSE, sep = "\t")
# ensembl.affyU133plus.conversion[,1] <- sub("\\..*","", ensembl.affyU133plus.conversion[,1])
# colnames(ensembl.affyU133plus.conversion) <- c("ensembl","AFFY_U133_PLUS")

# 1.3) TCGA RNAseq TPM for each tissue (I still don't have this for ENSEMBL)

# TCGA.tissues <- as.character(read.table(file = paste0(TCGA.names.path,paths[paths$folder_or_object=="TCGA_names","path_or_filename"]),sep = "\t")$V1)
#TCGA.tissues <- c("TCGA-LIHC")

# print("1) RNAseq TPM")
# RNAseq.sample.names <- list()
# RNAseq.TCGA.TPM.all <- c()
# for (i in seq(1:length(TCGA.tissues))) {
#   TCGA.cancer <- as.character(TCGA.tissues[i])
#   print(paste0(i," --> ",TCGA.cancer))
#   RNAseq.error <- FALSE
#   tryCatch( {   RNAseq.TCGA.TPM <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])),
#                                               header = TRUE, sep = "\t", row.names = 1)
#   colnames(RNAseq.TCGA.TPM) <- gsub("\\.","-",substr(colnames(RNAseq.TCGA.TPM),1,12))
#   # It needs to be multiplied by 1 million to be transformed to TPM values
#   #RNAseq.TCGA.TPM <- RNAseq.TCGA.TPM*1000000
#   # keep track of sample names for that cancer
#   RNAseq.sample.names[[TCGA.cancer]] <- colnames(RNAseq.TCGA.TPM)
#   # Join matrices in one
#   if (length(RNAseq.TCGA.TPM.all)==0) {RNAseq.TCGA.TPM.all <- RNAseq.TCGA.TPM}
#   else{RNAseq.TCGA.TPM.all <- cbind(RNAseq.TCGA.TPM.all,RNAseq.TCGA.TPM)}
#   # Check rownames are the same
#   print(paste0("ENSEMBL transcripts IDs are the same? ",table(rownames(RNAseq.TCGA.TPM)==rownames(RNAseq.TCGA.TPM.all))))
#   }
#   ,error = function(e) {
#     RNAseq.error <<- TRUE
#     print(e)
#     })
#   if (isTRUE(RNAseq.error)) {next}
# }
# print("Dimensions -->")
# print(dim(RNAseq.TCGA.TPM.all))

# 1.5) Pantissue TPM RNAseq
# 1.5.1) TCGA
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt",
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

RNAseq.TCGA.TPM.all <- TCGA_RNAseq_TPM_pancancer

# 1.4) NMD targets
# Gene Symbols and ENSEMBL genes IDs will be transformed to ENSEMBL transcripts IDs using the ENSEMBL v88 GTF file

# 1.4.1) Tani
# Updated gene symbols for UPF1

# Group A
Tani.UPF1.dir.unconf.gene.symbols <- read.table(file = paste0(NMD.targets.gene.symbols.updated,"Tani_UPF1_dir_unconf_gene_symbols_updated.txt"), 
                                                   header = FALSE, sep = "\t", colClasses = "vector")$V1
# Group B
Tani.UPF1.indir.conf.gene.symbols <- read.table(file = paste(NMD.targets.gene.symbols.updated,"Tani_UPF1_ind_conf_gene_symbols_updated.txt", sep = ""), 
                                                   header = FALSE, sep = "\t", colClasses = "vector")$V1
# Group C
Tani.UPF1.dir.conf.gene.symbols <- read.table(file = paste(NMD.targets.gene.symbols.updated, "Tani_UPF1_dir_conf_gene_symbols_updated.txt", sep = ""), 
                                                 header = FALSE, sep = "\t", colClasses = "vector")$V1
# All merged
Tani.UPF1.all.gene.symbols <- read.table(file = paste(NMD.targets.gene.symbols.updated, "Tani_UPF1_all_gene_symbols_updated.txt", sep = ""), 
                                                 header = FALSE, sep = "\t", colClasses = "vector")$V1

# Convert Gene Symbols to ENSEMBL genes ID
Tani.UPF1.all.gene.symbols[!Tani.UPF1.all.gene.symbols%in%ensembl.v88.genesymbol$gene_name]
Tani.UPF1.all.ensembl.gene <- ensembl.v88.genesymbol[ensembl.v88.genesymbol$gene_name%in%Tani.UPF1.all.gene.symbols,"gene_id"]
Tani.UPF1.dir.unconf.ensembl.gene <- ensembl.v88.genesymbol[ensembl.v88.genesymbol$gene_name%in%Tani.UPF1.dir.unconf.gene.symbols,"gene_id"]
Tani.UPF1.dir.conf.ensembl.gene <- ensembl.v88.genesymbol[ensembl.v88.genesymbol$gene_name%in%Tani.UPF1.dir.conf.gene.symbols,"gene_id"]
Tani.UPF1.indir.conf.ensembl.gene <- ensembl.v88.genesymbol[ensembl.v88.genesymbol$gene_name%in%Tani.UPF1.indir.conf.gene.symbols,"gene_id"]

# 1.4.2) Schmidt
# Updated gene symbols for SMG6

Schmidt.SMG6.gene.symbols <- read.table(file = paste(NMD.targets.gene.symbols.updated, "Schmidt_SMG6_gene_symbol_updated.txt", sep = ""), 
                                                 header = TRUE, sep = "\t", colClasses = "vector")

# Convert Schmidt Gene Symbols to ENSEMBL gene ID
Schmidt.SMG6.gene.symbols[!Schmidt.SMG6.gene.symbols$gene_symbol%in%ensembl.v88.genesymbol$gene_name,]
Schmidt.SMG6.ensembl.gene <- ensembl.v88.genesymbol[ensembl.v88.genesymbol$gene_name%in%Schmidt.SMG6.gene.symbols$gene_symbol,"gene_id"]

# 1.4.3) Colombo
# ENSEMBL genes IDs for SMG6, SMG7 and UPF1
Colombo.SMG6.ensembl.gene <- read.table(file = paste(NMD.targets.raw, "Colombo_SMG6_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
Colombo.SMG7.ensembl.gene <- read.table(file = paste(NMD.targets.raw, "Colombo_SMG7_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
Colombo.UPF1.ensembl.gene <- read.table(file = paste(NMD.targets.raw, "Colombo_UPF1_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")

# 1.4.4) Karousis
# ENSEMBL genes IDs for NMD global
Karousis.ensembl.gene <- read.table(file = paste(NMD.targets.raw, "Karousis_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
# 1.4.5) Courtney
# ENSEMBL genes IDs for NMD global
Courtney.ensembl.gene <- read.table(file = paste(NMD.targets.raw, "Courtney_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")

# 1.4.6) NMD targets from ensembl gtf
ensembl.v88.gtf.NMD.targets <- unique(ensembl.v88.gtf[ensembl.v88.gtf$transcript_type%in%"nonsense_mediated_decay","gene_id"])

# 2) NMD genesets

# Remove Tani.UPF1.dir.unconf.ensembl.gene from all the NMD genesets
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Colombo.SMG7.ensembl.gene$gene)
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Colombo.SMG6.ensembl.gene$gene)
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Colombo.UPF1.ensembl.gene$gene)
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Karousis.ensembl.gene$gene_id)
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Schmidt.SMG6.ensembl.gene)
table(Tani.UPF1.dir.unconf.ensembl.gene%in%Courtney.ensembl.gene$Ensembl_Gene_ID)

Tani.UPF1.all.ensembl.gene <- Tani.UPF1.all.ensembl.gene[!Tani.UPF1.all.ensembl.gene%in%Tani.UPF1.dir.unconf.ensembl.gene]
Colombo.SMG7.ensembl.gene <- Colombo.SMG7.ensembl.gene[!Colombo.SMG7.ensembl.gene$gene%in%Tani.UPF1.dir.unconf.ensembl.gene,]
Colombo.SMG6.ensembl.gene <- Colombo.SMG6.ensembl.gene[!Colombo.SMG6.ensembl.gene$gene%in%Tani.UPF1.dir.unconf.ensembl.gene,]
Colombo.UPF1.ensembl.gene <- Colombo.UPF1.ensembl.gene[!Colombo.UPF1.ensembl.gene$gene%in%Tani.UPF1.dir.unconf.ensembl.gene,]
Karousis.ensembl.gene <- Karousis.ensembl.gene[!Karousis.ensembl.gene$gene_id%in%Tani.UPF1.dir.unconf.ensembl.gene,]
Schmidt.SMG6.ensembl.gene <- Schmidt.SMG6.ensembl.gene[!Schmidt.SMG6.ensembl.gene%in%Tani.UPF1.dir.unconf.ensembl.gene]
Courtney.ensembl.gene <- Courtney.ensembl.gene[!Courtney.ensembl.gene$Ensembl_Gene_ID%in%Tani.UPF1.dir.unconf.ensembl.gene,]

# 2.1) SMG6 (genes shared in Colombo + Schmidt)

SMG6 <- Colombo.SMG6.ensembl.gene[Colombo.SMG6.ensembl.gene$gene%in%Schmidt.SMG6.ensembl.gene,"gene"]

# 2.2) SMG7 (genes in Colombo only)

SMG7 <- Colombo.SMG7.ensembl.gene$gene

# 2.3) NMD_tani (Tani UPF1, EJC indep)

NMD.Tani <- Tani.UPF1.all.ensembl.gene

# 2.4) NMD_colombo (Colombo UPF1 + SMG6 + SMG7, EJC indep/dep)

NMD.Colombo <- unique(c(Colombo.SMG6.ensembl.gene$gene, Colombo.SMG7.ensembl.gene$gene, Colombo.UPF1.ensembl.gene$gene))

# 2.5) NMD_karousis (Karousis UPF1 + SMG6 + SMG7, EJC indep/dep)

NMD.Karousis <- Karousis.ensembl.gene$gene_id

# 2.6) NMD_courtney (Karousis UPF1, EJC indep/dep)

NMD.Courtney <- Courtney.ensembl.gene$Ensembl_Gene_ID

# 2.7) NMD global (all genes except ensembl)

NMD.global <- unique(c(Tani.UPF1.all.ensembl.gene, Colombo.SMG6.ensembl.gene$gene, Colombo.SMG7.ensembl.gene$gene, Colombo.UPF1.ensembl.gene$gene, 
Schmidt.SMG6.ensembl.gene, Karousis.ensembl.gene$gene_id, Courtney.ensembl.gene$Ensembl_Gene_ID))

# 2.8) NMD_global_all_shared (all genes shared in all papers except Schmidt)
all.NMD.genesets <- list(NMD_Tani = NMD.Tani, NMD_Colombo = NMD.Colombo, NMD_Karousis = NMD.Karousis, NMD_Courtney = NMD.Courtney)
NMD.global.all.shared <- Reduce(intersect, all.NMD.genesets)

# 2.9) NMD_global_2_shared (genes shared in >= 2 papers)
all.NMD.genesets.vector <- unlist(all.NMD.genesets)
NMD.global.2.shared <- unique(all.NMD.genesets.vector [duplicated(all.NMD.genesets.vector )])

# 2.10) NMD ensembl

NMD.ensembl <- ensembl.v88.gtf.NMD.targets

# 2.10) For the Negative Control

all.NMD.genes <- unique(c(NMD.global,Tani.UPF1.dir.unconf.ensembl.gene,NMD.ensembl))
non.NMD.genes <- unique(ensembl.NMD.features[!ensembl.NMD.features$ensembl_gene_id%in%all.NMD.genes,"ensembl_gene_id"])

# 2.11) UPF3B independent genes symbols #####

# UPF3B.affymetrix <- read.table(file = paste(NMD.targets.raw, "UPF3B_affymetrix.txt", sep = ""),
                               # header = FALSE, sep = "\t", colClasses = "vector")$V1
# Convert to ensembl
# UPF3B.ensembl <- ensembl.affyU133plus.conversion[ensembl.affyU133plus.conversion$AFFY_U133_PLUS%in%UPF3B.affymetrix,"ensembl"]
# Convert to Gene Symbols
# UPF3B.gene.symbols <- unique(ensembl.gene.symbol.conversion.updated.updated[ensembl.gene.symbol.conversion.updated.updated$FeatureID%in%UPF3B.ensembl,"Approved.symbol"])

# 3) Ven Diagrams of shared NMD genesets

library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
        x = all.NMD.genesets,
        TACegory.names = names(all.NMD.genesets),
        filename = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_venn_diagram","path_or_filename"]),
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 1000 , 
        width = 1000 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .8,
        fontface = "bold",
        fontfamily = "sans",

        # Set names
        TAC.cex = 0.6,
        TAC.fontface = "bold",
        TAC.fontfamily = "sans",
)

# UPF1 intersect
png(file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_intersection_UPF1","path_or_filename"]), width = 4000, height = 3000, res = 300)
NMD.genesets.list = list(Tani_UPF1_dir_conf = Tani.UPF1.dir.conf.ensembl.gene, Tani_UPF1_indir_conf = Tani.UPF1.indir.conf.ensembl.gene,
                          Colombo_UPF1 = Colombo.UPF1.ensembl.gene$gene, NMD_Karousis = NMD.Karousis, NMD_Courtney = NMD.Courtney, NMD_ensembl = NMD.ensembl)
p <- upset(fromList(NMD.genesets.list), nsets = length(NMD.genesets.list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD.genesets.list), keep.order = TRUE)
print(p)
dev.off()

# SMG6/7 intersect
png(file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_intersection_SMG6_7","path_or_filename"]), width = 4000, height = 3000, res = 300)
NMD.genesets.list = list(Schmidt_SMG6 = Schmidt.SMG6.ensembl.gene, Colombo_SMG6 = Colombo.SMG6.ensembl.gene$gene, Colombo_SMG7 = Colombo.SMG7.ensembl.gene$gene,
                          Colombo_UPF1 = Colombo.UPF1.ensembl.gene$gene, NMD_Karousis = NMD.Karousis, NMD_Courtney = NMD.Courtney, NMD_ensembl = NMD.ensembl)
p <- upset(fromList(NMD.genesets.list), nsets = length(NMD.genesets.list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD.genesets.list), keep.order = TRUE)
print(p)
dev.off()

# All NMD genesets separately
png(file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_intersection_all","path_or_filename"]), width = 4000, height = 3000, res = 300)
NMD.genesets.list = list(Tani_dir_conf = Tani.UPF1.dir.conf.ensembl.gene, Tani_indir_conf = Tani.UPF1.indir.conf.ensembl.gene,
                          Schmidt_SMG6 = Schmidt.SMG6.ensembl.gene, Colombo_SMG6 = Colombo.SMG6.ensembl.gene$gene, Colombo_SMG7 = Colombo.SMG7.ensembl.gene$gene,
                          Colombo_UPF1 = Colombo.UPF1.ensembl.gene$gene, NMD_Karousis = NMD.Karousis, NMD_Courtney = NMD.Courtney, NMD_ensembl = NMD.ensembl)
p <- upset(fromList(NMD.genesets.list), nsets = length(NMD.genesets.list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD.genesets.list), keep.order = TRUE)
print(p)
dev.off()

# NMD genesets
png(file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_intersection_final","path_or_filename"]), width = 4000, height = 3000, res = 300)
# NMD.genesets.list = list(SMG6 = SMG6, SMG7 = SMG7, NMD_tani = NMD.Tani, NMD_Colombo = NMD.Colombo,
#                           NMD_global = NMD.global, NMD_global_all_shared = NMD.global.all.shared, NMD_global_2_shared = NMD.global.2.shared,
#                           NMD_Karousis = NMD.Karousis, NMD_Courtney = NMD.Courtney, NMD_ensembl = NMD.ensembl)
NMD.genesets.list = list(NMD_tani = NMD.Tani, NMD_Colombo = NMD.Colombo, NMD_Karousis = NMD.Karousis, 
                        NMD_Courtney = NMD.Courtney, NMD_ensembl = NMD.ensembl)
p <- upset(fromList(NMD.genesets.list), nsets = length(NMD.genesets.list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD.genesets.list), keep.order = TRUE)
print(p)
dev.off()

################################################################################################################################################

# 3) Filter NMD genesets using our NMD features

NMD.gene.sets.names <- c("all.ensembl.genes","non.NMD.genes","SMG6","SMG7","NMD.global", "NMD.ensembl",
                          "NMD.global.all.shared","NMD.global.2.shared","NMD.Karousis","NMD.Courtney","NMD.Tani","NMD.Colombo")
NMD.features <- c("uORFs","uORFs_1_30_bp","uORFs_2_30_bp","UTR3_EJC","UTR3_EJC_50nt","uORFs_UTR3_EJC","non_NMD","uORFs_num_translated","uORFs_num_trans_reinit","UTR3_GC_content",
"penultime_codon_GCG","penultime_codon_AGG","penultime_codon_CTG","penultime_codon_ACA_control","penultime_codon_CCA_control","penultime_codon_TAC_control")

# Add gene expression for all ENSEMBL transcripts IDs
RNAseq.matrix.TPM.filt <- RNAseq.TCGA.TPM.all[rownames(RNAseq.TCGA.TPM.all)%in%ensembl.NMD.features$ensembl_transcript_id,]
ensembl.transcripts <- rownames(RNAseq.matrix.TPM.filt)
gene.exp <- rowMedians(as.matrix(RNAseq.matrix.TPM.filt),na.rm = TRUE)
gene.exp.df <- data.frame(ensembl_id = ensembl.transcripts, gene_exp = gene.exp)
ensembl.NMD.features <- merge(ensembl.NMD.features,gene.exp.df, by.x = "ensembl_transcript_id", by.y = "ensembl_id", all.x = TRUE)
ensembl.NMD.features[,c("final_consensus","NMD_nonNMD_ratio")] <- NA
print(ensembl.NMD.features[1:3,])

ORFs_num_and_length <- function(ORFs.lengths, min.bp, min.ORFs) {
  
  ORFs.lengths.list <- strsplit(ORFs.lengths, ",")
  ORFs.bool <- lapply(ORFs.lengths.list, function(x) {
    # print(sum(as.numeric(x)))
    x <- as.numeric(x)
    if (length(x) == 0) {
      return(FALSE)
    } else if (is.na(x)) {
      return(FALSE)
    } else if ( sum(x >= min.bp) >= min.ORFs ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  return(unlist(ORFs.bool))
}

check_NMD_features <- function(ensembl_NMD_features, set = NULL, control, gene_exp_NMD_target = 1, gene_exp_control = 3) {

  if (!is.null(set)) {
    ensembl_NMD_features_set <- ensembl_NMD_features[ensembl_NMD_features$ensembl_transcript_id%in%set,]
  } else {
    ensembl_NMD_features_set <- ensembl_NMD_features
  }
  if (control == "yes") {
    # no Start or stop codons
    ensembl_NMD_features_set <- ensembl_NMD_features_set[grep(".*no stop.*|.*no start.*",ensembl_NMD_features_set$comment, invert = TRUE),]
    filt <- (ensembl_NMD_features_set$NMD_event_type == "") & (ensembl_NMD_features_set$gene_exp >= gene_exp_control) & (ensembl_NMD_features_set$uORFs==0 | ensembl_NMD_features_set$uORFs == "")
    controls.index <- which(filt)
    # If no controls
    if (length(controls.index) == 0 ) {
      #print("Gene with no control(s) or low expressed (<= 3), skipping")
      return("NA")
    }
    ensembl_NMD_features_set_filt <- ensembl_NMD_features_set[controls.index,]
    # Choose the control with highest expression
    ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt[which(ensembl_NMD_features_set_filt$gene_exp == max(ensembl_NMD_features_set_filt$gene_exp)),]
    # If more than 1 transcript is chosen (same gene expression), then take one random
    ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_gene_exp[sample(1:nrow(ensembl_NMD_features_set_gene_exp))[1],]
    control.transcript <- ensembl_NMD_features_set_gene_exp$ensembl_transcript_id
    return(control.transcript)
  } else if (control == "no") {
  # no Start or stop codons
    ensembl_NMD_features_set <- ensembl_NMD_features_set[grep(".*no stop.*|.*no start.*",ensembl_NMD_features_set$comment, invert = TRUE),]
    NMD.targets.index <- which( (ensembl_NMD_features_set$NMD_event_type != "") & (ensembl_NMD_features_set$gene_exp >= gene_exp_NMD_target) )
    # If there aren't NMD features or low expression, skip
    if (length(NMD.targets.index) == 0 ) {
      #print("Gene with no NMD target(s) or low expressed (< 1 TPM), skipping")
      return("NA")
  }
  # Criteria
  # 1) 3UTR EJC 50nt > uORF 
  # 2) highest 3'UTR GC content
  ensembl_NMD_features_set_filt <- ensembl_NMD_features_set[NMD.targets.index ,]
  ensembl_NMD_features_set_filt_3UTR_EJC <- ensembl_NMD_features_set_filt[grep("Intronic",ensembl_NMD_features_set_filt$NMD_event_type),]
  ensembl_NMD_features_set_filt_uORFs <- ensembl_NMD_features_set_filt[ensembl_NMD_features_set_filt$NMD_event_type == " | >=2 uORFs",]
  # 50nt
  ensembl_NMD_features_set_filt_3UTR_EJC <- ensembl_NMD_features_set_filt_3UTR_EJC[which(as.numeric(ensembl_NMD_features_set_filt_3UTR_EJC$stop_codon_to_3UTR_splice_site) >= 50),]

  if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) == 0 && nrow(ensembl_NMD_features_set_filt_uORFs) == 0) {
    return("NA")
  } else if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) >= 1 && nrow(ensembl_NMD_features_set_filt_uORFs) >=0 ) { # 3UTR EJC
    # Take the isoform with highest 3'UTR GC content
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_filt_3UTR_EJC[ensembl_NMD_features_set_filt_3UTR_EJC$UTR3_GC_content == max(ensembl_NMD_features_set_filt_3UTR_EJC$UTR3_GC_content),]
    #ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt_3UTR_EJC[which(ensembl_NMD_features_set_filt_3UTR_EJC$gene_exp == min(ensembl_NMD_features_set_filt_3UTR_EJC$gene_exp)),]
  } else if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) == 0 && nrow(ensembl_NMD_features_set_filt_uORFs) >= 1 ) { # uORFs
    # Take the isoform with highest 3'UTR GC content
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_filt_uORFs[ensembl_NMD_features_set_filt_uORFs$UTR3_GC_content == max(ensembl_NMD_features_set_filt_uORFs$UTR3_GC_content),]
    #ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt_uORFs[which(ensembl_NMD_features_set_filt_uORFs$gene_exp == min(ensembl_NMD_features_set_filt_uORFs$gene_exp)),]
  }
  # If more than 1 transcript is kept (same 3'UTR GC content)
  # then take the one with highest number of uORFs
  if (nrow(ensembl_NMD_features_set_3utr_gc) >1) {
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_3utr_gc[ensembl_NMD_features_set_3utr_gc$uORFs == max(ensembl_NMD_features_set_3utr_gc$uORFs),]
  }
  # If more than 1 transcript is kept (same #uORFs)
  # then take the one random
  if (nrow(ensembl_NMD_features_set_3utr_gc) >1) {
  ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_3utr_gc[sample(1:nrow(ensembl_NMD_features_set_3utr_gc))[1],]
  }
  NMD.target <- ensembl_NMD_features_set_3utr_gc$ensembl_transcript_id
  return(NMD.target)
  }
}

NMD_geneset_info_features <- function(ensembl_NMD_features, NMD_geneset, all_NMD_genes = all.NMD.genes, non_NMD_genes = non.NMD.genes) {

  # Keep info
  n <- nrow(ensembl_NMD_features)
  n.uORFs <- as.numeric(table((as.numeric(ensembl_NMD_features$uORFs) > 0))[2])
  n.genes <- length(unique(ensembl_NMD_features$ensembl_gene_id))
  NMD.features.df <- data.frame(NMD_features = c(NMD.features,"fisher_pval","fisher_OR"), value = NA, row.names = 1)
  colnames(NMD.features.df)[1] <- NMD_geneset

  NMD.features.df["transcripts_size",] <- n
  NMD.features.df["genes_size",] <- n.genes
  NMD.features.df["uORFs",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs" & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["uORFs_1_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_1_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["uORFs_2_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_2_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["UTR3_EJC",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR") / n,2)
  NMD.features.df["UTR3_EJC_50nt",] <- round(length(which(ensembl_NMD_features$stop_codon_to_3UTR_splice_site >= 50 & ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR")) / n,2)
  #NMD.features.df["uORFs_UTR3_intron",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs | Intronic-spliced 3UTR") / n,2)
  NMD.features.df["non_NMD",] <- round(sum(ensembl_NMD_features$NMD_event_type == "") / n,2)
  NMD.features.df["UTR3_GC_content",] <- median(as.numeric(ensembl_NMD_features$UTR3_GC_content), na.rm = TRUE)
  NMD.features.df["uORFs_num_translated",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_translated), na.rm = TRUE),2)
  NMD.features.df["uORFs_num_trans_reinit",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_trans_reinit), na.rm = TRUE),2)
  NMD.features.df["penultime_codon_GCG",] <- round(length(grep("GCG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Ala
  #NMD.genesets.res[filt,"penultime_codon_GTG_control"] <- round(length(grep("GTG",ensembl.NMD.features.filt$uORF_penultimate_codon)) / n,6) # Val
  NMD.features.df["penultime_codon_ACA_control",] <- round(length(grep("ACA",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Thr
  NMD.features.df["penultime_codon_AGG",] <- round(length(grep("AGG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Arg
  #NMD.genesets.res[filt,"penultime_codon_AAG_control"] <- round(length(grep("AAG",ensembl.NMD.features.filt$uORF_penultimate_codon)) / n,6) # Lys
  NMD.features.df["penultime_codon_CCA_control",] <- round(length(grep("CCA",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Pro
  NMD.features.df["penultime_codon_CTG",] <- round(length(grep("CTG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Leu
  #NMD.genesets.res[filt,"penultime_codon_CAT_control"] <- round(length(grep("CAT",ensembl.NMD.features.filt$uORF_penultimate_codon)) / n,6) # His
  NMD.features.df["penultime_codon_TAC_control",] <- round(length(grep("TAC",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Tyr
  
  # Fisher test
  ensembl.NMD.features[,"NMD_geneset"] <- NA
  ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%eval(parse(text=paste0(NMD_geneset))),"NMD_geneset"] <- NMD_geneset
  if (NMD_geneset == "non.NMD.genes") {
    ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%all_NMD_genes,"NMD_geneset"] <- "controls"
  } else {
    ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%non_NMD_genes,"NMD_geneset"] <- "controls"
  }

  print(table(ensembl.NMD.features$NMD_geneset))
  # NMD.feature <- names(table(ensembl.NMD.features$NMD_event_type))[4]
  NMD.noNMDfeat <- length(which(ensembl.NMD.features$NMD_event_type=="" & ensembl.NMD.features$NMD_geneset == NMD_geneset))
  NMD.NMDfeat <- length(which(ensembl.NMD.features$NMD_event_type!="" & ensembl.NMD.features$NMD_geneset == NMD_geneset))
  noNMD.noNMDfeat <- length(which(ensembl.NMD.features$NMD_event_type=="" & ensembl.NMD.features$NMD_geneset == "controls"))
  noNMD.NMDfeat <- length(which(ensembl.NMD.features$NMD_event_type!="" & ensembl.NMD.features$NMD_geneset == "controls"))
  fisher.table <- data.frame(NMD_feature = c(NMD.NMDfeat,noNMD.NMDfeat), nonNMD_feature = c(NMD.noNMDfeat,noNMD.noNMDfeat))
  rownames(fisher.table) <- c(NMD_geneset,"controls")
  fisher.res <- fisher.test(fisher.table)
  print(fisher.table)
  print(fisher.res)
  # chisq.test(fisher.table)$expected
  NMD.features.df["fisher_pval",]  <- fisher.res$p.value
  NMD.features.df["fisher_OR",]  <- fisher.res$estimate
  return(t(NMD.features.df))

}

ORFs.num.len.bool <- ORFs_num_and_length(ORFs.lengths = ensembl.NMD.features$uORFs_lengths, min.bp = 30, min.ORFs = 1)
ensembl.NMD.features[,"uORFs_1_30_bp"] <- ORFs.num.len.bool
ORFs.num.len.bool <- ORFs_num_and_length(ORFs.lengths = ensembl.NMD.features$uORFs_lengths, min.bp = 30, min.ORFs = 2)
ensembl.NMD.features[,"uORFs_2_30_bp"] <- ORFs.num.len.bool

all.ensembl.genes <- unique(ensembl.NMD.features$ensembl_gene_id)
NMD.genesets.res <- data.frame()
NMD.genesets.res.filt <- data.frame()

Karousis.ensembl.gene[,c("NMD_target","control")] <- NA



karousis_test <- function(ensembl.gene, ensembl.NMD.features.ensembl.gene) {
  # Check NMD features in all the transcript of the given gene and select the best NMD target
  NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
  # Check NMD features in all the transcript of the given gene and select the best nonNMD control one
  nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
  # Add our criteria to Karousis
  row <- which(Karousis.ensembl.gene$gene_id%in%ensembl.gene)
  Karousis.ensembl.gene[row,"NMD_target"] <<- NMD.target
  Karousis.ensembl.gene[row,"control"] <<- nonNMD.control
}

for (NMD.geneset in rev(NMD.gene.sets.names)) {
  
  print(paste0("_________________________ NMD geneset ---> ", NMD.geneset, "_________________________"))
  # Filter for NMD geneset
  ensembl.NMD.features.filt <- ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%eval(parse(text=paste0(NMD.geneset))),]

  NMD.genesets.res.tmp <- NMD_geneset_info_features(ensembl_NMD_features = ensembl.NMD.features.filt, NMD_geneset = NMD.geneset) 
  NMD.genesets.res <- rbind(NMD.genesets.res,NMD.genesets.res.tmp)
 
  # Check for each gene which isoform is NMD-target and which one could be used as control

  for (ensembl.gene in unique(ensembl.NMD.features.filt$ensembl_gene_id)) {

    NMD.nonNMD.pairs <- data.frame(NMD_target=NA,nonNMD_control=NA)
    ensembl.NMD.features.ensembl.gene <- ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_gene_id%in%ensembl.gene,]

    # Skip if the gene has only 1 isoform
    # For the negative control geneset, just take two random isoforms from the gene 
    if (length(unique(ensembl.NMD.features.ensembl.gene$ensembl_transcript_id)) == 1) {
      print("Gene with only 1 isoform, skipping...")
      next
    } else if (NMD.geneset == "non.NMD.genes") {
        isoforms.index <- which( (ensembl.NMD.features.ensembl.gene$gene_exp >= 3) )
        # If no controls
        if (length(isoforms.index) <= 1 ) {
          print("Gene with less than 2 isoforms expressed (<= 3), skipping")
          next
        } else {
          ensembl.NMD.features.ensembl.gene.filt <- ensembl.NMD.features.ensembl.gene[isoforms.index,]
          isoforms <- sample(ensembl.NMD.features.ensembl.gene.filt$ensembl_transcript_id)[1:2]
          NMD.target <- isoforms[1]
          nonNMD.control <- isoforms[2]
        }
    } else {
      # Check first if the gene is among the Karousis dataset, then obtain its corresponding NMD-nonNMD pair (if possible)
      if (ensembl.gene%in%Karousis.ensembl.gene$gene_id) {
        Karousis.ensembl.transcript <- Karousis.ensembl.gene[Karousis.ensembl.gene$gene_id%in%ensembl.gene,]
        # Check NMD isoform
        NMD.isoforms <- unique(Karousis.ensembl.transcript[Karousis.ensembl.transcript$NMD_isoforms%in%ensembl.NMD.features.ensembl.gene$ensembl_transcript_id,"NMD_isoforms"])
        # Check non_NMD_isoform
        nonNMD.isoforms <- unique(Karousis.ensembl.transcript[Karousis.ensembl.transcript$non_NMD_isoforms%in%ensembl.NMD.features.ensembl.gene$ensembl_transcript_id,"non_NMD_isoforms"])
        if (length(NMD.isoforms) >= 1) {
          # Check NMD features in a specific set of transcripts and select the best NMD target
          NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NMD.isoforms, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
        } else if (length(NMD.isoforms) == 0) {
          # Check NMD features in all the transcripts of the given gene and select the best NMD target
          NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
        }
        if (length(nonNMD.isoforms) >= 1) {
          # Check NMD features in a specific set of transcripts and select the best nonNMD control one
          nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = nonNMD.isoforms, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
        } else if (length(nonNMD.isoforms) == 0) {
          # Check NMD features in all the transcripts of the given gene and select the best nonNMD control one
          nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
        }
        # Check NMD features in all the transcripts of the given gene and select the best NMD target when the Karousis target has no NMD features
        if (NMD.target == "NA") {
          NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
        }
        # Karousis test
        karousis_test(ensembl.gene = ensembl.gene, ensembl.NMD.features.ensembl.gene = ensembl.NMD.features.ensembl.gene)
      } else {
          # Check NMD features in all the transcript of the given gene and select the best NMD target
          NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
          # Check NMD features in all the transcript of the given gene and select the best nonNMD control one
          nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
      }
    }
    # Final consensus
    NMD.nonNMD.pairs$NMD_target <- NMD.target
    NMD.nonNMD.pairs$nonNMD_control <- nonNMD.control

    # Take gene exp and make ratio
    if (NMD.target == "NA" | nonNMD.control == "NA") {
      next
    } else {
        NMD.target.gene.exp <- as.numeric(ensembl.NMD.features.ensembl.gene[ensembl.NMD.features.ensembl.gene$ensembl_transcript_id%in%NMD.target,"gene_exp"])
        control.gene.exp <- as.numeric(ensembl.NMD.features.ensembl.gene[ensembl.NMD.features.ensembl.gene$ensembl_transcript_id%in%nonNMD.control,"gene_exp"])
        NMD.nonNMD.ratio <- NMD.target.gene.exp/control.gene.exp 
      if ( ((NMD.nonNMD.ratio <= 0.9) && (NMD.geneset != "non.NMD.genes")) || (NMD.geneset == "non.NMD.genes") ) {
        print(NMD.nonNMD.pairs)
        ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_transcript_id%in%as.character(NMD.target),"final_consensus"] <- "NMD_target"
        ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_transcript_id%in%as.character(nonNMD.control ),"final_consensus"] <- "control"
      } else {
        print("NMD target has more gene expression than control transcript")
        next
      }
    }
    
    # Do NMD/non-NMD ratios
    #NMD.targets.index <- which(ensembl.NMD.features.ensembl.gene$NMD_event_type != "")
    # NMD.targets.index <- which( (ensembl.NMD.features.ensembl.gene$NMD_event_type != "") & (ensembl.NMD.features.ensembl.gene$gene_exp >=1) )
    # controls.index <- which( (ensembl.NMD.features.ensembl.gene$NMD_event_type == "") & (ensembl.NMD.features.ensembl.gene$gene_exp >=3) )
    # # If no NMD target or controls
    # if (length(NMD.targets.index) == 0 ) {
    #   print("Gene with no NMD target(s) or low expressed (<= 1), skipping")
    #   next
    # } else if (length(controls.index) == 0 ) {
    #   print("Gene with no control(s) or low expressed (<= 3), skipping")
    #   next
    # }
    # # Take the control with highest expression
    # controls.index <- which(ensembl.NMD.features.ensembl.gene$gene_exp == max(ensembl.NMD.features.ensembl.gene[controls.index,"gene_exp"]))
    # NMD.ratios <- expand.grid(NMD.targets = ensembl.NMD.features.ensembl.gene[NMD.targets.index,"ensembl_transcript_id"], controls = ensembl.NMD.features.ensembl.gene[controls.index,"ensembl_transcript_id"])
    # NMD.ratios$NMD_ratio <- NA
    # for (NMD.target in NMD.targets.index) {
    #   NMD.target.gene.exp <- ensembl.NMD.features.ensembl.gene[NMD.target,"gene_exp"]
    #   for (control in controls.index) {
    #     control.gene.exp <- ensembl.NMD.features.ensembl.gene[control,"gene_exp"]
    #     NMD.nonNMD.ratio <- NMD.target.gene.exp/control.gene.exp
    #     filter <- ( (NMD.ratios$NMD.targets%in%ensembl.NMD.features.ensembl.gene[NMD.target,"ensembl_transcript_id"]) & (NMD.ratios$controls%in%ensembl.NMD.features.ensembl.gene[control,"ensembl_transcript_id"]) )
    #     NMD.ratios[filter,"NMD_ratio"] <- NMD.nonNMD.ratio
    #   }
    # }
    # # Choose the best partner with lowest NMD/non-NMD ratio
    # filter <- NMD.ratios$NMD_ratio == min(NMD.ratios$NMD_ratio)
    # ensembl.partners <- NMD.ratios[filter,]
    # best.NMD.ratio <- min(NMD.ratios$NMD_ratio)
    # if ( is.na(best.NMD.ratio) | is.infinite(best.NMD.ratio) ) {
    #   print("NMD ratio non-numeric, skipping")
    #   next
    # }
    # # If ratio > 1, remove the gene
    # if ( best.NMD.ratio >= 1 ) {
    #   print("NMD ratio >=1, skipping")
    #   next
    # }
    # ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_transcript_id%in%as.character(ensembl.partners[,"NMD.targets"]),"final_consensus"] <- "NMD_target"
    # ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_transcript_id%in%as.character(ensembl.partners[,"controls"]),"final_consensus"] <- "control"
    # ensembl.NMD.features.filt[ensembl.NMD.features.filt$ensembl_transcript_id%in%as.character(ensembl.partners[,"NMD.targets"]),"NMD_nonNMD_ratio"] <- best.NMD.ratio
    #ensembl.NMD.features.filt <- ensembl.NMD.features.filt[order(ensembl.NMD.features.filt$ensembl_gene_id),]
  
  }
  # Save our NMD geneset
  ensembl.NMD.features.filt2 <- ensembl.NMD.features.filt[!is.na(ensembl.NMD.features.filt$final_consensus),]
  ensembl.NMD.features.filt2 <- ensembl.NMD.features.filt2[order(ensembl.NMD.features.filt2$gene_symbol),]

  # NMD features plots
  #ensembl.NMD.features.filt3 <- ensembl.NMD.features.filt2[ensembl.NMD.features.filt2$final_consensus=="NMD_target",]
  if (NMD.geneset != "non.NMD.genes") {
    NMD.genesets.res.tmp <- NMD_geneset_info_features(ensembl_NMD_features = ensembl.NMD.features.filt2[ensembl.NMD.features.filt2$final_consensus=="NMD_target",], NMD_geneset = NMD.geneset) 
  } else {
    NMD.genesets.res.tmp <- NMD_geneset_info_features(ensembl_NMD_features = ensembl.NMD.features.filt2[ensembl.NMD.features.filt2$final_consensus=="control",], NMD_geneset = NMD.geneset) 
  }
  #NMD.genesets.res.control <- NMD_geneset_info_features(ensembl_NMD_features = ensembl.NMD.features.filt2[ensembl.NMD.features.filt2$final_consensus=="control",], NMD_geneset = NMD.geneset) 
  NMD.genesets.res.filt <- rbind(NMD.genesets.res.filt,NMD.genesets.res.tmp)

  NMD.geneset <- gsub("\\.","_",NMD.geneset)
  write.table(ensembl.NMD.features.filt, file= paste0(NMD.targets.ensembl.path,NMD.geneset,"_ensembl.txt"), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
  write.table(ensembl.NMD.features.filt2, file= paste0(NMD.targets.ensembl.path,NMD.geneset,"_ensembl_filt.txt"), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
  
  # Plot gene expression of NMD targets and controls
  png(file = paste0(NMD.targets.ensembl.path,NMD.geneset,"_gene_exp_boxplot.png"), width = 4000, height = 3000, res = 300)
  p <- ggplot(data = ensembl.NMD.features.filt2, aes(x=final_consensus, y = log10(gene_exp), color = final_consensus)) +
    geom_boxplot() + ylim(c(0,1.5)) +
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
    ylab("log10(TPM gene exp)") + ggtitle(paste0("Gene expression of --> ",NMD.geneset)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=13, face="bold"),
          axis.title.y = element_text(color="black", size=13, face="bold"))
  print(p)
  dev.off()
}

# Raw NMD genesets

NMD.genesets.res$NMD_geneset <- rownames(NMD.genesets.res)
NMD.genesets.res$name <- paste0(NMD.genesets.res$NMD_geneset,"_",NMD.genesets.res$transcripts_size)

for (NMD.feature in NMD.features) {

  NMD.genesets.res[,paste0(NMD.feature,"_log2FC")] <- round(log2(NMD.genesets.res[,NMD.feature] / NMD.genesets.res[NMD.genesets.res$NMD_geneset=="non.NMD.genes",NMD.feature ]),2)
  
  NMD.genesets.res <- NMD.genesets.res[order(NMD.genesets.res[,NMD.feature]),]
  NMD.genesets.res$name <- factor(NMD.genesets.res$name, levels = NMD.genesets.res$name)
  png(file = gsub("\\[X\\]",NMD.feature,paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_barplot","path_or_filename"])), width = 4000, height = 3000, res = 300)
  p <- ggplot(data = NMD.genesets.res, aes(x=name, y=eval(parse(text=NMD.feature)), group = name, color = name, las = 2)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          axis.text.x = element_text(color="black", size=12, angle = 45, hjust = 1)) +
    geom_text(aes(label=eval(parse(text=NMD.feature))), vjust=1.6, color="white", size=3.5) + ylab(NMD.feature)
  print(p)
  dev.off()
  
}

write.table(NMD.genesets.res, file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_res","path_or_filename"]), quote=FALSE, sep='\t',row.names=TRUE, col.names = TRUE)

# Filtered NMD genesets by NMD features

NMD.genesets.res.filt$NMD_geneset <- rownames(NMD.genesets.res.filt)
NMD.genesets.res.filt$name <- paste0(NMD.genesets.res.filt$NMD_geneset,"_",NMD.genesets.res.filt$transcripts_size)

for (NMD.feature in NMD.features) {

  NMD.genesets.res.filt[,paste0(NMD.feature,"_log2FC")] <- round(log2(NMD.genesets.res.filt[,NMD.feature] / NMD.genesets.res.filt[NMD.genesets.res.filt$NMD_geneset=="non.NMD.genes",NMD.feature ]),2)
  
  NMD.genesets.res.filt <- NMD.genesets.res.filt[order(NMD.genesets.res.filt[,NMD.feature]),]
  NMD.genesets.res.filt$name <- factor(NMD.genesets.res.filt$name, levels = NMD.genesets.res.filt$name)
  png(file = gsub("\\[X\\]",NMD.feature,paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_filt_barplot","path_or_filename"])), width = 4000, height = 3000, res = 300)
  p <- ggplot(data = NMD.genesets.res.filt, aes(x=name, y=eval(parse(text=NMD.feature)), group = name, color = name, las = 2)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          axis.text.x = element_text(color="black", size=12, angle = 45, hjust = 1)) +
    geom_text(aes(label=eval(parse(text=NMD.feature))), vjust=1.6, color="white", size=3.5) + ylab(NMD.feature)
  print(p)
  dev.off()
  
}

write.table(NMD.genesets.res.filt, file = paste0(NMD.targets.res,paths[paths$folder_or_object=="NMD_genesets_res_filt","path_or_filename"]), quote=FALSE, sep='\t',row.names=TRUE, col.names = TRUE)

## Karousis test

Karousis.ensembl.gene.filt <- Karousis.ensembl.gene[grep("ENST",Karousis.ensembl.gene$NMD_isoforms),]
Karousis.ensembl.gene.filt <- Karousis.ensembl.gene.filt[grep("ENST",Karousis.ensembl.gene.filt$NMD_target),]
table(Karousis.ensembl.gene.filt$NMD_isoforms==Karousis.ensembl.gene.filt$NMD_target)
table(Karousis.ensembl.gene.filt$non_NMD_isoforms==Karousis.ensembl.gene.filt$NMD_target)


Karousis.ensembl.gene.filt <- Karousis.ensembl.gene[grep("ENST",Karousis.ensembl.gene$non_NMD_isoforms),]
Karousis.ensembl.gene.filt <- Karousis.ensembl.gene.filt[grep("ENST",Karousis.ensembl.gene.filt$control),]
table(Karousis.ensembl.gene.filt$non_NMD_isoforms==Karousis.ensembl.gene.filt$control)



Karousis.ensembl.gene.filt[30:40,]
ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id %in% "ENSG00000164054",c("ensembl_transcript_id","comment","NMD_event_type","gene_exp","UTR3_GC_content")]

# Expression of Karousis targets -> 75% is less than 1
# That's why most of my targets do not coincide
Karousis.ensembl.gene.filt <- Karousis.ensembl.gene[grep("ENST",Karousis.ensembl.gene$NMD_isoforms),]
table(ensembl.NMD.features[ensembl.NMD.features$ensembl_transcript_id%in%Karousis.ensembl.gene.filt$NMD_isoforms,"gene_exp"] < 1)
# If we remove those...
Karousis.filt <- ensembl.NMD.features[ensembl.NMD.features$ensembl_transcript_id%in%Karousis.ensembl.gene.filt$NMD_isoforms,]
Karousis.filt <- Karousis.filt[Karousis.filt$gene_exp >= 1,]

Karousis.filt <- Karousis.ensembl.gene.filt[Karousis.ensembl.gene.filt$NMD_isoforms%in%Karousis.filt$ensembl_transcript_id,]
Karousis.filt <- Karousis.filt[grep("ENST",Karousis.filt$NMD_target),]
# It's 50% haha
table(Karousis.filt$NMD_isoforms==Karousis.filt$NMD_target)

Karousis.filt[60:70,]
ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id %in% "ENSG00000143774",c("ensembl_transcript_id","comment","NMD_event_type","gene_exp","UTR3_GC_content")]



