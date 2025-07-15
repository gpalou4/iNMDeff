rm(list=ls())

################################ FUNCTIONS ################################
################################ FUNCTIONS ################################
################################ FUNCTIONS ################################

ORFs_num_and_length <- function(ORFs_lengths, min_bp, min_ORFs) {
  
  ORFs_lengths_list <- strsplit(ORFs_lengths, ",")
  ORFs_bool <- lapply(ORFs_lengths_list, function(x) {
    # print(sum(as.numeric(x)))
    x <- as.numeric(x)
    if (length(x) == 0) {
      return(FALSE)
    } else if (is.na(x)) {
      return(FALSE)
    } else if ( sum(x >= min_bp) >= min_ORFs ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  return(unlist(ORFs_bool))
}

NMD_geneset_info_features <- function(ensembl_NMD_features, NMD_geneset, all_NMD_genes = all_NMD_genes, non_NMD_genes = non_NMD_genes) {

  # Keep info
  n <- nrow(ensembl_NMD_features)
  n_uORFs <- as.numeric(table((as.numeric(ensembl_NMD_features$uORFs) > 0))[2])
  n_genes <- length(unique(ensembl_NMD_features$ensembl_gene_id))
  NMD_features_df <- data.frame(NMD_features = c(NMD_features,"fisher_pval","fisher_OR"), value = NA, row_names = 1)
  colnames(NMD_features_df)[1] <- NMD_geneset

  NMD_features_df["transcripts_size",] <- n
  NMD_features_df["genes_size",] <- n_genes
  NMD_features_df["uORFs",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs" & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD_features_df["uORFs_1_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_1_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD_features_df["uORFs_2_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_2_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD_features_df["UTR3_EJC",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR") / n,2)
  NMD_features_df["UTR3_EJC_50nt",] <- round(length(which(ensembl_NMD_features$stop_codon_to_3UTR_splice_site >= 50 & ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR")) / n,2)
  #NMD_features_df["uORFs_UTR3_intron",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs | Intronic-spliced 3UTR") / n,2)
  NMD_features_df["non_NMD",] <- round(sum(ensembl_NMD_features$NMD_event_type == "") / n,2)
  NMD_features_df["UTR3_GC_content",] <- median(as.numeric(ensembl_NMD_features$UTR3_GC_content), na_rm = TRUE)
  NMD_features_df["uORFs_num_translated",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_translated), na_rm = TRUE),2)
  NMD_features_df["uORFs_num_trans_reinit",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_trans_reinit), na_rm = TRUE),2)
  NMD_features_df["penultime_codon_GCG",] <- round(length(grep("GCG",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Ala
  #NMD_genesets_res[filt,"penultime_codon_GTG_control"] <- round(length(grep("GTG",ensembl_NMD_features_filt$uORF_penultimate_codon)) / n,6) # Val
  NMD_features_df["penultime_codon_ACA_control",] <- round(length(grep("ACA",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Thr
  NMD_features_df["penultime_codon_AGG",] <- round(length(grep("AGG",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Arg
  #NMD_genesets_res[filt,"penultime_codon_AAG_control"] <- round(length(grep("AAG",ensembl_NMD_features_filt$uORF_penultimate_codon)) / n,6) # Lys
  NMD_features_df["penultime_codon_CCA_control",] <- round(length(grep("CCA",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Pro
  NMD_features_df["penultime_codon_CTG",] <- round(length(grep("CTG",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Leu
  #NMD_genesets_res[filt,"penultime_codon_CAT_control"] <- round(length(grep("CAT",ensembl_NMD_features_filt$uORF_penultimate_codon)) / n,6) # His
  NMD_features_df["penultime_codon_TAC_control",] <- round(length(grep("TAC",ensembl_NMD_features$uORF_penultimate_codon)) / n_uORFs,6) # Tyr
  
  # Fisher test
  ensembl_NMD_features[,"NMD_geneset"] <- NA
  ensembl_NMD_features[ensembl_NMD_features$ensembl_gene_id%in%eval(parse(text=paste0(NMD_geneset))),"NMD_geneset"] <- NMD_geneset
  if (NMD_geneset == "non_NMD_genes") {
    ensembl_NMD_features[ensembl_NMD_features$ensembl_gene_id%in%all_NMD_genes,"NMD_geneset"] <- "controls"
  } else {
    ensembl_NMD_features[ensembl_NMD_features$ensembl_gene_id%in%non_NMD_genes,"NMD_geneset"] <- "controls"
  }

  print(table(ensembl_NMD_features$NMD_geneset))
  # NMD_feature <- names(table(ensembl_NMD_features$NMD_event_type))[4]
  NMD_noNMDfeat <- length(which(ensembl_NMD_features$NMD_event_type=="" & ensembl_NMD_features$NMD_geneset == NMD_geneset))
  NMD_NMDfeat <- length(which(ensembl_NMD_features$NMD_event_type!="" & ensembl_NMD_features$NMD_geneset == NMD_geneset))
  noNMD_noNMDfeat <- length(which(ensembl_NMD_features$NMD_event_type=="" & ensembl_NMD_features$NMD_geneset == "controls"))
  noNMD_NMDfeat <- length(which(ensembl_NMD_features$NMD_event_type!="" & ensembl_NMD_features$NMD_geneset == "controls"))
  fisher_table <- data.frame(NMD_feature = c(NMD_NMDfeat,noNMD_NMDfeat), nonNMD_feature = c(NMD_noNMDfeat,noNMD_noNMDfeat))
  rownames(fisher_table) <- c(NMD_geneset,"controls")
  fisher_res <- fisher_test(fisher_table)
  print(fisher_table)
  print(fisher_res)
  # chisq_test(fisher_table)$expected
  NMD_features_df["fisher_pval",]  <- fisher_res$p_value
  NMD_features_df["fisher_OR",]  <- fisher_res$estimate
  return(t(NMD_features_df))

}

################################ LIBRARIES ################################
################################ LIBRARIES ################################
################################ LIBRARIES ################################

library("matrixStats")
library("ggplot2")
library("dplyr")
# Read GTF
library("rtracklayer")
# Ven Diagrams
library("VennDiagram")
library("UpSetR")
library("tidyr")

################################ SCRIPT ################################
################################ SCRIPT ################################
################################ SCRIPT ################################

# 1) Data
paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/list_NMD_transcripts/list_ensembl_transcripts_NMD_features_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# args <- commandArgs(trailingOnly=TRUE)
# paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
conversor_tables_path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
# NMD targets
NMD_targets_gene_symbols_updated <- paths[paths$folder_or_object=="NMD_targets_gene_symbols_updated","path_or_filename"]
NMD_targets_raw <- paths[paths$folder_or_object=="NMD_targets_raw","path_or_filename"]
NMD_targets_ensembl_path <- paths[paths$folder_or_object=="NMD_targets_ensembl","path_or_filename"]
# Results
NMD_targets_res <- paths[paths$folder_or_object=="NMD_genesets_res_path","path_or_filename"]
TCGA_names_path <- paths[paths$folder_or_object=="TCGA_names_path","path_or_filename"]
RNAseq_path <- paths[paths$folder_or_object=="RNAseq_TPM_path","path_or_filename"]

# 1.1) ENSEMBL with NMD features
ensembl_NMD_features <- read.table(file = paste0(conversor_tables_path,paths[paths$folder_or_object=="ensembl_NMD_features","path_or_filename"]), 
                                    header = TRUE, sep = "\t", colClasses = "character")
# 1.2) ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
# Remove versions from IDs
ensembl_v88_gtf[,c("gene_id","transcript_id")] <- sapply(ensembl_v88_gtf[,c("gene_id","transcript_id")], function(col) { 
  gsub("(.*)\\..*","\\1", col)
}) 

# 1.3) Conversor Table for Gene Symbols - ENSEMBL gene ID

ensembl_v88_gene_transcript_genesymbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt", 
                                        header = TRUE, sep = "\t", colClasses = "character")

# X) Conversor Affymetrix human genome GeneChip U133 plus array - ensembl table
# ensembl_affyU133plus_conversion <- read.table(file = paste0(conversor_tables_path,paths[paths$folder_or_object=="ensemblToU133Plus2","path_or_filename"]), 
#                                            header = FALSE, sep = "\t")
# ensembl_affyU133plus_conversion[,1] <- sub("\\__*","", ensembl_affyU133plus_conversion[,1])
# colnames(ensembl_affyU133plus_conversion) <- c("ensembl","AFFY_U133_PLUS")

# 1.4) NMD targets
# Gene Symbols and ENSEMBL genes IDs will be transformed to ENSEMBL transcripts IDs using the ENSEMBL v88 GTF file

# 1.4.1) Tani
# Updated gene symbols for UPF1

# Group A
Tani_UPF1_dir_unconf_gene_symbols <- read.table(file = paste0(NMD_targets_gene_symbols_updated,"Tani_UPF1_dir_unconf_gene_symbols_updated.txt"), 
                                                   header = FALSE, sep = "\t", colClasses = "vector")$V1
# Group B
Tani_UPF1_indir_conf_gene_symbols <- read.table(file = paste(NMD_targets_gene_symbols_updated,"Tani_UPF1_ind_conf_gene_symbols_updated.txt", sep = ""), 
                                                   header = FALSE, sep = "\t", colClasses = "vector")$V1
# Group C
Tani_UPF1_dir_conf_gene_symbols <- read.table(file = paste(NMD_targets_gene_symbols_updated, "Tani_UPF1_dir_conf_gene_symbols_updated.txt", sep = ""), 
                                                 header = FALSE, sep = "\t", colClasses = "vector")$V1
# All merged
Tani_UPF1_all_gene_symbols <- read.table(file = paste(NMD_targets_gene_symbols_updated, "Tani_UPF1_all_gene_symbols_updated.txt", sep = ""), 
                                                 header = FALSE, sep = "\t", colClasses = "vector")$V1
# Convert Gene Symbols to ENSEMBL genes ID
Tani_UPF1_all_gene_symbols[!Tani_UPF1_all_gene_symbols%in%ensembl_v88_gene_transcript_genesymbol$gene_name]
Tani_UPF1_all_ensembl_gene <- unique(ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name%in%Tani_UPF1_all_gene_symbols,"gene_id"])
Tani_UPF1_dir_unconf_ensembl_gene <- unique(ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name%in%Tani_UPF1_dir_unconf_gene_symbols,"gene_id"])
Tani_UPF1_dir_conf_ensembl_gene <- unique(ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name%in%Tani_UPF1_dir_conf_gene_symbols,"gene_id"])
Tani_UPF1_indir_conf_ensembl_gene <- unique(ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name%in%Tani_UPF1_indir_conf_gene_symbols,"gene_id"])
# 1.4.2) Schmidt
# Updated gene symbols for SMG6
Schmidt_SMG6_gene_symbols <- read.table(file = paste(NMD_targets_gene_symbols_updated, "Schmidt_SMG6_gene_symbol_updated.txt", sep = ""), 
                                                 header = TRUE, sep = "\t", colClasses = "vector")
# Convert Schmidt Gene Symbols to ENSEMBL gene ID
Schmidt_SMG6_gene_symbols[!Schmidt_SMG6_gene_symbols$gene_symbol%in%ensembl_v88_gene_transcript_genesymbol$gene_name,]
Schmidt_SMG6_ensembl_gene <- ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name%in%Schmidt_SMG6_gene_symbols$gene_symbol,"gene_id"]
# 1.4.3) Colombo
# ENSEMBL genes IDs for SMG6, SMG7 and UPF1
Colombo_SMG6_ensembl_gene <- read.table(file = paste(NMD_targets_raw, "Colombo_SMG6_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
Colombo_SMG7_ensembl_gene <- read.table(file = paste(NMD_targets_raw, "Colombo_SMG7_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
Colombo_UPF1_ensembl_gene <- read.table(file = paste(NMD_targets_raw, "Colombo_UPF1_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
# 1.4.4) Karousis
# ENSEMBL genes IDs for NMD global
Karousis_ensembl_gene <- read.table(file = paste(NMD_targets_raw, "Karousis_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
# 1.4.5) Courtney
# ENSEMBL genes IDs for NMD global
Courtney_ensembl_gene <- read.table(file = paste(NMD_targets_raw, "Courtney_ensembl_gene_symbol.txt", sep = ""), 
                                                    header = TRUE, sep = "\t", colClasses = "vector")
# 1.4.6) NMD targets from ensembl gtf
ensembl_v88_gtf_NMD_targets <- unique(ensembl_v88_gtf[ensembl_v88_gtf$transcript_type%in%"nonsense_mediated_decay","gene_id"])
ensembl_v88_gtf_NMD_targets <- gsub("(.*)\\..*","\\1",ensembl_v88_gtf_NMD_targets)
ensembl_v88_gtf_NMD_targets <- gsub("PARY","",ensembl_v88_gtf_NMD_targets)

# 1.5) Pantissue TPM RNAseq
# 1.5.1) TCGA
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt",
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 1.5.2) GTEx
GTEx_RNAseq_TPM_pantissue <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct",
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)
print(dim(GTEx_RNAseq_TPM_pantissue))
# Fix Genes ID (Remove genes duplicated in Y chr)
GTEx_RNAseq_TPM_pantissue <- GTEx_RNAseq_TPM_pantissue[-grep("PAR",GTEx_RNAseq_TPM_pantissue$transcript_id),]
GTEx_RNAseq_TPM_pantissue$transcript_id <- gsub("(\\..*)","",GTEx_RNAseq_TPM_pantissue$transcript_id)
keep <- c(1,grep("GTEX",colnames(GTEx_RNAseq_TPM_pantissue)))
GTEx_RNAseq_TPM_pantissue <- GTEx_RNAseq_TPM_pantissue[,keep]
# 1.6) Transcripts TSS and TES
transcript_start_end_coords <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v92_transcript_start_end_coords.txt"
                                          , header = TRUE, sep = "\t")
# 1.7) MANE
ensembl_transcripts_MANE_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_MANE_v107.txt"
ensembl_transcripts_MANE <- read.table(file = ensembl_transcripts_MANE_path, header = TRUE)
# First the MANE Select if possible, if not, then randomly
ensembl_transcripts_MANE$MANE_main_isoform <- "yes"
# 1.8) Half-life
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome.csv"
half_life_1 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_2.csv"
half_life_2 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_bcells.csv"
half_life_3 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_PCA.csv"
half_life_4 <- read.table(file = file, header = TRUE, sep = ",")
# 1.9) Protein coding transcripts
ensembl_v88_coding_transcripts <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_coding_transcripts.txt",
                                  header = FALSE, sep = "\t",colClasses = "vector") 
# 1.10) Transcript length
ensembl_v88_transcripts_length <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_transcripts_length.txt",
                                  header = TRUE, sep = "\t") 

# 2) Create NMD genesets

# Remove Tani_UPF1_dir_unconf_ensembl_gene from all the NMD genesets
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Colombo_SMG7_ensembl_gene$gene)
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Colombo_SMG6_ensembl_gene$gene)
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Colombo_UPF1_ensembl_gene$gene)
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Karousis_ensembl_gene$gene_id)
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Schmidt_SMG6_ensembl_gene)
table(Tani_UPF1_dir_unconf_ensembl_gene%in%Courtney_ensembl_gene$Ensembl_Gene_ID)

Tani_UPF1_all_ensembl_gene <- Tani_UPF1_all_ensembl_gene[!Tani_UPF1_all_ensembl_gene%in%Tani_UPF1_dir_unconf_ensembl_gene]
Colombo_SMG7_ensembl_gene <- Colombo_SMG7_ensembl_gene[!Colombo_SMG7_ensembl_gene$gene%in%Tani_UPF1_dir_unconf_ensembl_gene,]
Colombo_SMG6_ensembl_gene <- Colombo_SMG6_ensembl_gene[!Colombo_SMG6_ensembl_gene$gene%in%Tani_UPF1_dir_unconf_ensembl_gene,]
Colombo_UPF1_ensembl_gene <- Colombo_UPF1_ensembl_gene[!Colombo_UPF1_ensembl_gene$gene%in%Tani_UPF1_dir_unconf_ensembl_gene,]
Karousis_ensembl_gene <- Karousis_ensembl_gene[!Karousis_ensembl_gene$gene_id%in%Tani_UPF1_dir_unconf_ensembl_gene,]
Schmidt_SMG6_ensembl_gene <- Schmidt_SMG6_ensembl_gene[!Schmidt_SMG6_ensembl_gene%in%Tani_UPF1_dir_unconf_ensembl_gene]
Courtney_ensembl_gene <- Courtney_ensembl_gene[!Courtney_ensembl_gene$Ensembl_Gene_ID%in%Tani_UPF1_dir_unconf_ensembl_gene,]

# 2.1) SMG6 (genes shared in Colombo + Schmidt)

SMG6 <- unique(Colombo_SMG6_ensembl_gene[Colombo_SMG6_ensembl_gene$gene%in%Schmidt_SMG6_ensembl_gene,"gene"])

# 2.2) SMG7 (genes in Colombo only)

SMG7 <- unique(Colombo_SMG7_ensembl_gene$gene)

# 2.3) NMD_tani (Tani UPF1, EJC indep)

NMD_Tani <- unique(Tani_UPF1_all_ensembl_gene)

# 2.4) NMD_colombo (Colombo UPF1 + SMG6 + SMG7, EJC indep/dep)

NMD_Colombo <- unique(c(Colombo_SMG6_ensembl_gene$gene, Colombo_SMG7_ensembl_gene$gene, Colombo_UPF1_ensembl_gene$gene))

# 2.5) NMD_karousis (Karousis UPF1 + SMG6 + SMG7, EJC indep/dep)
# We will accept genes that appears >= 2 times in Karousis as "shared in 2 papers"
# Reason? I'll have to find one LMAO
NMD_Karousis <- Karousis_ensembl_gene$gene_id

# 2.6) NMD_courtney (Karousis UPF1, EJC indep/dep)
# We will accept genes that appears >= 2 times in Courtney as "shared in 2 papers"
# These genes have 1 endogenous PTC in >= 2 transcripts or >= 2 PTCs in 1 transcript, which fives higher confidence
# (in reality, these genes improve NMDeff median tissues in my tests)
NMD_Courtney <- Courtney_ensembl_gene$Ensembl_Gene_ID

# 2.7) NMD global (all genes except ensembl)
NMD_global <- unique(c(Tani_UPF1_all_ensembl_gene, Colombo_SMG6_ensembl_gene$gene, Colombo_SMG7_ensembl_gene$gene, Colombo_UPF1_ensembl_gene$gene, 
Schmidt_SMG6_ensembl_gene, Karousis_ensembl_gene$gene_id, Courtney_ensembl_gene$Ensembl_Gene_ID))

# 2.8) NMD_global_all_shared (all genes shared in all papers except Schmidt)
all_NMD_genesets <- list(NMD_Tani = NMD_Tani, NMD_Colombo = NMD_Colombo, NMD_Karousis = NMD_Karousis, NMD_Courtney = NMD_Courtney)
NMD_global_all_shared <- Reduce(intersect, all_NMD_genesets)

# 2.9) NMD_global_2_shared (genes shared in >= 2 papers)
all_NMD_genesets_vector <- unlist(all_NMD_genesets)
#NMD_global_2_shared <- unique(all_NMD_genesets_vector [duplicated(all_NMD_genesets_vector )])
NMD_genes_frequency <- data.frame(table(all_NMD_genesets_vector))
colnames(NMD_genes_frequency)[1] <- "ensemb_gene_id" 
NMD_genes_frequency <- NMD_genes_frequency[order(NMD_genes_frequency$Freq, decreasing = TRUE),]
NMD_global_2_shared <- as.character(NMD_genes_frequency[NMD_genes_frequency$Freq >= 2,"ensemb_gene_id"])

# 2.10) NMD ensembl
NMD_ensembl <- ensembl_v88_gtf_NMD_targets

# 2.11) For the Negative Control

all_NMD_genes <- unique(c(NMD_global,Tani_UPF1_dir_unconf_ensembl_gene,NMD_ensembl))
non_NMD_genes <- unique(ensembl_NMD_features[!ensembl_NMD_features$ensembl_gene_id%in%all_NMD_genes,"ensembl_gene_id"])

# 2.12) UPF3B independent genes symbols #####

# UPF3B_affymetrix <- read.table(file = paste(NMD_targets_raw, "UPF3B_affymetrix.txt", sep = ""),
                               # header = FALSE, sep = "\t", colClasses = "vector")$V1
# Convert to ensembl
# UPF3B_ensembl <- ensembl_affyU133plus_conversion[ensembl_affyU133plus_conversion$AFFY_U133_PLUS%in%UPF3B_affymetrix,"ensembl"]
# Convert to Gene Symbols
# UPF3B_gene_symbols <- unique(ensembl_gene_symbol_conversion_updated_updated[ensembl_gene_symbol_conversion_updated_updated$FeatureID%in%UPF3B_ensembl,"Approved_symbol"])

# 3) Ven Diagrams of shared NMD genesets

# library(RColorBrewer)
# myCol <- brewer_pal(4, "Pastel2")

# venn_diagram(
#         x = all_NMD_genesets,
#         TACegory_names = names(all_NMD_genesets),
#         filename = paste0(NMD_targets_res,paths[paths$folder_or_object=="NMD_genesets_venn_diagram","path_or_filename"]),
#         output=TRUE,
        
#         # Output features
#         imagetype="png" ,
#         height = 1000 , 
#         width = 1000 , 
#         resolution = 300,
#         compression = "lzw",
        
#         # Circles
#         lwd = 2,
#         lty = 'blank',
#         fill = myCol,
        
#         # Numbers
#         cex = .8,
#         fontface = "bold",
#         fontfamily = "sans",

#         # Set names
#         TAC_cex = 0.6,
#         TAC_fontface = "bold",
#         TAC_fontfamily = "sans",
# )

# # UPF1 intersect
# png(file = paste0(NMD_targets_res,paths[paths$folder_or_object=="NMD_genesets_intersection_UPF1","path_or_filename"]), width = 4000, height = 3000, res = 300)
# NMD_genesets_list = list(Tani_UPF1_dir_conf = Tani_UPF1_dir_conf_ensembl_gene, Tani_UPF1_indir_conf = Tani_UPF1_indir_conf_ensembl_gene,
#                           Colombo_UPF1 = Colombo_UPF1_ensembl_gene$gene, NMD_Karousis = NMD_Karousis, NMD_Courtney = NMD_Courtney, NMD_ensembl = NMD_ensembl)
# p <- upset(fromList(NMD_genesets_list), nsets = length(NMD_genesets_list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD_genesets_list), keep.order = TRUE)
# print(p)
# dev.off()

# # SMG6/7 intersect
# png(file = paste0(NMD_targets_res,paths[paths$folder_or_object=="NMD_genesets_intersection_SMG6_7","path_or_filename"]), width = 4000, height = 3000, res = 300)
# NMD_genesets_list = list(Schmidt_SMG6 = Schmidt_SMG6_ensembl_gene, Colombo_SMG6 = Colombo_SMG6_ensembl_gene$gene, Colombo_SMG7 = Colombo_SMG7_ensembl_gene$gene,
#                           Colombo_UPF1 = Colombo_UPF1_ensembl_gene$gene, NMD_Karousis = NMD_Karousis, NMD_Courtney = NMD_Courtney, NMD_ensembl = NMD_ensembl)
# p <- upset(fromList(NMD_genesets_list), nsets = length(NMD_genesets_list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD_genesets_list), keep.order = TRUE)
# print(p)
# dev.off()

# # All NMD genesets separately
# png(file = paste0(NMD_targets_res,paths[paths$folder_or_object=="NMD_genesets_intersection_all","path_or_filename"]), width = 4000, height = 3000, res = 300)
# NMD_genesets_list = list(Tani_dir_conf = Tani_UPF1_dir_conf_ensembl_gene, Tani_indir_conf = Tani_UPF1_indir_conf_ensembl_gene,
#                           Schmidt_SMG6 = Schmidt_SMG6_ensembl_gene, Colombo_SMG6 = Colombo_SMG6_ensembl_gene$gene, Colombo_SMG7 = Colombo_SMG7_ensembl_gene$gene,
#                           Colombo_UPF1 = Colombo_UPF1_ensembl_gene$gene, NMD_Karousis = NMD_Karousis, NMD_Courtney = NMD_Courtney, NMD_ensembl = NMD_ensembl)
# p <- upset(fromList(NMD_genesets_list), nsets = length(NMD_genesets_list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD_genesets_list), keep.order = TRUE)
# print(p)
# dev.off()

# # NMD genesets
# png(file = paste0(NMD_targets_res,paths[paths$folder_or_object=="NMD_genesets_intersection_final","path_or_filename"]), width = 4000, height = 3000, res = 300)
# NMD_genesets_list = list(NMD_tani = NMD_Tani, NMD_Colombo = NMD_Colombo, NMD_Karousis = NMD_Karousis, 
#                         NMD_Courtney = NMD_Courtney, NMD_ensembl = NMD_ensembl)
# p <- upset(fromList(NMD_genesets_list), nsets = length(NMD_genesets_list), order.by = c("freq"), empty.intersections = "on", sets = names(NMD_genesets_list), keep.order = TRUE)
# print(p)
# dev.off()

################################################################################################################################################

# 3) Add additional transcript information 

# 3.1) Add gene expression for all ENSEMBL transcripts IDs
# TCGA
TCGA_RNAseq_TPM_pancancer_filt <- TCGA_RNAseq_TPM_pancancer %>%
                                  filter( rownames(TCGA_RNAseq_TPM_pancancer) %in% ensembl_NMD_features$ensembl_transcript_id)
ensembl_transcripts <- rownames(TCGA_RNAseq_TPM_pancancer_filt)
TCGA_pancancer_gene_exp <- rowMedians(as.matrix(TCGA_RNAseq_TPM_pancancer_filt),na.rm = TRUE)
TCGA_pancancer_gene_exp_df <- data.frame(ensembl_transcript_id = ensembl_transcripts, TCGA_pancancer_gene_exp = TCGA_pancancer_gene_exp)
ensembl_NMD_features <- merge(ensembl_NMD_features, TCGA_pancancer_gene_exp_df, by = "ensembl_transcript_id", all.x = TRUE)
# GTEx
GTEx_RNAseq_TPM_patissue_filt <- GTEx_RNAseq_TPM_pantissue %>%
                                  filter( GTEx_RNAseq_TPM_pantissue$transcript_id %in% ensembl_NMD_features$ensembl_transcript_id)
ensembl_transcripts <- GTEx_RNAseq_TPM_patissue_filt$transcript_id
GTEx_pantissue_gene_exp <- rowMedians(as.matrix(GTEx_RNAseq_TPM_patissue_filt[,-1]),na.rm = TRUE)
GTEx_pantissue_gene_exp_df <- data.frame(ensembl_transcript_id = ensembl_transcripts, GTEx_pantissue_gene_exp = GTEx_pantissue_gene_exp)
ensembl_NMD_features <- merge(ensembl_NMD_features, GTEx_pantissue_gene_exp_df, by = "ensembl_transcript_id", all.x = TRUE)

# 3.2) Check uORFs number and length
ORFs_num_len_bool <- ORFs_num_and_length(ORFs_lengths = ensembl_NMD_features$uORFs_lengths, min_bp = 30, min_ORFs = 1)
ensembl_NMD_features[,"uORFs_1_30_bp"] <- ORFs_num_len_bool
ORFs_num_len_bool <- ORFs_num_and_length(ORFs_lengths = ensembl_NMD_features$uORFs_lengths, min_bp = 30, min_ORFs = 2)
ensembl_NMD_features[,"uORFs_2_30_bp"] <- ORFs_num_len_bool

# 3.3) Found on Karousis dataset
ensembl_NMD_features <- ensembl_NMD_features %>%
                        mutate(Karousis_paper = case_when( 
                          ensembl_transcript_id %in% Karousis_ensembl_gene$NMD_isoforms ~ "NMD_target",
                          ensembl_transcript_id %in% Karousis_ensembl_gene$non_NMD_isoforms ~ "NMD_control"
                        ))

# 3.7) MANE transcript
ensembl_NMD_features <- ensembl_NMD_features %>%
                        mutate(MANE_Select = case_when( 
                          ensembl_transcript_id %in% ensembl_transcripts_MANE$ensembl_transcript_id ~ TRUE,
                        ))
ensembl_NMD_features$MANE_Select <- ifelse(is.na(ensembl_NMD_features$MANE_Select), FALSE, TRUE)

# 3.8) Half-Life
colnames(half_life_4)[3] <- "half_life_PC1"
ensembl_NMD_features <- merge(ensembl_NMD_features, half_life_4[,c("Ensembl.Gene.Id","half_life_PC1")], by.x = "ensembl_gene_id", by.y = "Ensembl.Gene.Id", all.x = TRUE)

# 3.9) GENCODE ENSEMBL information
# Tag (Including NMD decay), transcript_support_level, transcript type and Add Start and Stop codon positions (?)
ensembl_v88_gtf_filt <- ensembl_v88_gtf[,c("transcript_id","tag","transcript_support_level","transcript_type")]
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[!duplicated(ensembl_v88_gtf_filt),]
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[!is.na(ensembl_v88_gtf_filt$transcript_id),]
ensembl_NMD_features <- merge(ensembl_NMD_features, ensembl_v88_gtf_filt, by.x = "ensembl_transcript_id" , by.y = "transcript_id", all.x = TRUE)
# TSS and TES
ensembl_NMD_features <- merge(ensembl_NMD_features, transcript_start_end_coords[,c("ensembl_transcript_id","transcription_start_site","transcription_end_site")], by.x = "ensembl_transcript_id" , by.y = "ensembl_transcript_id", all.x = TRUE)

# 3.10) Transcript length
colnames(ensembl_v88_transcripts_length)[2] <- "transcript_length" 
ensembl_NMD_features <- merge(ensembl_NMD_features, ensembl_v88_transcripts_length, by.x = "ensembl_transcript_id" , by.y = "transcript_id", all.x = TRUE)

# 4) Filter NMD genesets using our NMD features
# Obtain potential NMD targets and NMD controls for all NMD genesets

ensembl_NMD_features$uORFs <- as.numeric(ensembl_NMD_features$uORFs)
ensembl_NMD_features$uORFs <- ifelse(is.na(ensembl_NMD_features$uORFs),0,ensembl_NMD_features$uORFs)

# Save ENSEMBL 
#write.table(ensembl_NMD_features, file = paste0(NMD_targets_ensembl_path,"ensembl_NMD_features_all.txt"), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
ensembl_NMD_features <- read.table(file = paste0(NMD_targets_ensembl_path,"ensembl_NMD_features_all.txt"),
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# NMD gene criteria
# 1) Is in one of our papers
# 2) Have more than one transcript

NMD_gene_sets_names <- c("all_NMD_genes","non_NMD_genes","SMG6","SMG7","NMD_global", "NMD_ensembl",
                          "NMD_global_all_shared","NMD_global_2_shared","NMD_Karousis","NMD_Courtney","NMD_Tani","NMD_Colombo")

for (NMD_geneset in NMD_gene_sets_names) {
  
  print(paste0("........................ NMD geneset ---> ", NMD_geneset, "........................"))
  # Filter for NMD geneset
  ensembl_NMD_features_filt <- ensembl_NMD_features[ensembl_NMD_features$ensembl_gene_id%in%eval(parse(text=paste0(NMD_geneset))),]
  # Fisher Tests
  # NMD_genesets_res_tmp <- NMD_geneset_info_features(ensembl_NMD_features = ensembl_NMD_features_filt, NMD_geneset = NMD_geneset) 
  # NMD_genesets_res <- rbind(NMD_genesets_res,NMD_genesets_res_tmp)

  # NMD target criteria
  # 1) Have a start and normal stop codon
  # 2) At least one NMD triggering feature 
  # 3) Have not 0 expression
  # EJC at 3'UTR with minimum 50nt from normal stop-codon or >= 2 uORF at least 30 bp each)
  NMD_targets <- ensembl_NMD_features_filt %>%
                              filter( ( NMD_event_type %in% c(" | Intronic-spliced 3UTR") 
                                      & stop_codon_to_3UTR_splice_site >= 50) | ( NMD_event_type %in% c(" | >=2 uORFs | Intronic-spliced 3UTR"," | >=2 uORFs" ) ) ) %>%
                              filter(!grepl(".*no stop.*|.*no start.*",comment)) %>%
                              filter(TCGA_pancancer_gene_exp > 0 | GTEx_pantissue_gene_exp > 0) %>%
                              pull(ensembl_transcript_id)
  # NMD control criteria
  # 1) Have a start and normal stop codon
  # 2) Do not have any NMD triggering feature
  # 3) Have not 0 expression
  NMD_controls <- ensembl_NMD_features_filt %>%
                              filter(NMD_event_type == "" & uORFs == 0) %>%
                              filter(!grepl(".*no stop.*|.*no start.*",comment)) %>%
                              filter(TCGA_pancancer_gene_exp > 0 | GTEx_pantissue_gene_exp > 0) %>%
                              pull(ensembl_transcript_id) 

  # Create NMD_type
  ensembl_NMD_features_filt$NMD_type <- NA
  ensembl_NMD_features_filt[ensembl_NMD_features_filt$ensembl_transcript_id %in% NMD_targets,"NMD_type"] <- "NMD_target"
  ensembl_NMD_features_filt[ensembl_NMD_features_filt$ensembl_transcript_id %in% NMD_controls,"NMD_type"] <- "NMD_control"

  # Save NMD_geneset
  write.table(ensembl_NMD_features_filt, file= paste0(NMD_targets_ensembl_path,NMD_geneset,"_ensembl.txt"), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
  
  median_quartile <- function(x) {
      out <- quantile(x, probs = c(0.25,0.5,0.75))
      names(out) <- c("ymin","y","ymax")
      return(out)
  }
  # Plot gene expression of NMD targets and controls
  df_long <- ensembl_NMD_features_filt[,c("NMD_type","TCGA_pancancer_gene_exp","GTEx_pantissue_gene_exp")] %>% 
              pivot_longer(cols = c(TCGA_pancancer_gene_exp, GTEx_pantissue_gene_exp), names_to = "dataset", values_to = "gene_exp")
  png(file = paste0(NMD_targets_ensembl_path,NMD_geneset,"_gene_exp_boxplot.png"), width = 4000, height = 3000, res = 300)
  p <- ggplot(df_long, aes(x = NMD_type, y = log10(gene_exp), color = NMD_type)) +
      geom_violin(alpha = 0.5,draw_quantiles = c(0.25, 0.5, 0.75)) +
      stat_summary(fun = median_quartile, geom = 'point') +
      geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
      facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
      theme_classic() +
      geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
      ylab(paste0("log10(TPM)")) + ggtitle(paste0("Pantissue gene expression of --> ",NMD_geneset)) +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(color="black", size=13, face="bold"),
      axis.title.y = element_text(color="black", size=13, face="bold"))
  print(p)
  dev.off()
}




