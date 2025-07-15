rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

VCF_sample_file <- function(TCGA_cancer, TCGA_sample, type) {
  if (type == "germline") {
    VCF_path <- gsub("\\[X\\]",gsub("\\.","-", strsplit(TCGA_cancer,"-")[[1]][2]), VCF_germline_path)
    tmp_paths <- list.dirs(path = paste0(VCF_path,gsub("\\.","-",TCGA_sample)), full.names = TRUE, recursive = TRUE)
    VCF_path_sample_germline <- paste0(tmp_paths[5],"/",strsplit(TCGA_cancer,"-")[[1]][2],"_annotated.HG_multianno.txt.gz")
    if (file.exists(VCF_path_sample_germline)) {
      print(VCF_path_sample_germline)
      VCF_sample <- read.table(file = VCF_path_sample_germline, header = TRUE, sep = "\t")
    } else {
      return(NULL)
    }
  }  else if (type == "somatic") {
    VCF_path <- gsub("\\[X\\]",gsub("\\.","-", strsplit(TCGA_cancer,"-")[[1]][2]), VCF_somatic_path)
    tmp_paths <- list.dirs(path = paste0(VCF_path,gsub("\\.","-",TCGA_sample)), full.names = TRUE, recursive = TRUE)
    VCF_path_sample_snvs <- paste0(tmp_paths[4],"/",TCGA_sample,"_snvs_annotated.HG_multianno.txt.gz")
    VCF_path_sample_indels <- paste0(tmp_paths[4],"/",TCGA_sample,"_indels_annotated.HG_multianno.txt.gz")
    existsSNVsFile <- file.exists(VCF_path_sample_snvs)
    existsIndelsFile <- file.exists(VCF_path_sample_indels)
    VCF_sample_snvs <- NULL
    VCF_sample_indels <- NULL
    if (existsSNVsFile) {
      print(VCF_path_sample_snvs)
      VCF_sample_snvs <- read.table(file = VCF_path_sample_snvs, header = TRUE, sep = "\t")
    } 
    if (existsIndelsFile) {
      print(VCF_path_sample_indels)
      VCF_sample_indels <- read.table(file = VCF_path_sample_indels, header = TRUE, sep = "\t")
      indels_number <- nrow(VCF_sample_indels)
    } 
    if ( !is.null(VCF_sample_snvs) & !is.null(VCF_sample_indels) ) {
      VCF_sample <- rbind(VCF_sample_snvs,VCF_sample_indels)
    } else if ( is.null(VCF_sample_snvs) & !is.null(VCF_sample_indels) ) {
      VCF_sample <- VCF_sample_indels 
    } else if ( !is.null(VCF_sample_snvs) & is.null(VCF_sample_indels) ) {
      VCF_sample <- VCF_sample_snvs
    } else {
      return(NULL)
    }
  }
  # Fix VAF numeric
  VAFs <- VCF_sample$AF
  VAFs <- ifelse(VAFs == "" | VAFs == ".",0,format(VAFs,scientific = TRUE))
  VAFs <- as.numeric(VAFs)
  VCF_sample$AF <- VAFs
  # Fix ENSEMBL Genes IDs
  # Some transcript swill have different Genes IDs (for whatever reason), but bc it's only used on the regression as gene.id it's not problem
  VCF_sample$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",VCF_sample$Gene.refGene)
  return(VCF_sample)
}

VCF_sample_filtering <- function(VCF_sample, mutation, VAF, filter = NULL) {
    # Filter variants by: nonsense (PTCs) and VAF or synonymous (negative control)
    VCF_sample_filt <- VCF_sample[VCF_sample$ExonicFunc.refGene %in% mutation &  VCF_sample$AF <= VAF,]
    # PASS variants
    if (!is.null(filter)) {
      VCF_sample_filt <- VCF_sample_filt[VCF_sample_filt$Otherinfo10 == filter,]
    }
    return(VCF_sample_filt)
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

# conda activate /home/dnaro/.conda/envs/NMD_regression_3
.libPaths( rev( .libPaths() ) )
# Read GTF
library("rtracklayer")
# To read fasta
library("seqinr")
# ORFs
library("ORFik")
# Random strings
library("stringi")
library("stringr")
# somatic VCF strelka VAF calculation
library("ICAMS")
library("dplyr")
library("purrr")

##############################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# This script checks for Multi-Nucleotide Variants (MNV) in nonsense mutations
# That means: another missense co-occuring in the same stop codon as the nonsense mutation
# converting the codon into a non-stop (rescued). Frequency is around 5% of nonsense

# 1) Load Data
# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA-UCS/TCGA-QN-A5NN/PTC_transcripts/"
TCGA_cancer <- args[2] # e.g. TCGA-UCEC
TCGA_sample <- args[3] # e.g. TCGA-FI-A2D6
VCF_type <- args[4] # somatic // germline
variant_type <- args[5] # stopgain // missense
print(args)

print(paste0("TCGA cancer --> ",TCGA_cancer))

conversor_tables_path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
VCF_germline_path <- paths[paths$folder_or_object=="VCF_germline_file_path","path_or_filename"]
VCF_somatic_path <- paths[paths$folder_or_object=="VCF_somatic_file_path","path_or_filename"]
transcripts_fasta_path <- paths[paths$folder_or_object=="transcripts_fasta_path","path_or_filename"]
PTC_transcripts_output_path <- paths[paths$folder_or_object=="PTC_transcripts_output_path","path_or_filename"]

# Output path
output_path <- gsub("\\[X1\\]",TCGA_cancer, PTC_transcripts_output_path)
output_path <- gsub("\\[X2\\]",gsub("\\.","-",TCGA_sample), output_path)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
output_path <- paste0(output_path,paste0(VCF_type,"_MNV_nonsense.txt"))
print(output_path)

if (!file.exists(output_path)) {
#if (file.exists(output_path)) {
  
  # 1.1) ENSEMBL transcripts IDs hg38 GTF
  ensembl_v88_gtf <- rtracklayer::import(paste0(conversor_tables_path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
  ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
  # Remove versions from IDs
  ensembl_v88_gtf[,c("gene_id","transcript_id")] <- sapply(ensembl_v88_gtf[,c("gene_id","transcript_id")], function(col) { 
    gsub("(.*)\\..*","\\1", col)
  })
  # 1.2) Transcripts fasta sequences
  transcripts_fasta_sequences <- read.fasta(file = paste0(transcripts_fasta_path,paths[paths$folder_or_object=="transcripts_fasta","path_or_filename"]))
  # 2) NMD rules
  print(paste0("Starting MNV search for TCGA sample --> ",TCGA_sample))
  print(paste0("Nonsense from ",VCF_type," variants"))
  # 2.1) VCF data
  # Obtain VCF germline data file
  # Check if file exists
  VCF_sample <- NULL
  VCF_sample <- VCF_sample_file(TCGA_sample = TCGA_sample, TCGA_cancer = TCGA_cancer, type = VCF_type)
  if ( is.null(VCF_sample) ) {
      print("VCF sample file not found, skipping...")
  } else {

    VCF_sample_filt <- VCF_sample %>%
        dplyr::filter(Func.refGene == "exonic") %>% data.frame()

    MNV_stopgain <- VCF_sample_filt %>%
        arrange(Chr, Start) %>%  # Sort by chromosome and start position
        group_by(Chr) %>%        # Group by chromosome
        mutate(
            next_start = lead(Start),        # Get the next Start position
            prev_start = lag(Start),         # Get the previous Start position
            distance_next = next_start - Start,  # Distance to the next variant
            distance_prev = Start - prev_start   # Distance to the previous variant
        ) %>%
        filter(
            distance_next %in% c(1, 2) |  # Include distances +1 and +2
            distance_prev %in% c(1, 2)    # Include distances -1 and -2
        ) %>%
        filter(ExonicFunc.refGene == "stopgain") %>%
        ungroup() %>% 
        data.frame()

    if (nrow(MNV_stopgain) != 0) {

      # Check if the 
      MNV_stopgain_final <- c()
      for (pos in unique(MNV_stopgain$Start)) {
          print(pos)
          MNV_pos <- VCF_sample_filt %>%
              filter(Start %in% c(pos,pos-1,pos-2,pos+1,pos+2)) %>%
              mutate(
                  # Extract codons from AAChange.refGene
                  codon_list = str_extract_all(AAChange.refGene, "p\\.\\w(\\d+)")
                  # codons = map(codon_list, ~ unique(as.integer(.x)))  # Extract unique codons
                  # is_same_codon = map_lgl(codons, ~ length(.x) == 1)   # Check if all codons are identical
              ) #%>% filter(is_same_codon) #%>%  # Keep only rows where mutations are in the same codon
          print(MNV_pos[,c("AAChange.refGene","codon_list")])
          # Missense and stopgain affect the same codon
          codons_1 <- unlist(MNV_pos$codon_list[1])
          codons_2 <- unlist(MNV_pos$codon_list[2])
          IsSameCodon <- FALSE
          if (length(codons_1) != length(codons_2)) {
              IsSameCodon <- FALSE
              print(IsSameCodon)
          } else {
              IsSameCodon <- sum(codons_1 == codons_2) >= 1
          }
          IsSameCodon
          if (isTRUE(IsSameCodon)) {
              MNV_stopgain_final <- rbind(MNV_stopgain_final,MNV_pos)
          }
      }

      # 2.2) Filter VCF and check for PTCs (or controls)
      # Filter VCF for variants (stopgain or synonymous)
      mutations <- c("stopgain","frameshift insertion","frameshift deletion")
      #mutations <- c("nonsynonymous SNV","synonymous SNV")
      MNV_stopgain_final_filt <- VCF_sample_filtering(VCF_sample = MNV_stopgain_final, 
                                  mutation = mutations, VAF = 1, filter = "PASS")
      if (!is.null(MNV_stopgain_final_filt)) {
        if (nrow(MNV_stopgain_final_filt) == 0) {
            print("No variants in this sample")
        } else {
            # 2.3) Check NMD rules for each PTC
            MNV_stopgain_final_filt$codon_list <- NULL
            write.table(MNV_stopgain_final_filt, file = output_path, sep = "\t", quote = FALSE,
                        col.names = TRUE, row.names = FALSE)
        }
      }
    }
  }
} else {
  print("MNV file already exist")
}
