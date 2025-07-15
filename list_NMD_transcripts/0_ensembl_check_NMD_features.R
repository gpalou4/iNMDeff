rm(list=ls())
.libPaths( rev( .libPaths() ) )
# Libraries
library("biomaRt")
# ORFs
library("ORFik")
# To read fasta
library("seqinr")
# Reverse complement sequences
library("spgs")
# Read GTF
library("rtracklayer")
# GC content
library("stringr")
# Read excel
library("readxl")

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/list_NMD_transcripts/ensembl_check_NMD_features_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
conversor_tables_path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]

# 1_1) ENSEMBL v88 hg38 GTF
ensembl_v88_gtf <- rtracklayer::import(paste0(conversor_tables_path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
unique(ensembl_v88_gtf$type)
# Remove versions from IDs
ensembl_v88_gtf[,c("gene_id","transcript_id")] <- sapply(ensembl_v88_gtf[,c("gene_id","transcript_id")], function(col) { 
  gsub("(.*)\\..*","\\1", col)
}) 

# 1_2) Riboseq data
uORFs_riboseq <- read_excel(paste0(conversor_tables_path,paths[paths$folder_or_object=="riboseq_Bohlen","path_or_filename"]), sheet = "Data")
uORFs_riboseq <- as.data.frame(uORFs_riboseq)
# 2) Classify transcripts based on NMD-features

ensembl_NMD_features <- ensembl_v88_gtf[,c("transcript_id","transcript_name","gene_id","gene_name")]
colnames(ensembl_NMD_features) <- c("ensembl_transcript_id","ensembl_transcript_name","ensembl_gene_id","gene_symbol")
## Remove NA in ENSEMBL transcript ID and duplicated rows
ensembl_NMD_features <- ensembl_NMD_features[!duplicated(ensembl_NMD_features),]
ensembl_NMD_features <- ensembl_NMD_features[!is.na(ensembl_NMD_features$ensembl_transcript_id),]
ensembl_NMD_features[,c("exons","CDS","UTR5","UTR3","UTR3_GC_content","UTR3_length","stop_codon_to_3UTR_splice_site",
                      "uORFs","uORFs_lengths","uORFs_start_pos","uORFs_translated","uORFs_num_translated","uORFs_trans_reinit","uORFs_num_trans_reinit","uORF_penultimate_codon","NMD_event_type","comment")] <- ""

n <- args[2]

for (i in n:nrow(ensembl_NMD_features)) {
  
  print(i)
  if ((i == 1) | (i == 10000) | (i == 30000) | (i == 60000) | (i == 100000) | (i == 140000) | (i == 180000)) {
     write.table(ensembl_NMD_features, file = paste0(conversor_tables_path,"ensembl_NMD_features_",i,".txt"), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
  }
  ensembl_transcript_id <- as.character(ensembl_NMD_features[i,"ensembl_transcript_id"])
  ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$transcript_id%in%ensembl_transcript_id,]
  if (nrow(ensembl_v88_gtf_filt)==0) {
    print("ENSEMBL transcript ID not found in GTF file")
    ensembl_NMD_features[i,"comment"] <- "ENSEMBL transcript ID not found in GTF file"
    next
  }
  # Strand
  strand <- as.character(unique(ensembl_v88_gtf_filt$strand))
  # Start / Stop codon
  start_codon_start <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$type == "start_codon","start"][1]
  stop_codon <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$type == "stop_codon",c("start","end")]
  # UTRs
  UTRs <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$type == "UTR",c("start","end")]
  # Obtain 3UTRs and 5UTRs
  if (strand == "+") {
    UTR5 <- UTRs[UTRs$start <= start_codon_start,]
    UTR3 <- UTRs[UTRs$end >= stop_codon$end,]
  } else if (strand == "-") {
    UTR5 <- UTRs[UTRs$start >= start_codon_start,]
    UTR3 <- UTRs[UTRs$end <= stop_codon$end,]
  }

  # If start or stop codons don't exist, then initialize UTRs to not enter in its conditions
  if (nrow(stop_codon)==0 || is.na(stop_codon)) {
    ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","no stop codon")
    UTR3 <- data.frame()
  } 
  if (length(start_codon_start)==0 || is.na(start_codon_start)) {
    ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","no start codon")
    UTR5 <- data.frame()
  }

  # 4.1) uORFs
  # ORF found between 5'UTR and start codon. Minimum of 2 ORFs to be NMD-targeted
  # First check if there are UTRs in the gene
  UTR5_num <- nrow(UTR5)
  ensembl_NMD_features[i,"UTR5"] <- as.numeric(UTR5_num)
  if (UTR5_num == 0 ) {
    ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","no 5UTR")
  } else {
    UTR5 <- UTR5[order(UTR5$start),]
    # If >1 UTR
    UTR5_length <- as.numeric(sum(UTR5["end"]-UTR5["start"]))+nrow(UTR5)
    # Prepare BED12 file
    # chrom
    chr <- as.character(unique(ensembl_v88_gtf_filt$seqnames))
    if ( length(grep("GL",chr)) == 1 ) {
      print("ENSEMBL transcript ID not in autossomal or sex chromosome")
      ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","ENSEMBL transcript ID not in autossomal or sex chromosome")
      next
    }
    # start
    start_pos <- UTR5$start[1]-1
    # end
    end_pos <- UTR5$end[nrow(UTR5)]+1
    # score
    score <- "1"
    # thickStart, thickEnd, itemRgb
    thick_start <- thick_end <- item_rgb <- "*"
    # blockCount (UTRs number), blockStarts (UTRs start positions) and blockSizes (UTRs size)
    UTRs_start <- paste0(UTR5$start-UTR5$start[1],",")
    UTRs_number <- length(UTRs_start)
    UTRs_start <- paste0(UTRs_start, collapse = "")
    UTRs_size <- UTR5$end-UTR5$start + 1
    UTRs_size_char <- paste0(UTRs_size, collapse = ",")
    # Create bedfile
    bed_file <- data.frame(chrom=chr, start=start_pos, end=end_pos, name=ensembl_transcript_id, score=score, strand=strand, thickStart=thick_start, 
                            thickEnd=thick_end, itemRgb=item_rgb,blockCount=UTRs_number, blockSizes=UTRs_size_char, blockStarts=UTRs_start)
    setwd("/scratch/1/gpalou/trash")
    write.table(bed_file, file= paste0(ensembl_transcript_id,"_UTR5.bed"), quote=FALSE, sep='\t',row.names=FALSE, col.names = F)
    # Get fasta sequence with UTRs concatenated (bash)
    system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa  -name -bed ",ensembl_transcript_id,"_UTR5.bed > ",ensembl_transcript_id,"_UTR5.fasta"))
    #system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_home/gpalou/data/human_genome/ensembl_87/Homo_sapiens_GRCh37_dna_chromosome_",chr,"_fa", " -name -bed ",ensembl_transcript_id,"_bed > ",ensembl_transcript_id,"_fasta"))

    # Read fasta
    bed12_fasta <- read.fasta(file = paste0(ensembl_transcript_id,"_UTR5.fasta"))
    # Delete fasta
    system(paste0("rm -rf ",ensembl_transcript_id,"_UTR5.bed ",ensembl_transcript_id,"_UTR5.fasta"))
    fasta_5UTR_sequence <- paste0(bed12_fasta[[1]],collapse="")
    # Then check the strand
    if (strand == "-") {
      fasta_5UTR_sequence <- tolower(as.character(reverseComplement(x = DNAString(paste0(fasta_5UTR_sequence, collapse = "")))))
    }
    ORFs <- findORFs(
      fasta_5UTR_sequence,
      startCodon = "atg",
      stopCodon = "taa|tag|tga",
      longestORF = FALSE,
      minimumLength = 0
    )
    ORFs_df <- data.frame(ORFs$`1`)
    num_ORFs <- nrow(ORFs_df)
    ensembl_NMD_features[i,"uORFs"] <- as.numeric(num_ORFs)

    if (num_ORFs >= 1) {
      if (num_ORFs >=2) {
        print(">=2 uORFs")
        ensembl_NMD_features[i,"NMD_event_type"] <- paste0(ensembl_NMD_features[i,"NMD_event_type"]," | ",">=2 uORFs")
      }  
      ORFs_df <- ORFs_df[order(ORFs_df$start),]
      # uORFs lengths
      length_ORFs <- ORFs_df[,"width"]
      length_ORFs <- paste0(length_ORFs, collapse =",")
      ensembl_NMD_features[i,"uORFs_lengths"] <- length_ORFs
      # Store the start positions aswell
      ORFs_start <- ORFs_df[,"start"]
      ORFs_start <- paste0(ORFs_start, collapse =",")
      ensembl_NMD_features[i,"uORFs_start_pos"] <- ORFs_start
      # Check if uORFs are translated using Ribo-seq data
      uORFs_riboseq_transcript <- uORFs_riboseq[uORFs_riboseq[,"Trascript ID"] %in% ensembl_transcript_id,c("Trascript ID","uORF Position","Codon -1","Z-vs-Z Analysis hit ?")]
      if (nrow(uORFs_riboseq_transcript) != 0) {
        uORFs_df_riboseq_merged <- merge(ORFs_df,uORFs_riboseq_transcript, by_x = "start", by_y = "uORF Position", all_x = TRUE)
        # Vector of potentially translated uORFs (0 = no,1 = yes)
        translated_uORFs <- ifelse(is.na(uORFs_df_riboseq_merged[,"Trascript ID"]),"NA",uORFs_df_riboseq_merged[,"Trascript ID"])
        translated_uORFs <- ifelse(translated_uORFs == "NA","NA",1)
        ensembl_NMD_features[i,"uORFs_translated"] <- paste0(translated_uORFs, collapse = ",")
        # Number of potentially translated uORFs
        ensembl_NMD_features[i,"uORFs_num_translated"] <- sum(str_count(translated_uORFs, "1"))
        # Vector of uORFs that promotes translation reinitiation (0 = no,1 = yes)
        trans_reinit_uORFs <- ifelse(is.na(uORFs_df_riboseq_merged[,"Z-vs-Z Analysis hit ?"]),"NA",uORFs_df_riboseq_merged[,"Z-vs-Z Analysis hit ?"])
        ensembl_NMD_features[i,"uORFs_trans_reinit"] <- paste0(trans_reinit_uORFs,collapse=",")
        # Number of uORFs that promotes translation reinitiation
        ensembl_NMD_features[i,"uORFs_num_trans_reinit"] <- sum(str_count(trans_reinit_uORFs, "1"))
        # Vector of penultime codon for uORFs
        penultime_codon_uORFs <- ifelse(is.na(uORFs_df_riboseq_merged[,"Codon -1"]),"NA",uORFs_df_riboseq_merged[,"Codon -1"])
        ensembl_NMD_features[i,"uORF_penultimate_codon"] <- paste0(penultime_codon_uORFs, collapse = ",")
      }
    }
  }
  
  # 4.2) Check the stop codon for each isoform of the same gene (NOT DONE)
  # If the stop codon coordinate doesn't match with the rest of the isoforms and is 50nt upstream of the last exon-exon junction (the longest isoform), it is a PTC
  CDS <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$type == "CDS",c("start","end")]
  CDS_num <- nrow(CDS)
  exons <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$type == "exon",c("start","end")]
  exons_num <- nrow(exons)
  ensembl_NMD_features[i,"exons"] <- as.numeric(exons_num)
  ensembl_NMD_features[i,"CDS"] <- as.numeric(CDS_num)
  if ( (CDS_num == 0) & (exons_num == 0) ) {
    ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","no CDS nor exons")
    
  }

  # 4.3) GC content in 3UTR
  UTR3_num <- nrow(UTR3)
  ensembl_NMD_features[i,"UTR3"] <- as.numeric(UTR3_num)
  if (UTR3_num == 0) {
    print("no 3UTR")
    ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ","no 3UTR")
  } else if ( UTR3_num >= 1 ) {
    UTR3 <- UTR3[order(UTR3$start),]
    # If >1 UTR
    UTR3_length <- as.numeric(sum(UTR3["end"]-UTR3["start"]))+nrow(UTR3)
    # Prepare BED12 file
    # chrom
    chr <- as.character(unique(ensembl_v88_gtf_filt$seqnames))
    # start
    start_pos <- UTR3$start[1]-1
    # end
    end_pos <- UTR3$end[nrow(UTR3)]+1
    # score
    score <- "1"
    # thickStart, thickEnd, itemRgb
    thick_start <- thick_end <- item_rgb <- "*"
    # blockCount (UTRs number), blockStarts (UTRs start positions) and blockSizes (UTRs size)
    UTRs_start <- paste0(UTR3$start-UTR3$start[1],",")
    UTRs_number <- length(UTRs_start)
    UTRs_start <- paste0(UTRs_start, collapse = "")
    UTRs_size <- UTR3$end-UTR3$start + 1
    UTRs_size_char <- paste0(UTRs_size, collapse = ",")
    # Create bedfile
    bed_file <- data.frame(chrom=chr, start=start_pos, end=end_pos, name=ensembl_transcript_id, score=score, strand=strand, thickStart=thick_start, 
                            thickEnd=thick_end, itemRgb=item_rgb,blockCount=UTRs_number, blockSizes=UTRs_size_char, blockStarts=UTRs_start)
    setwd("/scratch/1/gpalou/trash")
    write.table(bed_file, file= paste0(ensembl_transcript_id,"_UTR3.bed"), quote=FALSE, sep='\t',row.names=FALSE, col.names = F)
    # Get fasta sequence with UTRs concatenated (bash)
    system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa  -name -bed ",ensembl_transcript_id,"_UTR3.bed > ",ensembl_transcript_id,"_UTR3.fasta"))
    #system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_home/gpalou/data/human_genome/ensembl_87/Homo_sapiens_GRCh37_dna_chromosome_",chr,"_fa", " -name -bed ",ensembl_transcript_id,"_bed > ",ensembl_transcript_id,"_fasta"))

    # Read fasta
    bed12_fasta <- read.fasta(file = paste0(ensembl_transcript_id,"_UTR3.fasta"))
    # Delete fasta
    system(paste0("rm -rf ",ensembl_transcript_id,"_UTR3.bed ",ensembl_transcript_id,"_UTR3.fasta"))
    fasta_3UTR_sequence <- paste0(bed12_fasta[[1]],collapse="")
    # Then check the strand
    if (strand == "-") {
      fasta_3UTR_sequence <- tolower(as.character(reverseComplement(x = DNAString(paste0(fasta_3UTR_sequence, collapse = "")))))
    }
    # Calculate GC content
    G_count <- str_count(fasta_3UTR_sequence, "g")
    C_count <- str_count(fasta_3UTR_sequence, "c")
    GC_content <- (G_count + C_count) / str_length(fasta_3UTR_sequence) * 100 
    ensembl_NMD_features[i,"UTR3_GC_content"] <- as.numeric(round(GC_content,2))

    # 4.4) Intron-spliced site in 3'UTR
    # Splice site found in the 3'UTR or after the normal stop codon (basically >1 3'UTR coordinates)
    if (UTR3_num >= 2) {
      print("Intronic-spliced 3UTR")
      ensembl_NMD_features[i,"NMD_event_type"] <- paste0(ensembl_NMD_features[i,"NMD_event_type"]," | ","Intronic-spliced 3UTR")
      # Distance
      if (nrow(stop_codon) != 0) {
        if (strand == "+") {
          dist_stop_codon_to_UTR3 <- as.numeric(abs(min(stop_codon["start"])-UTR3[1,"end"]))
        } else if (strand == "-") {
          dist_stop_codon_to_UTR3 <- as.numeric(abs(max(stop_codon["end"])-UTR3[UTR3_num,"start"]))
        }
        if ((nrow(stop_codon) > 1)){
          ensembl_NMD_features[i,"comment"] <- paste0(ensembl_NMD_features[i,"comment"]," | ",">1 stop codon, I chose the most distant")
        }
        ensembl_NMD_features[i,"stop_codon_to_3UTR_splice_site"] <- as.numeric(dist_stop_codon_to_UTR3)
      }
    } 
    # 4.5) Check length of 3'UTR, those needs to be discarded (>500)
    UTR3_length <- sum((UTR3[2]-UTR3[1]))
    ensembl_NMD_features[i,"UTR3_length"] <- UTR3_length
  } 

}
 
write.table(ensembl_NMD_features, file = paste0(conversor_tables_path,paths[paths$folder_or_object=="ensembl_NMD_features","path_or_filename"]), quote=FALSE, sep='\t',row.names=FALSE, col.names = TRUE)
