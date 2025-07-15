# Libraries
library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library("data.table")

# 1) DATA
# 1.1) PRO CAP data
# John Lis, Cornell
# path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_CAP/CoPRO_srrori_Capped_Amplified_merged_1_unidirectional_peaks.bed"
# PRO_CAP <- read.table(file = path, header = FALSE, sep = "\t")

# colnames(PRO_CAP) <- c("chr","start","end","ID","Q_value","strand","peak_reads",
#     NA,NA,NA,NA,NA,NA,NA,NA,"TSS_summit_pos","Height_summit_reads")
# head(PRO_CAP)
# PRO_CAP <- PRO_CAP[,!is.na(colnames(PRO_CAP))]
# head(PRO_CAP)

# # Create an IRanges object with start and end positions
# ranges <- IRanges(start = PRO_CAP$start, end = PRO_CAP$end)

# # Optionally, create a GRanges object to add metadata (like strand)
# PRO_CAP_granges <- GRanges(seqnames = PRO_CAP$chr, ranges = ranges, strand = PRO_CAP$strand)

# # Add metadata columns (e.g., Q-value, TSS_summit_pos)
# mcols(PRO_CAP_granges) <- DataFrame(Q_value = PRO_CAP$Q_value,
#                                 peak_reads = PRO_CAP$peak_reads,
#                                 ID = PRO_CAP$ID,
#                                 TSS_summit_pos = PRO_CAP$TSS_summit_pos,
#                                 Height_summit_reads = PRO_CAP$Height_summit_reads)

# # View the resulting GRanges object
# PRO_CAP_granges

# #Haiyuan Yu, Cornell lab

# path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_CAP/brm_C1a_and_C1b_erm_1_unidirectional_peaks.bed"
# PRO_CAP <- read.table(file = path, header = FALSE, sep = "\t")

# colnames(PRO_CAP) <- c("chr","start","end","ID","Q_value","strand","peak_reads",
#     "TSS_summit_pos","Height_summit_reads")
# head(PRO_CAP)

# # Create an IRanges object with start and end positions
# ranges <- IRanges(start = PRO_CAP$start, end = PRO_CAP$end)

# # Optionally, create a GRanges object to add metadata (like strand)
# PRO_CAP_granges <- GRanges(seqnames = PRO_CAP$chr, ranges = ranges, strand = PRO_CAP$strand)

# # Add metadata columns (e.g., Q-value, TSS_summit_pos)
# mcols(PRO_CAP_granges) <- DataFrame(Q_value = PRO_CAP$Q_value,
#                                 peak_reads = PRO_CAP$peak_reads,
#                                 ID = PRO_CAP$ID,
#                                 TSS_summit_pos = PRO_CAP$TSS_summit_pos,
#                                 Height_summit_reads = PRO_CAP$Height_summit_reads)

# # View the resulting GRanges object
# PRO_CAP_granges

##### PRO-SEQ ####
# A549 --> https://www.encodeproject.org/experiments/ENCSR988BZM/
# 
path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_seq/A549_1.dREG.peak.full.bed"
PRO_seq_A549 <- read.table(file = path, header = FALSE, sep = "\t")

colnames(PRO_seq_A549) <- c("chr","start","end","dREG_score","p_value",
    "Height_summit_pos")
head(PRO_seq_A549)

# Create an IRanges object with start and end positions
ranges <- IRanges(start = PRO_seq_A549$start, end = PRO_seq_A549$end)

# Optionally, create a GRanges object to add metadata (like strand)
PRO_seq_A549_granges <- GRanges(seqnames = PRO_seq_A549$chr, ranges = ranges, strand = PRO_seq_A549$strand)

# Add metadata columns (e.g., Q-value, TSS_summit_pos)
mcols(PRO_seq_A549_granges) <- DataFrame(p_value = PRO_seq_A549$p_value,
                                dREG_score = PRO_seq_A549$dREG_score,
                                Height_summit_pos = PRO_seq_A549$Height_summit_pos)

# View the resulting GRanges object
PRO_seq_A549_granges

# U2OS --> https://www.encodeproject.org/experiments/ENCSR528HBJ/

path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_seq/U2OS_1.dREG.peak.full.bed"
PRO_seq_U2OS <- read.table(file = path, header = FALSE, sep = "\t")

colnames(PRO_seq_U2OS) <- c("chr","start","end","dREG_score","p_value",
    "Height_summit_pos")
head(PRO_seq_U2OS)

# Create an IRanges object with start and end positions
ranges <- IRanges(start = PRO_seq_U2OS$start, end = PRO_seq_U2OS$end)

# Optionally, create a GRanges object to add metadata (like strand)
PRO_seq_U2OS_granges <- GRanges(seqnames = PRO_seq_U2OS$chr, ranges = ranges, strand = PRO_seq_U2OS$strand)

# Add metadata columns (e.g., Q-value, TSS_summit_pos)
mcols(PRO_seq_U2OS_granges) <- DataFrame(p_value = PRO_seq_U2OS$p_value,
                                dREG_score = PRO_seq_U2OS$dREG_score,
                                Height_summit_pos = PRO_seq_U2OS$Height_summit_pos)

# View the resulting GRanges object
PRO_seq_U2OS_granges

summary(PRO_seq_U2OS_granges@ranges@width)

# Neuroblastoma SH-SY5Y

path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_seq/GSM6602003_Control_FWD.bigwig"
PRO_seq_SH_SY5Y_granges_1 <- import(path, format = "bigwig")
head(PRO_seq_SH_SY5Y_granges_1)
PRO_seq_SH_SY5Y_granges_1$score <- round(PRO_seq_SH_SY5Y_granges_1$score,2)
summary(PRO_seq_SH_SY5Y_granges_1$score)
PRO_seq_SH_SY5Y_granges_1 <- PRO_seq_SH_SY5Y_granges_1[PRO_seq_SH_SY5Y_granges_1$score > 15.99,]
length(PRO_seq_SH_SY5Y_granges_1$score)

path <- "/g/strcombio/fsupek_cancer1/gpalou/PRO/PRO_seq/GSM6602004_Control_REV.bigwig"
PRO_seq_SH_SY5Y_granges_2 <- import(path, format = "bigwig")
head(PRO_seq_SH_SY5Y_granges_2)
PRO_seq_SH_SY5Y_granges_2$score <- round(PRO_seq_SH_SY5Y_granges_2$score,2)
summary(PRO_seq_SH_SY5Y_granges_2$score)
PRO_seq_SH_SY5Y_granges_2 <- PRO_seq_SH_SY5Y_granges_2[PRO_seq_SH_SY5Y_granges_2$score > 17.18,]
length(PRO_seq_SH_SY5Y_granges_2$score)


summary(PRO_seq_SH_SY5Y_granges_1@ranges@width)
summary(PRO_seq_SH_SY5Y_granges_2@ranges@width)

# 1.2) NMD targets table

path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_old.txt"
NMD_targets <- read.table(file = path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Calculate the difference in transcription start sites between two transcripts for each gene
NMD_targets <- NMD_targets %>%
  group_by(gene_symbol) %>%
  filter(n() == 2) %>%  # Ensure there are exactly two transcripts for each gene
  mutate(diff_TSS = abs(transcription_start_site[2] - transcription_start_site[1])) %>%
  ungroup() %>% arrange(ensembl_gene_id) %>% as.data.frame()

# View the result
head(NMD_targets[, c("gene_symbol", "transcription_start_site", "diff_TSS")])

# 1.3) GFF used for gene annotation in RNA-seq

# Define the path to your GTF file
path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf"

# Read the GTF file into R
gtf_data <- import.gff(path)

# View the first few rows
gtf_data

# Filter by NMD targets
gtf_data$gene_id <- sub("\\..*", "", gtf_data$gene_id)
gtf_NMD_targets <- gtf_data[gtf_data$gene_id %in% unique(NMD_targets$ensembl_gene_id),]
gtf_NMD_targets <- gtf_NMD_targets[gtf_NMD_targets$type == "transcript",]

# 1.4) RNAseq matrix
# TCGA
TCGA_RNAseq_TPM_pancancer <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_transcript.txt", 
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print("Dimensions -->")
print(dim(TCGA_RNAseq_TPM_pancancer))
colnames(TCGA_RNAseq_TPM_pancancer) <- gsub("\\.","-",colnames(TCGA_RNAseq_TPM_pancancer))
# GTEx
GTEx_RNAseq_TPM_pantissue <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct",
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2)                         
print("Dimensions -->")
print(dim(GTEx_RNAseq_TPM_pantissue))
# Fix Genes ID (Remove genes duplicated in Y chr)
GTEx_RNAseq_TPM_pantissue <- GTEx_RNAseq_TPM_pantissue[-grep("PAR",GTEx_RNAseq_TPM_pantissue$transcript_id),]
keep <- c(1,grep("GTEX",colnames(GTEx_RNAseq_TPM_pantissue)))
GTEx_RNAseq_TPM_pantissue <- GTEx_RNAseq_TPM_pantissue[,keep]
rownames(GTEx_RNAseq_TPM_pantissue) <- gsub("(\\..*)","",GTEx_RNAseq_TPM_pantissue$transcript_id)
GTEx_RNAseq_TPM_pantissue$transcript_id <- NULL

# 1.5) FANTOM5 CAGE data

path <- "/g/strcombio/fsupek_cancer1/gpalou/FANTOM5_CAGE/promoter/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt"

# 2) Analysis

# 2.1) Find overlaps between PRO peaks and gene annotations

# These overlaps mean peaks that overlap TSS or are contained within the range of the transcript.
# Afterward we will classify those peaks into UTRs-only (do not overlap TSS), or TSS-proximal
# overlaps <- findOverlaps(PRO_CAP_granges, gtf_NMD_targets)

# Define proximity threshold in base pairs (e.g., 500 bp upstream of TSS)
proximity_threshold <- 0
# Extend the TSS ranges to include the proximity threshold
extended_gtf_ranges <- gtf_NMD_targets
start(extended_gtf_ranges) <- start(extended_gtf_ranges) - proximity_threshold
# Ensure no negative values (adjusting for chromosome starts)
# start(extended_gtf_ranges) <- pmax(start(extended_gtf_ranges), 1)
# Now find overlaps between the PRO_CAP peaks and the extended TSS regions

PRO_seq_peaks_detection <- function(cell_type) {

    if (cell_type == "A549") {
        Granges_object <- PRO_seq_A549_granges # change by one of our PRO-CAP or PRO-seq objects
        colnames(mcols(Granges_object))[colnames(mcols(Granges_object)) == "dREG_score"] <- "score"
    } else if (cell_type == "U2OS") {
        Granges_object <- PRO_seq_U2OS_granges
        colnames(mcols(Granges_object))[colnames(mcols(Granges_object)) == "dREG_score"] <- "score"
    } else if (cell_type == "SH_SY5Y_1") {
        Granges_object <- PRO_seq_SH_SY5Y_granges_1
        # Calculate the middle position and add it as a new column
        mcols(Granges_object)$Height_summit_pos <- start(Granges_object@ranges) + 
                                                floor((width(Granges_object@ranges) - 1) / 2)
    } else if (cell_type == "SH_SY5Y_2") {
        Granges_object <- PRO_seq_SH_SY5Y_granges_2
        # Calculate the middle position and add it as a new column
        mcols(Granges_object)$Height_summit_pos <- start(Granges_object@ranges) + 
                                                floor((width(Granges_object@ranges) - 1) / 2)
    }

    overlaps <- findOverlaps(Granges_object, extended_gtf_ranges)
    overlaps_df <- data.frame(overlaps)

    # Convert data to data.table for faster manipulation
    Granges_dt <- as.data.table(Granges_object)
    extended_gtf_ranges_dt <- as.data.table(extended_gtf_ranges)
    extended_gtf_ranges_dt[,c("score","width")] <- NULL
    Granges_dt[,c("strand")] <- NULL

    # 1) Overlap PRO-seq peaks with NMD targets
    PRO_with_ann <- overlaps_df %>%
        mutate(
            peak_row = queryHits,
            gene_row = subjectHits
        ) %>%
        rowwise() %>%
        do({
            peak_data <- Granges_dt[.data$peak_row]
            gene_data <- extended_gtf_ranges_dt[.data$gene_row]
            gene_data <- rename(gene_data, chr = seqnames, start_gene = start, end_gene = end)
            cbind(peak_data, gene_data)
        }) %>%
        ungroup() %>%
        bind_rows()

    PRO_with_ann_filt <- PRO_with_ann[,c("seqnames","start","end","score", "Height_summit_pos",
                    "start_gene","end_gene","type","gene_id","gene_type","gene_name",
                    "transcript_id","strand")]
    PRO_with_ann_filt$transcript_id <- sub("\\..*", "", PRO_with_ann_filt$transcript_id)
    # Merge
    PRO_with_ann_NMD_targets <- merge(NMD_targets,PRO_with_ann_filt, 
            by.x = "ensembl_transcript_id", 
            by.y = "transcript_id", all.x = TRUE)
    table(is.na(PRO_with_ann_NMD_targets$score))
    PRO_with_ann_NMD_targets <- PRO_with_ann_NMD_targets %>%
        arrange(ensembl_gene_id,ensembl_transcript_id)

    # 2) Classify peaks if they are 5'UTR or TSS overlapping

    classify_peak_with_distance <- function(transcription_start_site, transcription_end_site, peak_start, peak_end, UTR_data, strand) {
    
        # Check if the peak overlaps TSS (start before TSS and end after)
        if (strand == "+") {
            overlaps_TSS_only <- peak_start <= transcription_start_site & peak_end >= transcription_start_site
        } else if (strand == "-") {
            overlaps_TSS_only <- peak_end >= transcription_start_site & peak_start <= transcription_start_site
        }
    
        overlaps_5UTR_only <- FALSE
        overlaps_TSS_and_5UTR <- FALSE

        for (i in 1:nrow(UTR_data)) {
            # overlaps_5UTR_only with 1000bp distance from peak end to UTR end
            if (strand == "+" && peak_start >= UTR_data$start[i] && peak_end <= UTR_data$end[i]) {
                overlaps_5UTR_only <- TRUE
            } else if (strand == "-" && peak_end <= UTR_data$end[i] && ( abs(peak_start - UTR_data$start[i]) < 1000) ) {
                overlaps_5UTR_only <- TRUE
            }
            # overlaps_TSS_and_5UTR with 1000bp distance from peak end to UTR end
            if (strand == "+" && peak_start <= transcription_start_site && ( abs(peak_end - UTR_data$end[i]) < 1000) ) {
                overlaps_TSS_and_5UTR <- TRUE
            } else if (strand == "-" && peak_end >= transcription_start_site && ( abs(peak_start - UTR_data$start[i]) < 1000) )  {
                overlaps_TSS_and_5UTR <- TRUE
            }        
        }

        # Classify peak
        if (overlaps_TSS_and_5UTR) {
            return("Both TSS and 5'UTR")
        } else if (overlaps_TSS_only) {
            return("Only TSS")
        } else if (overlaps_5UTR_only) {
            return("Only 5'UTR")
        } else {
            return("No Overlap")
        }
    }

    # UTR information
    gtf_NMD_targets_UTR <- gtf_data[gtf_data$gene_id %in% unique(NMD_targets$ensembl_gene_id),]
    gtf_NMD_targets_UTR <- gtf_NMD_targets_UTR[gtf_NMD_targets_UTR$type == "UTR",]
    gtf_NMD_targets_UTR$transcript_id <- sub("\\..*", "", gtf_NMD_targets_UTR$transcript_id)

    PRO_with_ann_NMD_targets$PRO_seq_peak_class <- NA

    for (i in 1:nrow(PRO_with_ann_NMD_targets)) {
        print(round(i/nrow(PRO_with_ann_NMD_targets),2)*100)
        PRO_peaks_tmp <- PRO_with_ann_NMD_targets[i,]
        if (is.na(PRO_peaks_tmp$start)) { next }
        transcript <-  PRO_peaks_tmp$ensembl_transcript_id
        transcript_UTR_tmp <- data.frame(gtf_NMD_targets_UTR[which(gtf_NMD_targets_UTR$transcript_id %in% transcript),])[,c("start","end","width","strand")]
        # Is the peak CONTAINED (not overlapping) within any of the UTRs from the gene? --> UTR peak
        # Else --> TSS overlapping
        # Apply classification function
        strand <- unique(as.character(transcript_UTR_tmp$strand))
        PRO_seq_peak_class <- mapply(classify_peak_with_distance,
                                                transcription_start_site = PRO_peaks_tmp$transcription_start_site,
                                                transcription_end_site = PRO_peaks_tmp$transcription_end_site,
                                                peak_start = PRO_peaks_tmp$start,
                                                peak_end = PRO_peaks_tmp$end,
                                                UTR_data = list(transcript_UTR_tmp),  # Provide UTR data as list for each row
                                                strand = strand)
        PRO_with_ann_NMD_targets[i,"PRO_seq_peak_class"] <- PRO_seq_peak_class
        PRO_with_ann_NMD_targets[i,"strand"] <- strand
    }
    head(PRO_with_ann_NMD_targets)

    # 3) Use the Summit of the peak to decide to which transcript does the peak corresponds

    # Calculate the absolute difference between each 'Height_summit_pos' and every 'transcription_start_site'
    PRO_with_ann_NMD_targets <- PRO_with_ann_NMD_targets %>%
        group_by(ensembl_gene_id) %>%
          mutate(PRO_seq_closest_TSS_diff = abs(Height_summit_pos - transcription_start_site)) %>%
        # mutate(PRO_seq_closest_TSS_diff = abs(start - transcription_start_site)) %>%
        ungroup()
    # 
    # For each gene, identify the "main" peak (i.e., the one with the smallest difference to the TSS)
    PRO_with_ann_NMD_targets <- PRO_with_ann_NMD_targets %>%
        #   group_by(ensembl_gene_id,Height_summit_pos) %>%
        group_by(ensembl_gene_id,start) %>%
        mutate(peak_status = ifelse(PRO_seq_closest_TSS_diff == min(PRO_seq_closest_TSS_diff), "main", "not_main")) %>%
        ungroup() %>% data.frame()


    # 4) Add gene expression from brain-only samples

    # Calculate the median of expression across Brain vs rest of tissues for each NMD-target/control transcript
    calculate_medians <- function(transcript, brain_samples, RNAseq_data) {
        # Extract data for the given transcript and brain vs non-brain samples
        transcript_data <- RNAseq_data[which(rownames(RNAseq_data) %in% transcript), ]
        # Get expression data for brain and non-brain samples
        brain_data <- transcript_data[,brain_samples]
        non_brain_data <- transcript_data[,!colnames(transcript_data) %in% brain_samples]
        # Calculate median for brain samples
        transcript_median_brain <- median(as.numeric(brain_data), na.rm = TRUE)
        # Calculate median for non-brain samples
        transcript_median_rest <- median(as.numeric(non_brain_data), na.rm = TRUE)
        # Perform a t-test to compare the means between brain and non-brain groups
        wilcox_result <- wilcox.test(as.numeric(brain_data), as.numeric(non_brain_data), alternative = "greater")
        # shapiro_test <- shapiro.test(as.numeric(brain_data))
        # shapiro_test$p.value
        return(c(median_brain = transcript_median_brain, median_rest = transcript_median_rest, wilcox_pval = wilcox_result$p.value))
    }

    # 4.1) TCGA
    sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
    sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
    brain_samples_TCGA <- sample_NMD_efficiencies_TCGA %>%
                    filter(cancer_type %in% c("TCGA-LGG","TCGA-GBM")) %>%
                    pull(sample) %>% as.character()

    table(colnames(TCGA_RNAseq_TPM_pancancer) %in% brain_samples_TCGA)

    # Apply the function across all rows of PRO_with_ann_NMD_targets using `apply()` or `sapply()`
    TCGA_medians <- t(sapply(unique(PRO_with_ann_NMD_targets$ensembl_transcript_id), function(transcript) {
        calculate_medians(transcript, brain_samples_TCGA, TCGA_RNAseq_TPM_pancancer)
    }))

    PRO_with_ann_NMD_targets <- merge(PRO_with_ann_NMD_targets,TCGA_medians, by.x = "ensembl_transcript_id", by.y = "row.names")
    PRO_with_ann_NMD_targets$TCGA_brain_enriched <- PRO_with_ann_NMD_targets$wilcox_pval < 0.05
    colnames(PRO_with_ann_NMD_targets)[colnames(PRO_with_ann_NMD_targets) %in% colnames(TCGA_medians)] <- c("TCGA_TPM_median_exp_in_brain","TCGA_TPM_median_exp_in_rest","TCGA_wilcox_pval")
    head(PRO_with_ann_NMD_targets)

    table(PRO_with_ann_NMD_targets$TCGA_brain_enriched,
            PRO_with_ann_NMD_targets$final_consensus)

    # 4.2) GTEx
    sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
    sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    brain_samples_GTEx <- sample_NMD_efficiencies_GTEx %>%
                    filter(acronyms %in% c("BRNACC","BRNAMY","BRNCDT","BRNCHA","BRNCHB","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC",
                            "BRNPTM","BRNSNG","BRNSPC")) %>%
                    pull(sample_full_barcode) %>% as.character()
    colnames(GTEx_RNAseq_TPM_pantissue) <- gsub("\\.", "-", colnames(GTEx_RNAseq_TPM_pantissue))

    table(colnames(GTEx_RNAseq_TPM_pantissue) %in% brain_samples_GTEx)

    # Apply the function across all rows of PRO_with_ann_NMD_targets using `apply()` or `sapply()`
    GTEx_medians <- t(sapply(unique(PRO_with_ann_NMD_targets$ensembl_transcript_id), function(transcript) {
        calculate_medians(transcript, brain_samples_GTEx, GTEx_RNAseq_TPM_pantissue)
    }))

    PRO_with_ann_NMD_targets <- merge(PRO_with_ann_NMD_targets,GTEx_medians, by.x = "ensembl_transcript_id", by.y = "row.names")
    PRO_with_ann_NMD_targets$GTEx_brain_enriched <- PRO_with_ann_NMD_targets$wilcox_pval < 0.05
    colnames(PRO_with_ann_NMD_targets)[colnames(PRO_with_ann_NMD_targets) %in% colnames(GTEx_medians)] <- c("GTEx_TPM_median_exp_in_brain","GTEx_TPM_median_exp_in_rest","GTEx_wilcox_pval")
    head(PRO_with_ann_NMD_targets)

    table(PRO_with_ann_NMD_targets$GTEx_brain_enriched,
            PRO_with_ann_NMD_targets$final_consensus)

    # 5) Save
    path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PRO_seq/PRO_seq_",cell_type,"_NMD_targets.txt")
    # PRO_with_ann_NMD_targets <- read.table(file = path, header = TRUE, sep = "\t")
    write.table(PRO_with_ann_NMD_targets, file = path, 
                sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)

    return(PRO_with_ann_NMD_targets)

}

PRO_seq_SH_SY5Y_1 <- PRO_seq_peaks_detection(cell_type = "SH_SY5Y_1")
PRO_seq_SH_SY5Y_1$cell_type <- "SH_SY5Y_1"
PRO_seq_SH_SY5Y_2 <- PRO_seq_peaks_detection(cell_type = "SH_SY5Y_2")
PRO_seq_SH_SY5Y_2$cell_type <- "SH_SY5Y_2"
PRO_seq_U2OS <- PRO_seq_peaks_detection(cell_type = "U2OS")
PRO_seq_U2OS$cell_type <- "U2OS"
PRO_seq_A549 <- PRO_seq_peaks_detection(cell_type = "A549")
PRO_seq_A549$cell_type <- "A549"

# 3) Filters
PRO_with_ann_NMD_targets_all <- rbind(PRO_seq_SH_SY5Y_1, PRO_seq_SH_SY5Y_2, PRO_seq_U2OS, PRO_seq_A549)
path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PRO_seq/PRO_seq_all_cell_lines_NMD_targets.txt")
PRO_with_ann_NMD_targets_all <- read.table(file = path, header = TRUE, sep = "\t")
# write.table(PRO_with_ann_NMD_targets_all, file = path, 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)

PRO_with_ann_NMD_targets_final <- c()

# combined_scores <- PRO_with_ann_NMD_targets_all %>%
#   filter(PRO_seq_peak_class != "No Overlap", PRO_seq_closest_TSS_diff <= 1000) %>%
#   pull(score)  # Extract all scores
# target_vector <- sort(combined_scores, decreasing = FALSE)  # Use the sorted values as the target

# library(preprocessCore)  # For quantile normalization

for (cell in c("SH_SY5Y_1","SH_SY5Y_2","U2OS","A549")) {

    PRO_with_ann_NMD_targets_tmp <- PRO_with_ann_NMD_targets_all %>%
                        filter(cell_type == cell)
    # 3.1) Filter Peaks that are inside the gene but are not in the 5'UTR
    PRO_with_ann_NMD_targets_filt <- PRO_with_ann_NMD_targets_tmp %>%
        filter(PRO_seq_peak_class != "No Overlap") %>%
        # filter(peak_status == "main" & PRO_seq_closest_TSS_diff <= 100) %>% # Distance between Peak and TSS < 50bp
        filter(PRO_seq_closest_TSS_diff <= 1000) %>% # Distance between Peak and TSS < 50bp
        # filter(score > quantile(score, 0.25, na.rm = TRUE)) %>% # Remove rows with score <= 1st quantile
        arrange(ensembl_gene_id,final_consensus)
    dim(PRO_with_ann_NMD_targets_filt)
    # head(PRO_with_ann_NMD_targets_filt[,c("ensembl_transcript_id","PRO_seq_closest_TSS_diff","score","Height_summit_pos")])
    # Filter the row with the lowest PRO_seq_closest_TSS_diff and highest score
    PRO_with_ann_NMD_targets_filt <- PRO_with_ann_NMD_targets_filt %>%
        group_by(ensembl_transcript_id) %>%
        filter(PRO_seq_closest_TSS_diff == min(PRO_seq_closest_TSS_diff)) %>%
        #   slice_max(score, with_ties = FALSE) %>%
        ungroup()
    dim(PRO_with_ann_NMD_targets_filt)
    # Check for each gene if the 'Height_summit_pos' is the same for both 'control' and 'NMD_target'
    PRO_with_ann_NMD_targets_filt <- PRO_with_ann_NMD_targets_filt %>%
        group_by(ensembl_gene_id) %>%
        filter(!any(duplicated(Height_summit_pos) & final_consensus %in% c("control", "NMD_target"))) %>%
        ungroup() %>% data.frame()
    dim(PRO_with_ann_NMD_targets_filt)
   
    # Scale
    PRO_with_ann_NMD_targets_filt$score <- scale(log10(PRO_with_ann_NMD_targets_filt$score))

    # Quantile normalize the 'score' column
    # PRO_with_ann_NMD_targets_filt$score <- scale(normalize.quantiles.use.target(as.matrix(PRO_with_ann_NMD_targets_filt$score),
    #                                                                       target = target_vector)

    dim(PRO_with_ann_NMD_targets_filt)
    unique(PRO_with_ann_NMD_targets_filt$ensembl_gene_id)
    table(PRO_with_ann_NMD_targets_filt$final_consensus)

    # Re-merge
    PRO_with_ann_NMD_targets_final <- rbind(PRO_with_ann_NMD_targets_final,PRO_with_ann_NMD_targets_filt)
 
}

table(PRO_with_ann_NMD_targets_final$cell_type)

# Plot difference of Pro-Seq peaks between Brain-enriched and not Brain-enriched
PRO_with_ann_NMD_targets_final$TCGA_brain_enriched <- factor(PRO_with_ann_NMD_targets_final$TCGA_brain_enriched, 
                                                           levels = c(TRUE, FALSE), 
                                                           labels = c("Brain-enriched", "Non-brain-enriched"))
PRO_with_ann_NMD_targets_final$GTEx_brain_enriched <- factor(PRO_with_ann_NMD_targets_final$GTEx_brain_enriched, 
                                                           levels = c(TRUE, FALSE), 
                                                           labels = c("Brain-enriched", "Non-brain-enriched"))

# Plot dREG score density for TCGA (using TCGA_brain_enriched for fill)
plot <- ggplot(PRO_with_ann_NMD_targets_final, 
               aes(y = score, x = TCGA_brain_enriched, fill = factor(TCGA_brain_enriched),
                   color = final_consensus)) +
  # Boxplot with transparency
  geom_boxplot(alpha = 0.7) +
  # Jittered points
  geom_jitter(aes(color = final_consensus), 
              position = position_jitterdodge(), 
              size = 2, alpha = 0.5) +
  # Modern color palette for fill and color
  scale_fill_brewer(palette = "Set2") +  # Updated softer palette
  scale_color_brewer(palette = "Dark2") +  # Updated softer palette for consensus
  # Facet customization
  facet_wrap(. ~ cell_type, 
             scales = "free", 
             strip.position = "bottom") +  # Move panel titles below
  theme(strip.background = element_blank(),  # Remove panel background
        strip.text = element_text(size = 10, face = "bold")) +  # Keep titles
  # Add horizontal line and p-values
  geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3, label.y = 2 ) +
  # Labels and axis adjustments
  labs(title = "TCGA", fill = "", color = "Transcript",
       x = "", y = expression(log[10]("PRO-Seq score"))) +
  theme_classic() +
  guides(fill = "none") +
  theme(
    # legend.title = element_blank(), 
    plot.title = element_text(hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5), # Rotate X-axis labels
    strip.background = element_blank(),  # Remove strip background
    strip.placement = "outside"          # Keep facet titles as outer labels
  ) +
  # Reorder facets
  facet_wrap(~factor(cell_type, levels = c("SH_SY5Y_1", "SH_SY5Y_2", "U2OS", "A549")))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PRO_seq/TCGA_PRO_seq_diff_Brain_enriched_vs_not.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

# Plot dREG score density for TCGA (using TCGA_brain_enriched for fill)
plot <- ggplot(PRO_with_ann_NMD_targets_final, 
               aes(y = score, x = GTEx_brain_enriched, fill = factor(GTEx_brain_enriched),
                   color = final_consensus)) +
  # Boxplot with transparency
  geom_boxplot(alpha = 0.7) +
  # Jittered points
  geom_jitter(aes(color = final_consensus), 
              position = position_jitterdodge(), 
              size = 2, alpha = 0.5) +
  # Modern color palette for fill and color
  scale_fill_brewer(palette = "Set2") +  # Updated softer palette
  scale_color_brewer(palette = "Dark2") +  # Updated softer palette for consensus
  # Facet customization
  facet_wrap(. ~ cell_type, 
             scales = "free", 
             strip.position = "bottom") +  # Move panel titles below
  theme(strip.background = element_blank(),  # Remove panel background
        strip.text = element_text(size = 10, face = "bold")) +  # Keep titles
  # Add horizontal line and p-values
  geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3, label.y = 2 ) +
  # Labels and axis adjustments
  labs(title = "GTex", fill = "", color = "Transcript",
       x = "", y = expression(log[10]("PRO-Seq score"))) +
  theme_classic() +
  guides(fill = "none") +
  theme(
    # legend.title = element_blank(), 
    plot.title = element_text(hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5), # Rotate X-axis labels
    strip.background = element_blank(),  # Remove strip background
    strip.placement = "outside"          # Keep facet titles as outer labels
  ) +
  # Reorder facets
  facet_wrap(~factor(cell_type, levels = c("SH_SY5Y_1", "SH_SY5Y_2", "U2OS", "A549")))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PRO_seq/GTEx_PRO_seq_diff_Brain_enriched_vs_not.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

# save R object
write.table(PRO_with_ann_NMD_targets_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3_new/SuppFig3B_C.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(PRO_with_ann_NMD_targets_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3_new/SuppFig3B_C.RData")
