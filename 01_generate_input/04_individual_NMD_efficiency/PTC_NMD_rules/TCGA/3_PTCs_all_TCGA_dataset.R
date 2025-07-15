library(dplyr)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
type <- args[1]

TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/[X1]/"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1

# 1) Merge all PTCs from all TCGA samples into one dataframe
#TCGA_cancers <- TCGA_cancers[!TCGA_cancers %in% "TCGA-SKCM"]

#assign("last.warning", NULL, envir = baseenv())

PTC_transcripts_all_TCGA <- c()
for (TCGA_cancer in TCGA_cancers) {
    print(TCGA_cancer)
    TCGA_cancer_path <- gsub("\\[X1\\]",TCGA_cancer,TCGA_cancers_path)
    TCGA_samples <- list.files(TCGA_cancer_path)
    PTC_transcripts_cancer <- c()
    for (TCGA_sample in TCGA_samples) {
        # Check if PTC file exists
        PTC_transcripts_file_path <- paste0(TCGA_cancer_path,TCGA_sample,"/PTC_transcripts/",type,"_PTC_transcripts_metadata.txt")
        #print(PTC_transcripts_file_path)
        print(TCGA_sample)
        #print(warnings())
        if (file.exists(PTC_transcripts_file_path)) {
            # Open PTCs file
            PTC_transcripts_file <- read.table(file = PTC_transcripts_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            # For the time being I will remove columns related to PCA subtype, as the number is different among different cancer types
            # We don't need the nucleotide sequence either
            #PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue|fasta_sequence_wt|fasta_sequence_mut",colnames(PTC_transcripts_file))]
            PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue",colnames(PTC_transcripts_file))]
            # Add TCGA barcode as a new variable
            PTC_transcripts_file_filt$TCGA_barcode <- TCGA_sample
            # Add TCGA cancer as new variable
            PTC_transcripts_file_filt$TCGA_cancer <- TCGA_cancer
            if ( is.null(PTC_transcripts_file_filt$germline_SNV) ) {
                PTC_transcripts_file_filt$germline_SNV <- "no"
            }
            #print(dim(PTC_transcripts_file))
            # Join data frame 
            PTC_transcripts_cancer <- rbind(PTC_transcripts_cancer,PTC_transcripts_file_filt)
        }
    }
    print(warnings())
    PTC_transcripts_all_TCGA <- rbind(PTC_transcripts_all_TCGA,PTC_transcripts_cancer)
}

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_all_TCGA_seq.txt")
write.table(PTC_transcripts_all_TCGA, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
#PTC_transcripts_all_TCGA <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)

#table(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM == "Inf")[2]/nrow(PTC_transcripts_all_TCGA_filt)
#table(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM == "-Inf")[2]/nrow(PTC_transcripts_all_TCGA_filt)

# 2) Perform filterings to obtain a confident set of PTCs
median(PTC_transcripts_all_TCGA$NMD_efficiency_TPM,na.rm = TRUE)
# 2.1) Remove NAs in NMDeff
dim(PTC_transcripts_all_TCGA)
PTC_transcripts_all_TCGA$NMD_efficiency_TPM <- ifelse(PTC_transcripts_all_TCGA$NMD_efficiency_TPM %in% c("Inf","-Inf"),NA,PTC_transcripts_all_TCGA$NMD_efficiency_TPM)
table(is.na(PTC_transcripts_all_TCGA$NMD_efficiency_TPM))
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA[!is.na(PTC_transcripts_all_TCGA$NMD_efficiency_TPM),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.2) Remove FS indel mutations that do not have a predicted downstream PTC
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[-which(is.na(PTC_transcripts_all_TCGA_filt$PTC_CDS_pos)),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.3) Coefficient of variation (mean/SD) > 0.5 --> This reduces a lot the NMDeff values (from 1.6 raw to 0.3)
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$coeff_var <= 1),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.4) At least 10 samples per WT-transcript per RNAseq NMF cluster
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$TPM_WT_samples >= 5),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.5) Ovelapping additional germline nonsense, FS or splice site variants
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$germline_SNV == "no"),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.6) Overlapping somatic CNV/mut 
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$somatic_CNV_SNV == "no"),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.7) Genes with 1 exon and without splice site at 3'UTR
filter <- which(PTC_transcripts_all_TCGA_filt$transcript_CDS_exon_num == 1 & PTC_transcripts_all_TCGA_filt$splice_site_3UTR == "no")
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[-filter,]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.8) LOEUF score (first decile are the most neg selected ones)
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[-which(PTC_transcripts_all_TCGA_filt$LOEUF_decile == 0),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.9) Positive Selected Genes (PSG)
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$PSG == "no"),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.10) Keep only Heterozygous variants
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$genotype %in% c("0/1","0|1","1|0")),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.11) Low expressed transcripts (<5 TPM)
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$median_TPM_exp_transcript >= 5),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.12) Variant Allele Frequency (VAF) of 0.2 
if (type == "somatic") {
    PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[which(PTC_transcripts_all_TCGA_filt$VAF >= 0.2),]
}
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
summary(PTC_transcripts_all_TCGA_filt$VAF)

# 2.13) Add turn_over datasets
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome.csv"
half_life_1 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_2.csv"
half_life_2 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_bcells.csv"
half_life_3 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_PCA.csv"
half_life_4 <- read.table(file = file, header = TRUE, sep = ",")
# merge
vec1 <- unlist(lapply(strsplit(PTC_transcripts_all_TCGA_filt$gene_id,";"), function(gene){gene[1]}))
vec2 <- unlist(lapply(strsplit(PTC_transcripts_all_TCGA_filt$gene_id,";"), function(gene){gene[2]}))
PTC_transcripts_all_TCGA_filt$gene_id_1 <- vec1
PTC_transcripts_all_TCGA_filt$gene_id_2 <- vec2
PTC_transcripts_all_TCGA_filt <- merge(PTC_transcripts_all_TCGA_filt,half_life_4[,c("Ensembl.Gene.Id","half.life..PC1.")], by.x ="gene_id_1", by.y ="Ensembl.Gene.Id", all.x = TRUE)
PTC_transcripts_all_TCGA_filt <- merge(PTC_transcripts_all_TCGA_filt,half_life_4[,c("Ensembl.Gene.Id","half.life..PC1.")], by.x ="gene_id_2", by.y ="Ensembl.Gene.Id", all.x = TRUE)
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt %>% 
    mutate(half_life = coalesce(half.life..PC1..x,half.life..PC1..y))
PTC_transcripts_all_TCGA_filt[,c("half.life..PC1..x","half.life..PC1..y")] <- NULL

# 2.14) Remove repeated PTCs
# 1) PTCs can overlap different transcripts within the same sample
# For each same PTC start position and SAME sample let's take one transcript
ensembl_transcripts_MANE_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_MANE_v107.txt"
ensembl_transcripts_MANE <- read.table(file = ensembl_transcripts_MANE_path, header = TRUE)
# First the MANE Select if possible, if not, then randomly
ensembl_transcripts_MANE$MANE_main_isoform <- "yes"
PTC_transcripts_all_TCGA_filt <- merge(PTC_transcripts_all_TCGA_filt,ensembl_transcripts_MANE[,c("ensembl_transcript_id","MANE_main_isoform")],
        by.x="transcript_id", by.y ="ensembl_transcript_id", all.x = TRUE)
PTC_transcripts_all_TCGA_filt$MANE_main_isoform <- ifelse(is.na(PTC_transcripts_all_TCGA_filt$MANE_main_isoform),"no","yes")
PTC_transcripts_all_TCGA_clean <- aggregate(NMD_efficiency_TPM ~ TCGA_cancer + TCGA_barcode + PTC_CDS_pos + start_pos + MANE_main_isoform, data = PTC_transcripts_all_TCGA_filt, FUN=function(x) x[sample(seq_along(x), 1)])
# Remove those transcripts with both MANE + random (keep onle the MANE ones) ~ 664 transcripts
PTC_transcripts_all_TCGA_clean <- PTC_transcripts_all_TCGA_clean[-which(duplicated(PTC_transcripts_all_TCGA_clean[,1:4], fromLast = TRUE)),]
# Merge again
common_col <- c("TCGA_cancer","TCGA_barcode","NMD_efficiency_TPM","PTC_CDS_pos","start_pos","MANE_main_isoform")
PTC_transcripts_all_TCGA_clean <- merge(PTC_transcripts_all_TCGA_clean,PTC_transcripts_all_TCGA_filt, by.x=common_col, by.y=common_col, all.x = TRUE)
PTC_transcripts_all_TCGA_clean <- PTC_transcripts_all_TCGA_clean[!duplicated(PTC_transcripts_all_TCGA_clean),]
dim(PTC_transcripts_all_TCGA_clean)
median(PTC_transcripts_all_TCGA_clean$NMD_efficiency_TPM,na.rm = TRUE)
# 2) PTCs can be repeated across samples within the same tissue (interindividual variability) or across tissues
# 2.1) Add the "duplicated within cancer" column
tmp <- ddply(PTC_transcripts_all_TCGA_clean,.(start_pos,PTC_CDS_pos,transcript_id,TCGA_cancer),nrow)
tmp <- tmp[which(tmp$V1 >= 2),]
PTCs_duplicated <- (PTC_transcripts_all_TCGA_clean$start_pos %in% tmp$start_pos) & (PTC_transcripts_all_TCGA_clean$PTC_CDS_pos %in% tmp$PTC_CDS_pos) & (PTC_transcripts_all_TCGA_clean$transcript_id %in% tmp$transcript_id) & (PTC_transcripts_all_TCGA_clean$TCGA_cancer %in% tmp$TCGA_cancer)
PTC_transcripts_all_TCGA_clean$duplicated_within_cancer <- "no"
PTC_transcripts_all_TCGA_clean[PTCs_duplicated,"duplicated_within_cancer"] <- "yes"
# 2.2) Add the "duplicated across cancer" column
tmp <- ddply(PTC_transcripts_all_TCGA_clean,.(start_pos,PTC_CDS_pos,transcript_id),nrow)
tmp <- tmp[which(tmp$V1 >= 2),]
PTCs_duplicated <- (PTC_transcripts_all_TCGA_clean$start_pos %in% tmp$start_pos) & (PTC_transcripts_all_TCGA_clean$PTC_CDS_pos %in% tmp$PTC_CDS_pos) & (PTC_transcripts_all_TCGA_clean$transcript_id %in% tmp$transcript_id)
PTC_transcripts_all_TCGA_clean$duplicated_across_cancers <- "no"
PTC_transcripts_all_TCGA_clean[PTCs_duplicated,"duplicated_across_cancers"] <- "yes"
dim(PTC_transcripts_all_TCGA_clean)
median(PTC_transcripts_all_TCGA_clean$NMD_efficiency_TPM,na.rm = TRUE)
table(PTC_transcripts_all_TCGA_clean$duplicated_within_cancer)
table(PTC_transcripts_all_TCGA_clean$duplicated_across_cancers)

# Save filtered dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_all_TCGA_confident_seq.txt")
write.table(PTC_transcripts_all_TCGA_clean, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
# PTC_transcripts_all_TCGA_clean <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)

# 3) ASE
# Merge all PTCs from all TCGA samples into one dataframe

ASE_transcripts_all_TCGA <- c()
for (TCGA_cancer in TCGA_cancers) {
    print(TCGA_cancer)
    TCGA_cancer_path <- gsub("\\[X1\\]",TCGA_cancer,TCGA_cancers_path)
    TCGA_samples <- list.files(TCGA_cancer_path)
    ASE_transcripts_cancer <- c()
    for (TCGA_sample in TCGA_samples) {
        # Check if PTC file exists
        PTC_transcripts_file_path <- paste0(TCGA_cancer_path,TCGA_sample,"/PTC_transcripts/",type,"_PTCs_ASE_transcripts_metadata.txt")
        #print(PTC_transcripts_file_path)
        if (file.exists(PTC_transcripts_file_path)) {
            # Open PTCs file
            PTC_transcripts_file <- read.table(file = PTC_transcripts_file_path, header = TRUE, sep = "\t")
            # For the time being I will remove columns related to PCA subtype, as the number is different among different cancer types
            # We don't need the nucleotide sequence either
            #PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue|fasta_sequence_wt|fasta_sequence_mut",colnames(PTC_transcripts_file))]
            PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue",colnames(PTC_transcripts_file))]
            # Add TCGA barcode as a new variable
            PTC_transcripts_file_filt$TCGA_barcode <- TCGA_sample
            # Add TCGA cancer as new variable
            PTC_transcripts_file_filt$TCGA_cancer <- TCGA_cancer
            if ( is.null(PTC_transcripts_file_filt$germline_SNV) ) {
                PTC_transcripts_file_filt$germline_SNV <- "no"
            }
            #print(dim(PTC_transcripts_file))
            # Join data frame 
            ASE_transcripts_cancer <- rbind(ASE_transcripts_cancer,PTC_transcripts_file_filt)
        }
    }
    ASE_transcripts_all_TCGA <- rbind(ASE_transcripts_all_TCGA,ASE_transcripts_cancer)
}

ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM <- -log2(ASE_transcripts_all_TCGA$altCount_ASE / ASE_transcripts_all_TCGA$refCount_ASE)

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_ASE_all_TCGA_seq.txt")
write.table(ASE_transcripts_all_TCGA, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
#ASE_transcripts_all_TCGA <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)


# 4) Perform filterings to obtain a confident set of PTCs
median(ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM,na.rm = TRUE)
# 4.1) Remove NAs in NMDeff
dim(ASE_transcripts_all_TCGA)
ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM <- ifelse(ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM == "Inf",NA,ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM)
ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM <- ifelse(ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM == "-Inf",NA,ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM)
table(is.na(ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM))
ASE_transcripts_all_TCGA_filt <- ASE_transcripts_all_TCGA[!is.na(ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM),]
dim(ASE_transcripts_all_TCGA_filt)
median(ASE_transcripts_all_TCGA_filt$ASE_NMD_efficiency_TPM,na.rm = TRUE)
# 4.2) Remove FS indel mutations that do not have a predicted downstream PTC
ASE_transcripts_all_TCGA_filt <- ASE_transcripts_all_TCGA_filt[-which(is.na(ASE_transcripts_all_TCGA_filt$PTC_CDS_pos)),]
dim(ASE_transcripts_all_TCGA_filt)
median(ASE_transcripts_all_TCGA_filt$ASE_NMD_efficiency_TPM,na.rm = TRUE)

# Merge
common_cols <- intersect(colnames(PTC_transcripts_all_TCGA_clean),colnames(ASE_transcripts_all_TCGA_filt))
PTC_ASE_transcripts_all_TCGA <- merge(PTC_transcripts_all_TCGA_clean, ASE_transcripts_all_TCGA_filt, by = common_cols, all.x = TRUE)
PTC_ASE_transcripts_all_TCGA <- PTC_ASE_transcripts_all_TCGA[!is.na(PTC_ASE_transcripts_all_TCGA$ASE_NMD_efficiency_TPM),]

# Save confident dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_ASE_all_TCGA_confident_seq.txt")
write.table(PTC_ASE_transcripts_all_TCGA, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)







################## TEMP LOEUF test ################



cols <- c("TPM_transcript_expression","TPM_transcript_expression_WT","NMD_efficiency_TPM","altCount_ASE.x","refCount_ASE.x","ASE_NMD_efficiency_TPM","LOEUF_decile","LOEUF_score","TPM_WT_samples")

PTCs_ASE_NMD_efficiencies_TCGA_filt <- PTC_ASE_transcripts_all_TCGA[PTC_ASE_transcripts_all_TCGA$altCount_ASE > 3,]
#PTCs_ASE_NMD_efficiencies_TCGA_filt2 <- PTCs_ASE_NMD_efficiencies_TCGA_filt[which(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_decile < 2),]
#PTCs_ASE_NMD_efficiencies_TCGA_filt2 <- PTCs_ASE_NMD_efficiencies_TCGA_filt2[PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM <= quantile(PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM)[3],]
#PTCs_ASE_NMD_efficiencies_TCGA_filt2[60:70,cols]

# Low LOEUF, low PTC NMDeff --> should have high ASE NMDeff than expected (what is our expectation???)
# median(PTCs_ASE_NMD_efficiencies_TCGA_filt2$ASE_NMD_efficiency_TPM)
# median(PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM)

# How many times is ASE NMDeff higher than PTC NMDeff in LOW LOEUF vs HIGH LOEUF?
PTCs_ASE_NMD_efficiencies_TCGA_filt[,c("ASE_NMD_efficiency_TPM","NMD_efficiency_TPM")] <- scale(PTCs_ASE_NMD_efficiencies_TCGA_filt[,c("ASE_NMD_efficiency_TPM","NMD_efficiency_TPM")])



for (LOEUF_decile in c(1:9)) {
    if (LOEUF_decile == 1) {print(paste0("---------- % of PTCs where ASE_NMDeff is higher than PTC_NMDeff ----------"))}
    PTCs_ASE_NMD_efficiencies_TCGA_filt2 <- PTCs_ASE_NMD_efficiencies_TCGA_filt[which(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_decile == LOEUF_decile),]
    perc <- sum(PTCs_ASE_NMD_efficiencies_TCGA_filt2$ASE_NMD_efficiency_TPM > PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM) / nrow(PTCs_ASE_NMD_efficiencies_TCGA_filt2)
    print(paste0("LOEUF decile --> ",LOEUF_decile))
    print(paste0("Percentage: ",round(perc,2)))

}
LOEUF_median <- median(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_score,na.rm=TRUE)
PTCs_ASE_NMD_efficiencies_TCGA_filt2 <- PTCs_ASE_NMD_efficiencies_TCGA_filt[which(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_decile < 2),]
sum(PTCs_ASE_NMD_efficiencies_TCGA_filt2$ASE_NMD_efficiency_TPM > PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM) / nrow(PTCs_ASE_NMD_efficiencies_TCGA_filt2)
PTCs_ASE_NMD_efficiencies_TCGA_filt3 <- PTCs_ASE_NMD_efficiencies_TCGA_filt[which(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_decile >= 2),]
sum(PTCs_ASE_NMD_efficiencies_TCGA_filt3$ASE_NMD_efficiency_TPM > PTCs_ASE_NMD_efficiencies_TCGA_filt3$NMD_efficiency_TPM) / nrow(PTCs_ASE_NMD_efficiencies_TCGA_filt3)


#cor.test(PTCs_ASE_NMD_efficiencies_TCGA_filt2$NMD_efficiency_TPM,PTCs_ASE_NMD_efficiencies_TCGA_filt2$ASE_NMD_efficiency_TPM)



# ME sale positiva la correlacon, pero tengo que mirar que la diferencia entre NMDeffs sea muy grande en este caso en particular debido al dosage compensatin.



filter <- which(PTCs_ASE_NMD_efficiencies_TCGA_filt$NMD_efficiency_TPM < 0 & PTCs_ASE_NMD_efficiencies_TCGA_filt$ASE_NMD_efficiency_TPM > 0)
df <- PTCs_ASE_NMD_efficiencies_TCGA_filt[filter,]
summary(df$LOEUF_score)

df <- PTCs_ASE_NMD_efficiencies_TCGA_filt[which(PTCs_ASE_NMD_efficiencies_TCGA_filt$LOEUF_decile < 3),]
df[1:10,]


