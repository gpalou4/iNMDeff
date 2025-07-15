library(dplyr)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
type <- args[1]

GTEx_tissue_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt"
GTEx_tissues_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/[X1]/"
GTEx_tissues <- read.table(file = GTEx_tissue_names_path, stringsAsFactors = FALSE)$V1

# 1) Merge all PTCs from all GTEx samples into one dataframe

PTC_transcripts_all_GTEx <- c()
for (GTEx_tissue in GTEx_tissues) {
    print(GTEx_tissue)
    GTEx_tissue_path <- gsub("\\[X1\\]",GTEx_tissue,GTEx_tissues_path)
    GTEx_samples <- list.files(GTEx_tissue_path)
    GTEx_samples <- GTEx_samples[grep("GTEX-\\w*_germline_PTC_transcripts_metadata.txt",GTEx_samples)]
    GTEx_samples <- gsub("(GTEX\\-\\w{4,5})\\_.*","\\1",GTEx_samples)
    PTC_transcripts_tissue <- c()
    for (GTEx_sample in GTEx_samples) {
        # Check if PTC file exists
        PTC_transcripts_file_path <- paste0(GTEx_tissue_path,GTEx_sample,"_",type,"_PTC_transcripts_metadata.txt")
        #print(PTC_transcripts_file_path)
        if (file.exists(PTC_transcripts_file_path)) {
            # Open PTCs file
            PTC_transcripts_file <- read.table(file = PTC_transcripts_file_path, header = TRUE, sep = "\t")
            # For the time being I will remove columns related to PCA subtype, as the number is different among different tissue types
            # We don't need the nucleotide sequence either
            #PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue|fasta_sequence_wt|fasta_sequence_mut",colnames(PTC_transcripts_file))]
            remove <- grep("subtissue",colnames(PTC_transcripts_file))
            if (length(remove) != 0) {
                PTC_transcripts_file_filt <- PTC_transcripts_file[,-remove]
            } else {
                PTC_transcripts_file_filt <- PTC_transcripts_file
            }
            # Add GTEx barcode as a new variable
            PTC_transcripts_file_filt$GTEx_sample <- GTEx_sample
            # Add GTEx tissue as new variable
            PTC_transcripts_file_filt$GTEx_tissue <- GTEx_tissue
            if ( is.null(PTC_transcripts_file_filt$germline_SNV) ) {
                PTC_transcripts_file_filt$germline_SNV <- "no"
            }
            #print(dim(PTC_transcripts_file))
            # Join data frame 
            PTC_transcripts_tissue <- rbind(PTC_transcripts_tissue,PTC_transcripts_file_filt)
        }
    }
    # Fix tissues with <4 PCs
    if (length(grep("tissue_PC*",colnames(PTC_transcripts_tissue))) == 3) {
        PTC_transcripts_tissue$tissue_PC4 <- NA
    }
    PTC_transcripts_all_GTEx <- rbind(PTC_transcripts_all_GTEx,PTC_transcripts_tissue)
}

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_all_GTEx_seq.txt")
write.table(PTC_transcripts_all_GTEx, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
#PTC_transcripts_all_GTEx <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)

# 2) Perform filterings to obtain a confident set of PTCs
median(PTC_transcripts_all_GTEx$NMD_efficiency_TPM,na.rm = TRUE)
# 2.1) Remove NAs in NMDeff
dim(PTC_transcripts_all_GTEx)
PTC_transcripts_all_GTEx$NMD_efficiency_TPM <- ifelse(PTC_transcripts_all_GTEx$NMD_efficiency_TPM %in% c("Inf","-Inf"),NA,PTC_transcripts_all_GTEx$NMD_efficiency_TPM)
table(is.na(PTC_transcripts_all_GTEx$NMD_efficiency_TPM))
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx[!is.na(PTC_transcripts_all_GTEx$NMD_efficiency_TPM),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.2) Remove FS indel mutations that do not have a predicted downstream PTC
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[-which(is.na(PTC_transcripts_all_GTEx_filt$PTC_CDS_pos)),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.2) Ovelapping additional germline nonsense, FS or splice site variants
# PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$germline_SNV == "no"),]
# dim(PTC_transcripts_all_GTEx_filt)
# median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.4) Genes with 1 exon and without splice site at 3'UTR
filter <- which(PTC_transcripts_all_GTEx_filt$transcript_CDS_exon_num == 1 & PTC_transcripts_all_GTEx_filt$splice_site_3UTR == "no")
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[-filter,]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.5) Positive Selected Genes (PSG)
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$PSG == "no"),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.6) At least 3 samples per WT-transcript per RNAseq NMF cluster
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$TPM_WT_samples >= 3),]
print(dim(PTC_transcripts_all_GTEx_filt))
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.7) Coefficient of variation (mean/SD) > 0.5 --> This reduces a lot the NMDeff values (from 1.6 raw to 0.3)
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$coeff_var <= 2),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.8) Remove PTCs located in NMD-evading regions
# 1) last exon or 55 nt
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[!PTC_transcripts_all_GTEx_filt$X55_nt_last_exon == "NMD-evading",]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2) First 200nt
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[PTC_transcripts_all_GTEx_filt$TSS_PTC_dist > 200,]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.9) LOEUF score (first decile are the most neg selected ones)
# We don't remove NA's tho
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[-which(PTC_transcripts_all_GTEx_filt$LOEUF_decile <= 2),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.10) Keep only Heterozygous variants
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$genotype %in% c("0/1","0|1","1|0")),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.11) Low expressed transcripts (<5 TPM)
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt[which(PTC_transcripts_all_GTEx_filt$median_TPM_exp_transcript >= 5),]
dim(PTC_transcripts_all_GTEx_filt)
median(PTC_transcripts_all_GTEx_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.12) Add turn_over datasets
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome.csv"
half_life_1 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_2.csv"
half_life_2 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_transcriptome_bcells.csv"
half_life_3 <- read.table(file = file, header = TRUE, sep = ",")
file <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/halflifes_PCA.csv"
half_life_4 <- read.table(file = file, header = TRUE, sep = ",")
# merge
vec1 <- unlist(lapply(strsplit(PTC_transcripts_all_GTEx_filt$gene_id,";"), function(gene){gene[1]}))
vec2 <- unlist(lapply(strsplit(PTC_transcripts_all_GTEx_filt$gene_id,";"), function(gene){gene[2]}))
PTC_transcripts_all_GTEx_filt$gene_id_1 <- vec1
PTC_transcripts_all_GTEx_filt$gene_id_2 <- vec2
PTC_transcripts_all_GTEx_filt <- merge(PTC_transcripts_all_GTEx_filt,half_life_4[,c("Ensembl.Gene.Id","half.life..PC1.")], by.x ="gene_id_1", by.y ="Ensembl.Gene.Id", all.x = TRUE)
PTC_transcripts_all_GTEx_filt <- merge(PTC_transcripts_all_GTEx_filt,half_life_4[,c("Ensembl.Gene.Id","half.life..PC1.")], by.x ="gene_id_2", by.y ="Ensembl.Gene.Id", all.x = TRUE)
PTC_transcripts_all_GTEx_filt <- PTC_transcripts_all_GTEx_filt %>% 
    mutate(half_life = coalesce(half.life..PC1..x,half.life..PC1..y))
PTC_transcripts_all_GTEx_filt[,c("half.life..PC1..x","half.life..PC1..y")] <- NULL

# 2.12) Remove repeated PTCs
# 1) PTCs can overlap different transcripts within the same sample
# For each same PTC start position and SAME sample let's take one transcript
ensembl_transcripts_MANE_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_MANE_v107.txt"
ensembl_transcripts_MANE <- read.table(file = ensembl_transcripts_MANE_path, header = TRUE)
# First the MANE Select if possible, if not, then randomly
ensembl_transcripts_MANE$MANE_main_isoform <- "yes"
PTC_transcripts_all_GTEx_filt <- merge(PTC_transcripts_all_GTEx_filt,ensembl_transcripts_MANE[,c("ensembl_transcript_id","MANE_main_isoform")],
        by.x="transcript_id", by.y ="ensembl_transcript_id", all.x = TRUE)
PTC_transcripts_all_GTEx_filt$MANE_main_isoform <- ifelse(is.na(PTC_transcripts_all_GTEx_filt$MANE_main_isoform),"no","yes")
PTC_transcripts_all_GTEx_clean <- aggregate(NMD_efficiency_TPM ~ GTEx_tissue + GTEx_sample + PTC_CDS_pos + start_pos + MANE_main_isoform, data = PTC_transcripts_all_GTEx_filt, FUN=function(x) x[sample(seq_along(x), 1)])
# Remove those transcripts with both MANE + random (keep onle the MANE ones) ~ 664 transcripts
PTC_transcripts_all_GTEx_clean <- PTC_transcripts_all_GTEx_clean[-which(duplicated(PTC_transcripts_all_GTEx_clean[,1:4], fromLast = TRUE)),]
# Merge again
common_col <- c("GTEx_tissue","GTEx_sample","NMD_efficiency_TPM","PTC_CDS_pos","start_pos","MANE_main_isoform")
PTC_transcripts_all_GTEx_clean <- merge(PTC_transcripts_all_GTEx_clean,PTC_transcripts_all_GTEx_filt, by.x=common_col, by.y=common_col, all.x = TRUE)
PTC_transcripts_all_GTEx_clean <- PTC_transcripts_all_GTEx_clean[!duplicated(PTC_transcripts_all_GTEx_clean),]
dim(PTC_transcripts_all_GTEx_clean)
median(PTC_transcripts_all_GTEx_clean$NMD_efficiency_TPM,na.rm = TRUE)

# 2) PTCs can be repeated across samples within the same tissue (interindividual variability, we need to correct it for this analysis)
tmp <- ddply(PTC_transcripts_all_GTEx_clean,.(start_pos,PTC_CDS_pos,transcript_id),nrow)
tmp <- tmp[which(tmp$V1 <= 5),]
PTCs_to_keep <- (PTC_transcripts_all_GTEx_clean$start_pos %in% tmp$start_pos) & (PTC_transcripts_all_GTEx_clean$PTC_CDS_pos %in% tmp$PTC_CDS_pos) & (PTC_transcripts_all_GTEx_clean$transcript_id %in% tmp$transcript_id)
PTC_transcripts_all_GTEx_clean <- PTC_transcripts_all_GTEx_clean[PTCs_to_keep,]
dim(PTC_transcripts_all_GTEx_clean)
median(PTC_transcripts_all_GTEx_clean$NMD_efficiency_TPM,na.rm = TRUE)

# Save filtered dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_all_GTEx_confident_seq.txt")
write.table(PTC_transcripts_all_GTEx_clean, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)

# 3) ASE
# Merge all PTCs from all GTEx samples into one dataframe

ASE_transcripts_all_GTEx <- c()
for (GTEx_tissue in GTEx_tissues) {
    print(GTEx_tissue)
    GTEx_tissue_path <- gsub("\\[X1\\]",GTEx_tissue,GTEx_tissues_path)
    GTEx_file_names <- list.files(GTEx_tissue_path)
    GTEx_samples <- gsub("(.*)_germline.*","\\1",GTEx_file_names)
    ASE_transcripts_tissue <- c()
    for (GTEx_sample in GTEx_samples) {
        # Check if PTC file exists
        PTC_transcripts_file_path <- paste0(GTEx_tissue_path,"/",GTEx_sample,"_",type,"_PTC_ASE_transcripts_metadata.txt")
        #print(PTC_transcripts_file_path)
        if (file.exists(PTC_transcripts_file_path)) {
            # Open PTCs file
            PTC_transcripts_file <- read.table(file = PTC_transcripts_file_path, header = TRUE, sep = "\t")
            # For the time being I will remove columns related to PCA subtype, as the number is different among different tissue types
            # We don't need the nucleotide sequence either
            #PTC_transcripts_file_filt <- PTC_transcripts_file[,-grep("subtissue|fasta_sequence_wt|fasta_sequence_mut",colnames(PTC_transcripts_file))]
            # Add GTEx barcode as a new variable
            PTC_transcripts_file$GTEx_sample <- GTEx_sample
            # Add GTEx tissue as new variable
            PTC_transcripts_file$GTEx_tissue <- GTEx_tissue
            if ( is.null(PTC_transcripts_file$germline_SNV) ) {
                PTC_transcripts_file$germline_SNV <- "no"
            }
            #print(dim(PTC_transcripts_file))
            # Join data frame 
            ASE_transcripts_tissue <- rbind(ASE_transcripts_tissue,PTC_transcripts_file)
        }
    }
    ASE_transcripts_all_GTEx <- rbind(ASE_transcripts_all_GTEx,ASE_transcripts_tissue)
}

ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM <- -log2(ASE_transcripts_all_GTEx$altCount_ASE / ASE_transcripts_all_GTEx$refCount_ASE)

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_ASE_all_GTEx_seq.txt")
write.table(ASE_transcripts_all_GTEx, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
#ASE_transcripts_all_GTEx <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)



# 4) Perform filterings to obtain a confident set of PTCs
median(ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM,na.rm = TRUE)
# 4.1) Remove NAs in NMDeff
dim(ASE_transcripts_all_GTEx)
ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM <- ifelse(ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM == "Inf",NA,ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM)
ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM <- ifelse(ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM == "-Inf",NA,ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM)
table(is.na(ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM))
ASE_transcripts_all_GTEx_filt <- ASE_transcripts_all_GTEx[!is.na(ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM),]
dim(ASE_transcripts_all_GTEx_filt)
median(ASE_transcripts_all_GTEx_filt$ASE_NMD_efficiency_TPM,na.rm = TRUE)
# 4.2) Remove FS indel mutations that do not have a predicted downstream PTC
ASE_transcripts_all_GTEx_filt <- ASE_transcripts_all_GTEx_filt[-which(is.na(ASE_transcripts_all_GTEx_filt$PTC_CDS_pos)),]
dim(ASE_transcripts_all_GTEx_filt)
median(ASE_transcripts_all_GTEx_filt$ASE_NMD_efficiency_TPM,na.rm = TRUE)

# Merge
PTC_ASE_transcripts_all_GTEx <- merge(PTC_transcripts_all_GTEx_clean, ASE_transcripts_all_GTEx_filt, by = colnames(ASE_transcripts_all_GTEx_filt)[c(1:60,64,65)], all.x = TRUE)
PTC_ASE_transcripts_all_GTEx <- PTC_ASE_transcripts_all_GTEx[!is.na(PTC_ASE_transcripts_all_GTEx$ASE_NMD_efficiency_TPM),]

# Save confident dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_PTCs_ASE_all_GTEx_confident.txt")
write.table(PTC_ASE_transcripts_all_GTEx, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)