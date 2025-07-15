#conda activate /home/gpalou/anaconda3_envs/general

# Enter in fsupeksvr2 (it crushes elsewhere)
library(tidyverse)

# 1) Open NMDeff table
type <- "germline"
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_TCGA_dataset/",type,"_PTCs_all_TCGA.txt")
PTC_transcripts_all_TCGA <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors =FALSE)
# 2) Filterings
median(PTC_transcripts_all_TCGA$NMD_efficiency_TPM,na.rm = TRUE)
# 2.1) Remove NAs in NMDeff
dim(PTC_transcripts_all_TCGA)
PTC_transcripts_all_TCGA$NMD_efficiency_TPM <- ifelse(PTC_transcripts_all_TCGA$NMD_efficiency_TPM == "Inf",NA,PTC_transcripts_all_TCGA$NMD_efficiency_TPM)
PTC_transcripts_all_TCGA$NMD_efficiency_TPM <- ifelse(PTC_transcripts_all_TCGA$NMD_efficiency_TPM == "-Inf",NA,PTC_transcripts_all_TCGA$NMD_efficiency_TPM)
table(is.na(PTC_transcripts_all_TCGA$NMD_efficiency_TPM))
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA[!is.na(PTC_transcripts_all_TCGA$NMD_efficiency_TPM),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)
# 2.2) Remove FS indel mutations that do not have a predicted downstream PTC
PTC_transcripts_all_TCGA_filt <- PTC_transcripts_all_TCGA_filt[-which(is.na(PTC_transcripts_all_TCGA_filt$PTC_CDS_pos)),]
dim(PTC_transcripts_all_TCGA_filt)
median(PTC_transcripts_all_TCGA_filt$NMD_efficiency_TPM,na.rm = TRUE)

# 3) Difference between NMDeff NMD-triggering vs NMD-evading
NMDeff_PTCs_transcript <- aggregate(NMD_efficiency_TPM ~ X55_nt_last_exon + transcript_id, data = PTC_transcripts_all_TCGA_filt, median)
NMDeff_PTCs_transcript$X55_nt_last_exon <- ifelse(NMDeff_PTCs_transcript$X55_nt_last_exon == "NMD-triggering","NMD_triggering","NMD_evading")
PTC_transcripts_NMDeff_diff <- NMDeff_PTCs_transcript %>%
  pivot_wider(names_from = X55_nt_last_exon,
              values_from = NMD_efficiency_TPM,
              id_cols = transcript_id) %>%
  mutate(NMDdiff = NMD_triggering - NMD_evading) %>%
  arrange(NMDdiff)
# 4) Transcripts with bad NMDeff ratio
# filter <- which(PTC_transcripts_NMDeff_diff$NMDdiff <= 0)
# transcripts_to_filter <- PTC_transcripts_NMDeff_diff[filter,"transcript_id"]
# Save
output_path <- "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/TCGA_PTC_transcripts_NMDeff_diff.txt"
write.table(PTC_transcripts_NMDeff_diff, file = output_path, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
