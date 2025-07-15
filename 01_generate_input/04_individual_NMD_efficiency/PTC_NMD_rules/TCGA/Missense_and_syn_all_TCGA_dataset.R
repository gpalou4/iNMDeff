library(dplyr)
library(plyr)

type <- "somatic"

TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/[X1]/"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1

# 1) Merge all MisSyn from all TCGA samples into one dataframe

MisSyn_transcripts_all_TCGA <- c()
for (TCGA_cancer in TCGA_cancers) {
    print(TCGA_cancer)
    TCGA_cancer_path <- gsub("\\[X1\\]",TCGA_cancer,TCGA_cancers_path)
    TCGA_samples <- list.files(TCGA_cancer_path)
    MisSyn_transcripts_cancer <- c()
    for (TCGA_sample in TCGA_samples) {
        # Check if MisSyn file exists
        MisSyn_transcripts_file_path <- paste0(TCGA_cancer_path,TCGA_sample,"/PTC_transcripts/",type,"_missense_syn_transcripts_metadata.txt")
        #print(MisSyn_transcripts_file_path)
        if (file.exists(MisSyn_transcripts_file_path)) {
            # Open MisSyns file
            MisSyn_transcripts_file <- read.table(file = MisSyn_transcripts_file_path, header = TRUE, sep = "\t")
            # For the time being I will remove columns related to PCA subtype, as the number is different among different cancer types
            # We don't need the nucleotide sequence either
            MisSyn_transcripts_file_filt <- MisSyn_transcripts_file[,-grep("subtissue|fasta_sequence_wt|fasta_sequence_mut",colnames(MisSyn_transcripts_file))]
            #MisSyn_transcripts_file_filt <- MisSyn_transcripts_file[,-grep("subtissue",colnames(MisSyn_transcripts_file))]
            # Add TCGA barcode as a new variable
            MisSyn_transcripts_file_filt$TCGA_barcode <- TCGA_sample
            # Add TCGA cancer as new variable
            MisSyn_transcripts_file_filt$TCGA_cancer <- TCGA_cancer
            if ( is.null(MisSyn_transcripts_file_filt$germline_SNV) ) {
                MisSyn_transcripts_file_filt$germline_SNV <- "no"
            }
            #print(dim(MisSyn_transcripts_file))
            # Join data frame 
            MisSyn_transcripts_cancer <- rbind(MisSyn_transcripts_cancer,MisSyn_transcripts_file_filt)
        }
    }
    MisSyn_transcripts_all_TCGA <- rbind(MisSyn_transcripts_all_TCGA,MisSyn_transcripts_cancer)
}

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/",type,"_MisSyn_all_TCGA.txt")
write.table(MisSyn_transcripts_all_TCGA, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)

