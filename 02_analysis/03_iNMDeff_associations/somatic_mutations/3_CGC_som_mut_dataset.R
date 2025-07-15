# Open files

TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/[X1]/"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)
common_col_merge <- c("Gene.Symbol","Synonyms")
common_col <- colnames(CGC)[-which(colnames(CGC) %in% common_col_merge)]
mutations_non_syn <- c("stopgain","nonsynonymous SNV","startloss","stoploss","frameshift deletion",
                "frameshift insertion","nonframeshift deletion","nonframeshift insertion","splicing")
mutations_non_syn <- paste(mutations_non_syn,collapse="|")

# 1) Merge all CGC files from all TCGA samples into one dataframe
cols <- c("Gene.Symbol","Synonyms","somatic_mut_transcript","somatic_CNV")
CGC_all_cancers <- c()
for (TCGA_cancer in TCGA_cancers) {
    print(TCGA_cancer)
    TCGA_cancer_path <- gsub("\\[X1\\]",gsub("TCGA-","",TCGA_cancer),TCGA_cancers_path)
    TCGA_samples <- list.files(TCGA_cancer_path)
    CGC_cancer <- CGC
    for (TCGA_sample in TCGA_samples) {
        # Check if CGC file exists
        CGC_somatic_mut_CNV_file_path <- paste0(TCGA_cancer_path,TCGA_sample,"/",TCGA_sample,"_CGC_somatic_mut_CNV.txt")
        if (file.exists(CGC_somatic_mut_CNV_file_path)) {
            # Open sample CGC file
            CGC_somatic_mut_CNV <- read.csv(file = CGC_somatic_mut_CNV_file_path, header = TRUE, sep = "\t")
            CGC_somatic_mut_CNV$somatic_mut_synonymous <- "no"
            # Remove columns
            CGC_somatic_mut_CNV <- CGC_somatic_mut_CNV[,!colnames(CGC_somatic_mut_CNV)%in%common_col]
            # 1) Check if sample has somatic mutations
            hasSomaticMut <- grep("somatic_mut_transcript",colnames(CGC_somatic_mut_CNV))
            hasSomaticMutTruncating <- grep("somatic_mut_truncating",colnames(CGC_somatic_mut_CNV))
            if (length(hasSomaticMutTruncating) == 0) {
                CGC_somatic_mut_CNV$somatic_mut_truncating <- "no"
            }
            hasSomaticMutMissense <- grep("somatic_mut_missense",colnames(CGC_somatic_mut_CNV))
            if (length(hasSomaticMutMissense) == 0) {
                CGC_somatic_mut_CNV$somatic_mut_missense <- "no"
            }
            if (length(hasSomaticMut) == 1) {
                # Filters
                syn_index <- grep("synonymous SNV",CGC_somatic_mut_CNV$somatic_variant_classification)
                syn_mut <- CGC_somatic_mut_CNV[syn_index,]
                syn_mut_only <- syn_mut[grep(mutations_non_syn,syn_mut$somatic_variant_classification, invert = TRUE),]
                syn_mut_only <- syn_mut_only[is.na(syn_mut_only$somatic_CNV),]
                if (nrow(syn_mut_only) >= 1) {
                    syn_index <- which(CGC_somatic_mut_CNV$Gene.Symbol %in% as.character(syn_mut_only$Gene.Symbol))
                    #CGC_somatic_mut_CNV_filt <- CGC_somatic_mut_CNV[-syn_index,]
                    # Add a new column for Synonymous mut, separately
                    CGC_somatic_mut_CNV[syn_index,"somatic_mut_synonymous"] <- "yes"
                }
                CGC_somatic_mut_CNV_filt <- CGC_somatic_mut_CNV
                # Remove unnecessary cols
                CGC_somatic_mut_CNV_filt[,c("somatic_filter","somatic_EVS","somatic_VAF","somatic_mut_transcript","somatic_variant_classification")] <- NULL
                if (nrow(CGC_somatic_mut_CNV_filt) >= 1 ){
                    hasSomaticMut <- 1
                } else {
                    hasSomaticMut <- c()
                }
            } else {
                CGC_somatic_mut_CNV_filt <- CGC_somatic_mut_CNV
            }
            # Merge
            CGC_cancer <- merge(CGC_cancer, CGC_somatic_mut_CNV_filt, by.x = common_col_merge, by.y = common_col_merge, all.x = TRUE)
            # 2) Check if sample has somatic CNV
            hasSomaticCNV <- grep("somatic_CNV",colnames(CGC_somatic_mut_CNV))
            # Modify somatic_CNV amd somatic_mut_missense somatic_mut_truncating columns
            if (length(hasSomaticCNV) == 1) {
                CGC_cancer$somatic_CNV <- ifelse(is.na(as.character(CGC_cancer$somatic_CNV)),"no",as.character(CGC_cancer$somatic_CNV))
            } 
            if (length(hasSomaticMut) == 1) {
                CGC_cancer$somatic_mut_missense <- ifelse(is.na(as.character(CGC_cancer$somatic_mut_missense)),"no",as.character(CGC_cancer$somatic_mut_missense))
                CGC_cancer$somatic_mut_truncating <- ifelse(is.na(as.character(CGC_cancer$somatic_mut_truncating)),"no",as.character(CGC_cancer$somatic_mut_truncating))
                CGC_cancer$somatic_mut_synonymous <- ifelse(is.na(as.character(CGC_cancer$somatic_mut_synonymous)),"no",as.character(CGC_cancer$somatic_mut_synonymous))
            }
            colnames(CGC_cancer)[colnames(CGC_cancer) == "somatic_CNV"] <- paste0(TCGA_sample,"_somatic_CNV")
            colnames(CGC_cancer)[colnames(CGC_cancer) == "somatic_mut_missense"] <- paste0(TCGA_sample,"_somatic_mut_missense")
            colnames(CGC_cancer)[colnames(CGC_cancer) == "somatic_mut_truncating"] <- paste0(TCGA_sample,"_somatic_mut_truncating")
            colnames(CGC_cancer)[colnames(CGC_cancer) == "somatic_mut_synonymous"] <- paste0(TCGA_sample,"_somatic_mut_synonymous")
        }
    }
    # Save CGC for the cancer
    output_path <- paste0(TCGA_cancer_path,TCGA_cancer,"_CGC_somatic_mut_CNV.txt")
    write.table(CGC_cancer, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
    # Merge with rest of cancers
    if (length(CGC_all_cancers) == 0 ) {
        CGC_all_cancers <- CGC_cancer
    } else {
        CGC_cancer <- CGC_cancer[,!colnames(CGC_cancer)%in%common_col]
        CGC_all_cancers <- merge(CGC_all_cancers, CGC_cancer, by.x = common_col_merge, by.y = common_col_merge, all.x = TRUE)
    }
}

# Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/TCGA_CGC_somatic_mut_CNV.txt")
write.table(CGC_all_cancers, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
