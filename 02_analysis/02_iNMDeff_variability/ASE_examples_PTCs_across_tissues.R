library(ggplot2)
library(dplyr)

TCGA_ASE_germline_PTCs <- function() {

    # Obtain ASE file
    file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_ASE_all_TCGA_seq.txt")
    TCGA_ASE_PTCs <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Filters
    
    # 5) Filterings for the NB regression
    # 5.1) Coverage 5 for nonsense and 2 for FS (filter for ASE PTCs)
    filter_coverage_nonsense <- which(TCGA_ASE_PTCs$totalCount_ASE >= 5 & TCGA_ASE_PTCs$stopgain == "nonsense")
    filter_coverage_FS <- which(TCGA_ASE_PTCs$totalCount_ASE >= 2 & TCGA_ASE_PTCs$stopgain != "nonsense")
    TCGA_ASE_PTCs_filt <- TCGA_ASE_PTCs[c(filter_coverage_nonsense,filter_coverage_FS),]
    # 5.2) Remove ASE-PTCs from FS with no predicted PTC
    if (length(filter_coverage_FS) != 0) {
        filter_FS_no_PTC <- which(TCGA_ASE_PTCs_filt$stopgain != "nonsense" & !is.na(TCGA_ASE_PTCs_filt$PTC_CDS_pos)) 
    } else {
        filter_FS_no_PTC <- NULL
    }
    # 5.3) Remove ASE-PTCs with specified VAF
    filter_VAF <- which(TCGA_ASE_PTCs_filt$VAF > 0.2)
    # 5.4) Remove transcripts with 1 exon
    filter_single_exon <- which(TCGA_ASE_PTCs_filt$transcript_CDS_exon_num == 1 & TCGA_ASE_PTCs_filt$splice_site_3UTR == "no")
    # 5.5) LOEUF score 1
    LOEUF_score <- 0
    filter_LOEUF <- which(TCGA_ASE_PTCs_filt$LOEUF_decile <= LOEUF_score)
    # 5.6) PSG genes
    filter_PSG <- which(TCGA_ASE_PTCs_filt$PSG == "yes")
    # 5.6) Transcripts overlapping somatic mut/CNV
    filter_somatic_mut_CNV <- which(TCGA_ASE_PTCs_filt$somatic_CNV_SNV == "yes")
    # 5.8) NMD-triggering vs NMD-evading
    filter_NMDtype <- which(
        TCGA_ASE_PTCs_filt$X55_nt_last_exon == "NMD-evading" |
        TCGA_ASE_PTCs_filt$TSS_PTC_dist <= 200
            )
    # 5.9) Non protein-coding
    filter_nonCoding <- which(TCGA_ASE_PTCs_filt$protein_coding == "no")
    # 5.10) Overlapping germline mutations...
    # Add up filters
    all_filters <- unique(c(filter_FS_no_PTC,filter_VAF,filter_single_exon,filter_LOEUF,filter_somatic_mut_CNV,
                    filter_NMDtype,filter_nonCoding,filter_PSG, filter_homVariants))
    TCGA_ASE_PTCs_filt <- TCGA_ASE_PTCs_filt[-all_filters,]

    NMD_final_df <- TCGA_ASE_PTCs_filt[,c("stopgain","variant_type","genotype","start_pos","gene_id",
                "refCount_ASE","altCount_ASE","totalCount_ASE","VAF","TCGA_barcode","TCGA_cancer")]
    colnames(NMD_final_df) <- c("stopgain","variant_type","genotype","start_pos","gene_id","refCount","altCount","totalCount","MAF",
                                "TCGA_barcode","TCGA_cancer")
    # Remove redundant transcripts with same ASE counts
    NMD_final_df <- NMD_final_df[which(!duplicated(NMD_final_df)),]
    # Return
    return(NMD_final_df)
}

GTEx_ASE_germline_PTCs <- function() {

    # Obtain ASE file
    file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_ASE_all_GTEx_seq.txt")
    GTEx_ASE_PTCs <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Filters
    
    # 5) Filterings for the NB regression
    # 5.1) Coverage 5 for nonsense and 2 for FS (filter for ASE PTCs)
    filter_coverage_nonsense <- which(GTEx_ASE_PTCs$totalCount_ASE >= 5 & GTEx_ASE_PTCs$stopgain == "nonsense")
    filter_coverage_FS <- which(GTEx_ASE_PTCs$totalCount_ASE >= 2 & GTEx_ASE_PTCs$stopgain != "nonsense")
    GTEx_ASE_PTCs_filt <- GTEx_ASE_PTCs[c(filter_coverage_nonsense,filter_coverage_FS),]
    # 5.2) Remove ASE-PTCs from FS with no predicted PTC
    if (length(filter_coverage_FS) != 0) {
        filter_FS_no_PTC <- which(GTEx_ASE_PTCs_filt$stopgain != "nonsense" & !is.na(GTEx_ASE_PTCs_filt$PTC_CDS_pos)) 
    } else {
        filter_FS_no_PTC <- NULL
    }
    # 5.3) Remove ASE-PTCs with specified VAF
    filter_VAF <- which(GTEx_ASE_PTCs_filt$VAF > 0.2)
    # 5.4) Remove transcripts with 1 exon
    filter_single_exon <- which(GTEx_ASE_PTCs_filt$transcript_CDS_exon_num == 1 & GTEx_ASE_PTCs_filt$splice_site_3UTR == "no")
    # 5.5) LOEUF score 1
    LOEUF_score <- 0
    filter_LOEUF <- which(GTEx_ASE_PTCs_filt$LOEUF_decile <= LOEUF_score)
    # 5.6) PSG genes
    filter_PSG <- which(GTEx_ASE_PTCs_filt$PSG == "yes")
    # 5.6) Transcripts overlapping somatic mut/CNV
    # filter_somatic_mut_CNV <- which(GTEx_ASE_PTCs_filt$somatic_CNV_SNV == "yes")
    # 5.8) NMD-triggering vs NMD-evading
    filter_NMDtype <- which(
        GTEx_ASE_PTCs_filt$X55_nt_last_exon == "NMD-evading" |
        GTEx_ASE_PTCs_filt$TSS_PTC_dist <= 200
            )
    # 5.9) Non protein-coding
    filter_nonCoding <- which(GTEx_ASE_PTCs_filt$protein_coding == "no")
    # 5.10) Overlapping germline mutations...
    # 5.11) Take only Het variants
    filter_homVariants <- which(!GTEx_ASE_PTCs_filt$genotype %in% c("0/1","1|0","0|1"))
    # Add up filters
    all_filters <- unique(c(filter_FS_no_PTC,filter_VAF,filter_single_exon,filter_LOEUF,
                    filter_NMDtype,filter_nonCoding,filter_PSG, filter_homVariants))
    GTEx_ASE_PTCs_filt <- GTEx_ASE_PTCs_filt[-all_filters,]

    NMD_final_df <- GTEx_ASE_PTCs_filt[,c("stopgain","variant_type","genotype","start_pos","gene_id",
                "refCount_ASE","altCount_ASE","totalCount_ASE","VAF","acronyms","GTEx_sample")]
    colnames(NMD_final_df) <- c("stopgain","variant_type","genotype","start_pos","gene_id","refCount","altCount","totalCount","MAF",
                                "GTEx_tissue","GTEx_sample")
    # Remove redundant transcripts with same ASE counts
    NMD_final_df <- NMD_final_df[which(!duplicated(NMD_final_df)),]
    # Return
    return(NMD_final_df)
}

library("ggplot2")
library("dplyr")
library("viridis")

# 1) Data

# 1.1) TCGA sample NMDeff
# endogenous_NMD_genesets <-  c("NMD Colombo","NMD Karousis","NMD Tani","NMD Courtney","NMD Ensembl",
#                       "NMD All","NMD Consensus","NMD SMG6","NMD SMG7",
#                       "RandomGenes without NMD features","RandomGenes with NMD features")
# ASE_NMD_genesets <- c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01","Synonymous 0.01",
#                       "PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2")

# # PTC // ASE // Endogenosample_NMD_efficiencies_TCGAus
# sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
# sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# 2) Select recurrent PTCs examples from GTEx data

# 2.1) GTEx
GTEx_ASE_PTCs_df <- GTEx_ASE_germline_PTCs()
dim(GTEx_ASE_PTCs_df)
# Calculate ASE (ALT / (REF + ALT) * 100)
GTEx_ASE_PTCs_df$ASE_percentage <- (GTEx_ASE_PTCs_df$altCount / 
                                   (GTEx_ASE_PTCs_df$refCount + GTEx_ASE_PTCs_df$altCount)) * 100
# 2.2) TCGA
TCGA_ASE_PTCs_df <- TCGA_ASE_germline_PTCs()
dim(TCGA_ASE_PTCs_df)
# Calculate ASE (ALT / (REF + ALT) * 100)
TCGA_ASE_PTCs_df$ASE_percentage <- (TCGA_ASE_PTCs_df$altCount / 
                                   (TCGA_ASE_PTCs_df$refCount + TCGA_ASE_PTCs_df$altCount)) * 100

# 2.3) Obtain a PTC that is common between TCGA and GTEx

# pos <- rev(sort(table(TCGA_ASE_PTCs_df$start_pos)))[44]
TCGA_pos <- rev(sort(table(TCGA_ASE_PTCs_df$start_pos)))[1:200]
GTEx_pos <- rev(sort(table(GTEx_ASE_PTCs_df$start_pos)))[1:200]
common_pos <- names(TCGA_pos)[names(TCGA_pos) %in% names(GTEx_pos)][75]

# 2.3.1) Plot GTEx tissues 
# 73504978
GTEx_ASE_PTCs_df_filt <- GTEx_ASE_PTCs_df[which(GTEx_ASE_PTCs_df$start_pos %in%  common_pos),]
# GTEx_ASE_PTCs_df_filt

# Order tissues by median ASE
GTEx_ASE_PTCs_df_filt <- GTEx_ASE_PTCs_df_filt %>%
  group_by(GTEx_tissue) %>%
  mutate(median_ASE = median(ASE_percentage, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(median_ASE)) %>%
  mutate(GTEx_tissue = factor(GTEx_tissue, levels = unique(GTEx_tissue)))

# Create a palette with 54 distinct colors
num_tissues <- length(unique(GTEx_ASE_PTCs_df_filt$GTEx_tissue))
color_palette <- viridis::viridis(num_tissues, option = "plasma")

# Create the plot
plot <- ggplot(GTEx_ASE_PTCs_df_filt, aes(x = GTEx_tissue, y = ASE_percentage, color = GTEx_tissue)) +
    geom_boxplot(outlier.shape = NA) +  # Avoid overlapping outliers with points
    geom_jitter(aes(alpha = 0.5), width = 0.2, size = 2) +  # Add points with alpha
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +  # Add 50% line
    theme_classic() +
    labs(
    x = "GTEx Tissue",
    y = "ASE (%)",
    title = paste0("ASE for PTC variant in position ", common_pos)
    ) +
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100)) +  # Set Y limits to 100
#   scale_color_brewer(palette = "Set3")  # Use a color palette for 14 tissues  
    scale_color_manual(values = color_palette)  # Use the custom color palette


final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/ASE_PTC_examples_across_GTEx_tissues.png")
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

# 2.3.2) TCGA cancers
TCGA_ASE_PTCs_df_filt <- TCGA_ASE_PTCs_df[which(TCGA_ASE_PTCs_df$start_pos %in%  common_pos),]

# Order tissues by median ASE
TCGA_ASE_PTCs_df_filt <- TCGA_ASE_PTCs_df_filt %>%
  group_by(TCGA_cancer) %>%
  mutate(median_ASE = median(ASE_percentage, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(median_ASE)) %>%
  mutate(TCGA_cancer = factor(TCGA_cancer, levels = unique(TCGA_cancer)))

# Create a palette with 54 distinct colors
num_tissues <- length(unique(TCGA_ASE_PTCs_df_filt$TCGA_cancer))
color_palette <- viridis::viridis(num_tissues, option = "plasma")

# Create the plot
plot <- ggplot(TCGA_ASE_PTCs_df_filt, aes(x = TCGA_cancer, y = ASE_percentage, color = TCGA_cancer)) +
    geom_boxplot(outlier.shape = NA) +  # Avoid overlapping outliers with points
    geom_jitter(aes(alpha = 0.5), width = 0.2, size = 2) +  # Add points with alpha
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +  # Add 50% line
    theme_classic() +
    labs(
    x = "",
    y = "ASE (%)",
    title = paste0("ASE for PTC variant in position ", pos)
    ) +
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100)) +  # Set Y limits to 100
#   scale_color_brewer(palette = "Set3")  # Use a color palette for 14 tissues  
    scale_color_manual(values = color_palette)  # Use the custom color palette


final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/ASE_PTC_examples_across_TCGA_cancers.png")
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 



