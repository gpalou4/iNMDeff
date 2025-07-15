modify_NMDeff_dataframe <- function(sample_NMDeff, dataset, scale = FALSE) {
  # Convert some columns to factors
  if (dataset == "TCGA") {
    factor_cols <- c("cancer_type","cancer_type_MSI","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                    "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
    # Remove "TCGA" from the cancer type
    sample_NMDeff$cancer_type <- gsub("TCGA-","",sample_NMDeff$cancer_type)
    sample_NMDeff$cancer_type_strat <- gsub("TCGA-","",sample_NMDeff$cancer_type_strat)
    sample_NMDeff$cancer_type_MSI <- gsub("TCGA-","",sample_NMDeff$cancer_type_MSI)
  } else if (dataset == "GTEx") {
    factor_cols <- c("tissue","sample")
  }
  sample_NMDeff[factor_cols] <- lapply(sample_NMDeff[factor_cols], factor) 
  # Rename NMD genesets for the 3 methods
  all_NMD_genesets <- c(endogenous_NMD_genesets,ASE_NMD_genesets,"NMDeff_mean")
  # Endogenous
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global_2_shared","endogenous_NMD_global_2_shared_randomized")] <- c("NMD Consensus","NMD Consensus randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global","endogenous_NMD_global_randomized")] <- c("NMD All","NMD All randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Karousis","endogenous_NMD_Karousis_randomized")] <- c("NMD Karousis","NMD Karousis randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Colombo","endogenous_NMD_Colombo_randomized")] <- c("NMD Colombo","NMD Colombo randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_non_NMD_neg_control_with_NMD_features","endogenous_non_NMD_neg_control_with_NMD_features_randomized")] <- c("RandomGenes with NMD features","RandomGenes with NMD features - randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_randomized")] <- c("RandomGenes without NMD features","RandomGenes without NMD features - randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Courtney","endogenous_NMD_Courtney_randomized")] <- c("NMD Courtney","NMD Courtney randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_ensembl","endogenous_NMD_ensembl_randomized")] <- c("NMD Ensembl","NMD Ensembl randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_SMG6","endogenous_SMG6_randomized")] <- c("NMD SMG6","NMD SMG6 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_SMG7","endogenous_SMG7_randomized")] <- c("NMD SMG7","NMD SMG7 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Tani","endogenous_NMD_Tani_randomized")] <- c("NMD Tani","NMD Tani randomized")
  # ASE
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.2","ASE_stopgain_0.2_randomized")] <- c("PTC NMD-triggering 0.2","PTC NMD-triggering 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.01","ASE_stopgain_0.01_randomized")] <- c("PTC NMD-triggering 0.01","PTC NMD-triggering 0.01 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.2","ASE_stopgain_NMD_evading_0.2_randomized")] <- c("PTC NMD-evading 0.2","PTC NMD-evading 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.01","ASE_stopgain_NMD_evading_0.01_randomized")] <- c("PTC NMD-evading 0.01","PTC NMD-evading 0.01 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_synonymous_0.2","ASE_synonymous_0.2_randomized")] <- c("Synonymous 0.2","Synonymous 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_synonymous_0.01","ASE_synonymous_0.01_randomized")] <- c("Synonymous 0.01","Synonymous 0.01 randomized")
  # Scale NMD genesets for the three methods
  # Change the sign (coefficients are reversed), so higher values means high NMDeff
  sample_NMDeff[,all_NMD_genesets] <- -sample_NMDeff[,all_NMD_genesets]
  # Scale
  if (isTRUE(scale)) {
      sample_NMDeff[,all_NMD_genesets] <- scale(sample_NMDeff[,all_NMD_genesets])
  }
  # Filter samples with low PTC number in ASE
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.2 < 3),c("PTC NMD-triggering 0.2","PTC NMD-evading 0.2","NMDeff_mean")] <- NA
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.01 < 3),c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01")] <- NA

  return(sample_NMDeff)
}

ASE_variants_filtering <- function(TCGA_sample) {

    # Obtain ASE file
    file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA-COAD/",TCGA_sample,"/PTC_transcripts/germline_PTCs_ASE_transcripts_metadata.txt")
    ASE_sample <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Filters
    LOEUF_score <- 0
    # 5) Filterings for the NB regression
    # 5.1) Coverage 5 for nonsense and 2 for FS (filter for ASE PTCs)
    filter_coverage_nonsense <- which(ASE_sample$totalCount_ASE >= 5 & ASE_sample$stopgain == "nonsense")
    filter_coverage_FS <- which(ASE_sample$totalCount_ASE >= 2 & ASE_sample$stopgain != "nonsense")
    ASE_sample_filt <- ASE_sample[c(filter_coverage_nonsense,filter_coverage_FS),]
    # 5.2) Remove ASE-PTCs from FS with no predicted PTC
    if (length(filter_coverage_FS) != 0) {
        filter_FS_no_PTC <- which(ASE_sample_filt$stopgain != "nonsense" & !is.na(ASE_sample_filt$PTC_CDS_pos)) 
    } else {
        filter_FS_no_PTC <- NULL
    }
    # 5.3) Remove ASE-PTCs with specified VAF
    filter_VAF <- which(ASE_sample_filt$VAF > 0.2)
    # 5.4) Remove transcripts with 1 exon
    filter_single_exon <- which(ASE_sample_filt$transcript_CDS_exon_num == 1 & ASE_sample_filt$splice_site_3UTR == "no")
    # 5.5) LOEUF score 1
    filter_LOEUF <- which(ASE_sample_filt$LOEUF_decile <= LOEUF_score)
    # 5.6) PSG genes
    filter_PSG <- which(ASE_sample_filt$PSG == "yes")
    # 5.6) Transcripts overlapping somatic mut/CNV
    filter_somatic_mut_CNV <- which(ASE_sample_filt$somatic_CNV_SNV == "yes")
    # 5.8) NMD-triggering vs NMD-evading
    filter_NMDtype <- which(ASE_sample_filt$X55_nt_last_exon == "NMD-evading")
    # 5.9) Non protein-coding
    filter_nonCoding <- which(ASE_sample_filt$protein_coding == "no")
    # 5.10) Overlapping germline mutations...
    # Add up filters
    all_filters <- unique(c(filter_FS_no_PTC,filter_VAF,filter_single_exon,filter_LOEUF,filter_somatic_mut_CNV,filter_NMDtype,filter_nonCoding,filter_PSG))
    ASE_sample_filt <- ASE_sample_filt[-all_filters,]

    NMD_final_df <- ASE_sample_filt[,c("stopgain","variant_type","genotype","start_pos","gene_id","refCount_ASE","altCount_ASE","totalCount_ASE","VAF")]
    colnames(NMD_final_df) <- c("stopgain","variant_type","genotype","start_pos","gene_id","refCount","altCount","totalCount","MAF")
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
endogenous_NMD_genesets <-  c("NMD Colombo","NMD Karousis","NMD Tani","NMD Courtney","NMD Ensembl",
                      "NMD All","NMD Consensus","NMD SMG6","NMD SMG7",
                      "RandomGenes without NMD features","RandomGenes with NMD features")
ASE_NMD_genesets <- c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01","Synonymous 0.01",
                      "PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2")

# PTC // ASE // Endogenosample_NMD_efficiencies_TCGAus
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# 2) Select two MSI samples with ASE-NMDeff high and ASE-NMDeff low
# TCGA-COAD
NMDeff_filt <- sample_NMD_efficiencies_TCGA %>%
            filter(cancer_type == "COAD")
summary(NMDeff_filt[,"PTC NMD-triggering 0.2"])
set.seed(234)
# set.seed(444) #This has both NMDeff high and low in both methods agreeing
filt <- which(NMDeff_filt[,"PTC NMD-triggering 0.2"] < 0 & NMDeff_filt[,"MSI_status"] == "MSI-H")
NMDeff_low <- NMDeff_filt[sample(filt)[1],]
NMDeff_low[,c("PTC NMD-triggering 0.2", "NMD Consensus")]
filt <- which(NMDeff_filt[,"PTC NMD-triggering 0.2"] > 1.5 & NMDeff_filt[,"MSI_status"] == "MSI-H")
NMDeff_high <- NMDeff_filt[sample(filt)[1],]
NMDeff_high[,c("PTC NMD-triggering 0.2", "NMD Consensus")]

# 3) Obtain their ASE data
# 3.1) Low NMDeff sample
TCGA_sample_low <- as.character(NMDeff_low$sample)
NMDeff_low_ASE_df <- ASE_variants_filtering(TCGA_sample = TCGA_sample_low)
NMDeff_low_ASE_df$TCGA_sample <- TCGA_sample_low
NMDeff_low_ASE_df$MSI_status <- "MSI-H"
# 3.2) High NMDeff sample
TCGA_sample_high <- as.character(NMDeff_high$sample)
NMDeff_high_ASE_df <- ASE_variants_filtering(TCGA_sample = TCGA_sample_high)
NMDeff_high_ASE_df$TCGA_sample <- TCGA_sample_high
NMDeff_high_ASE_df$MSI_status <- "MSI-H"
# 3.3) Merge
NMDeff_ASE_df <- rbind(NMDeff_high_ASE_df,NMDeff_low_ASE_df)
NMDeff_ASE_df$refCount_perc <- round(NMDeff_ASE_df$refCount / NMDeff_ASE_df$totalCount,2)
NMDeff_ASE_df$altCount_perc <- round(NMDeff_ASE_df$altCount / NMDeff_ASE_df$totalCount,2)
NMDeff_ASE_df <- NMDeff_ASE_df %>%
                group_by(TCGA_sample) %>%
                arrange(desc(altCount_perc)) %>%
                arrange(desc(TCGA_sample))
data.frame(NMDeff_ASE_df)

NMDeff_ASE_stacked <- stack(NMDeff_ASE_df[,c("refCount_perc","altCount_perc")])
NMDeff_ASE_stacked$MSI_status <- rep(NMDeff_ASE_df$MSI_status,2)
NMDeff_ASE_stacked$TCGA_sample <- rep(NMDeff_ASE_df$TCGA_sample,2)
NMDeff_ASE_stacked$variant_type <- rep(NMDeff_ASE_df$variant_type,2)
NMDeff_ASE_stacked$variant <- rep(1:nrow(NMDeff_ASE_df),2)
colnames(NMDeff_ASE_stacked) <- c("ASE","allele","MSI_status","TCGA_sample","variant_type","variant")
NMDeff_ASE_stacked$TCGA_sample_type <- ifelse(NMDeff_ASE_stacked$TCGA_sample == TCGA_sample_high,"High","Low")

write.table(NMDeff_ASE_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1C.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)

p <- NMDeff_ASE_stacked %>%
    ggplot(aes(y = ASE, x = variant, fill = allele, color = TCGA_sample)) + 
    geom_bar(position="stack", stat="identity", size = 1.5)  + theme_classic() + #ylim(c(0,100)) +
    labs(title = "High and Low ASE iNMDeff in Two Individuals", x = "", y = "ASE", fill = paste0("PTC"), color = "iNMDeff") + 
    scale_fill_brewer(labels = c("REF", "ALT"), palette = "Paired", direction = -1) +
    #scale_color_brewer(labels = c("REF", "ALT"), palette = "Dark2", direction = 1) +
    scale_color_manual(  
            values = c("Low" = "#2D3263", "High" = "#F07626"), 
            labels = c("High", "Low")
            ) +
    #scale_color_viridis(labels = c("Low", "High"),discrete=TRUE, option="viridis") +
    theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.title.x = element_text(color="black", size=40),
        axis.title.y = element_text(color="black", size=45),
        axis.text.y = element_text(color="black", size=40),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 30),         # Increase legend label font size
        legend.title = element_text(size = 35),
        legend.box = "vertical",
        legend.key.width = unit(1.5, "cm"),  # Adjust the key width here
        legend.key.height = unit(1.5, "cm"),
        legend.title.align = 0.5,
        legend.position = "right") +
        guides(fill = guide_legend(override.aes = list(size = 11)),
            color = guide_legend(override.aes = list(size = 11)) )

png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/schematics/ASE_NMDeffHighvsLow.png"), 
    width = 4000, height = 3500, res = 300)
print(p)
dev.off()


