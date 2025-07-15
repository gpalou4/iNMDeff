# P-values
library("psych") 
library("ggpmisc")
library("ggplot2")
library("cowplot")
library("ggpubr")
# Correlation matrix
#library("ggcorrplot")
library("corrplot")
library("dplyr")
library("viridis")
library(RColorBrewer)
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("ggpmisc")
library(readxl)
library(dplyr)
library(ggplot2)
library(coin)
library(tidyr)

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
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global_2_shared","endogenous_NMD_global_2_shared_randomized")] <- c("endogenous_NMD_Consensus","endogenous_NMD_Consensus_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global","endogenous_NMD_global_randomized")] <- c("endogenous_NMD_all","endogenous_NMD_all_randomized")
  # ASE
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.2","ASE_stopgain_0.2_randomized")] <- c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_triggering_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.01","ASE_stopgain_0.01_randomized")] <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_triggering_0.01_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.2","ASE_stopgain_NMD_evading_0.2_randomized")] <- c("ASE_PTC_NMD_evading_0.2","ASE_PTC_NMD_evading_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.01","ASE_stopgain_NMD_evading_0.01_randomized")] <- c("ASE_PTC_NMD_evading_0.01","ASE_PTC_NMD_evading_0.01_randomized")

  # Scale NMD genesets for the three methods
  # Change the sign (coefficients are reversed), so higher values means high NMDeff
  sample_NMDeff[,all_NMD_genesets] <- -sample_NMDeff[,all_NMD_genesets]
  # Scale
  if (isTRUE(scale)) {
      sample_NMDeff[,all_NMD_genesets] <- scale(sample_NMDeff[,all_NMD_genesets])
  }
  # Filter samples with low PTC number in ASE
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.2 < 3),c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","NMDeff_mean")] <- NA
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.01 < 3),c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01")] <- NA

  return(sample_NMDeff)
}

############################
######### TCGA #############
############################

# 1) Data
NMD_genesets <- c("endogenous_NMD_global_2_shared","ASE_stopgain_0.01","ASE_stopgain_0.2","PTCs_stopgain_NMD_triggering")

# 1.1) sample NMD efficiencies TCGA
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")

# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

# 1.2) TCGA metadata
TCGA_metadata <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_clinical_all.tsv", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TCGA_metadata <- TCGA_metadata[,colnames(TCGA_metadata) %in% c("submitter_id", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")]
colnames(TCGA_metadata) <- c("sample", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")

# 1.3) TCGA SKCM immunotherapy data

# Merge and add survival variables
sample_NMD_efficiencies_TCGA[,c("vital_status","days_to_death","days_to_last_follow_up")] <- NULL
TCGA_NMDeff <- merge(sample_NMD_efficiencies_TCGA, TCGA_metadata, all.x = TRUE)
TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$days_to_last_follow_up <- TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$days_to_death
TCGA_NMDeff$status <- NA
TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$status <- 1
TCGA_NMDeff[TCGA_NMDeff$vital_status == "alive",]$status <- 0
TCGA_NMDeff$days_to_last_follow_up <- as.numeric(TCGA_NMDeff$days_to_last_follow_up)/365.25
TCGA_NMDeff$age_at_diagnosis <- as.numeric(TCGA_NMDeff$age_at_diagnosis)/365.25

#Drug therapy
TCGA_cancer_metadata_drug_therapy_all <- c()
cancer_types_original <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))
for (cancer in cancer_types_original) {
  print(cancer)
  cancer <- tolower(cancer)
  # Drug therapy metadata
  if (cancer != "laml") {
    TCGA_cancer_metadata_drug_therapy <- read.csv(paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_gdc_clinical/nationwidechildrens.org_clinical_drug_",cancer,".txt"), 
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cols_keep <- c("bcr_patient_barcode","pharmaceutical_therapy_drug_name","pharmaceutical_therapy_type")
    if ("treatment_best_response" %in% colnames(TCGA_cancer_metadata_drug_therapy)) {
      cols_keep <- c(cols_keep,"treatment_best_response")
      TCGA_cancer_metadata_drug_therapy <- TCGA_cancer_metadata_drug_therapy[-1,cols_keep]
    } else {
      TCGA_cancer_metadata_drug_therapy <- TCGA_cancer_metadata_drug_therapy[-1,cols_keep]
      TCGA_cancer_metadata_drug_therapy$treatment_best_response <- NA
    }
    colnames(TCGA_cancer_metadata_drug_therapy)[colnames(TCGA_cancer_metadata_drug_therapy) %in% "bcr_patient_barcode"] <- "sample"
    if (length(TCGA_cancer_metadata_drug_therapy_all) == 0) {
      TCGA_cancer_metadata_drug_therapy_all <- TCGA_cancer_metadata_drug_therapy
    } else {
      TCGA_cancer_metadata_drug_therapy_all <- rbind(TCGA_cancer_metadata_drug_therapy_all,TCGA_cancer_metadata_drug_therapy)
    }
  }
  if (cancer == "coad") {next} else if (cancer == "read") { cancer <- "coadread"}
}

# TCGA_df <- TCGA_NMDeff[TCGA_NMDeff$cancer_type == "SKCM",]
TCGA_df <- TCGA_NMDeff

# 1.4) TCGA cell type fractions
TCIA <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCIA/TCIA-CellTypeFractionsData.tsv",
          header = FALSE, sep = "\t")
TCIA <- TCIA[!duplicated(TCIA),]
colnames(TCIA) <- c("TCGA_sample","TCGA_cancer","cell_type","quanTIseq_lsei_TIL10","Cibersort_LM22")
TCIA <- TCIA[which(TCIA$cell_type == "CD8 T cells"),]

TCGA_TCIA <- merge(TCGA_df, TCIA, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)
# TCGA_TCIA <- TCGA_TCIA[,c("cancer_type_strat","treatment_response","cell_type","endogenous_NMD_Consensus","quanTIseq_lsei_TIL10","Cibersort_LM22")]
# cibertsort_cell_type <- c("CD8 T cells", "Regulatory T cells","Monocyte","Macrophage M1","Neutrophil","Macrophage M2")
# quanTIseq_lsei_TIL10 <- c("Uncharacterized cells", "CD8 T cells", "Natural killer cells", "Regulatory T cells", 
#                         "Monocyte", "Macrophage M1", "Neutrophil", "Dendritic cells", "Macrophage M2", "B cells", "CD4 T cells")

# 2) Plots with all samples (no split by treatment)

TCGA_TCIA <- TCGA_TCIA %>%
  mutate(cancer_type_new = case_when(
    cancer_type %in% c("KICH","KIRC","KIRP") ~ "Pan-kidney",
    cancer_type %in% "BRCA" ~ "BRCA",
    cancer_type %in% "SKCM" ~ "SKCM",
    TRUE ~ "Other" 
  ))

for (percentile in c(0.1,0.2,0.3,0.4,0.5)) {

  print(percentile)
  TCGA_TCIA <- TCGA_TCIA %>%
    group_by(cancer_type) %>%
    mutate(
      NMD_type = case_when(
        endogenous_NMD_Consensus >= quantile(endogenous_NMD_Consensus, 1-percentile, na.rm = TRUE) ~ "High",
        endogenous_NMD_Consensus <= quantile(endogenous_NMD_Consensus, percentile, na.rm = TRUE) ~ "Low",
        TRUE ~ NA_character_ # Assign NA to middle 40%
      ),
      TIB_type = case_when(
        TIB >= quantile(TIB, 0.5, na.rm = TRUE) ~ "TIB-high",
        TIB <= quantile(TIB, 0.5, na.rm = TRUE) ~ "TIB-low",
        TRUE ~ NA_character_ # Assign NA to middle 40%
      )
    )

  # head(TCGA_TCIA)
  # table(TCGA_TCIA$cancer_type,TCGA_TCIA$NMD_type)

  # 2) Boxplots
  # Filter the data to ensure no missing values in proportions
  filtered_data <- TCGA_TCIA %>%
    filter(!is.na(Cibersort_LM22), !is.na(quanTIseq_lsei_TIL10), !is.na(NMD_type)) #%>%
    # filter(Cibersort_LM22 != 0 | quanTIseq_lsei_TIL10 != 0) #%>%
  #   filter(cancer_type_strat %in% c("SKCM","UCEC","UCEC_MSI","UCEC_POLE","STAD","PAAD","READ","LUAD","LUSC","OV","LGG"))

  # Boxplot for Cibersort_LM22
  plot <- ggplot(filtered_data, aes(x = NMD_type, y = Cibersort_LM22*100, fill = NMD_type)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
    geom_point(
      aes(color = NMD_type), # Map points to the NMD_type group
      position = position_dodge(0.8), # Ensure points align with the boxplots
      size = 2, # Adjust point size
      alpha = 0.7 # Adjust transparency of points
    ) +
    facet_wrap(.~cancer_type_new, scale = "free_y") +
    labs(
      title = "Cibersort_LM22",
      x = "",
      y = "CD8 T %"
    ) +
    theme_classic() +
  #   coord_cartesian(ylim=c(0,0.2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      # aes(group = interaction(TIB_type, NMD_type)),
      method = "wilcox.test",
      # method.args = list(alternative = "less"),
      label = "p.format"
    )
  #   scale_fill_manual(values = c("High" = "blue", "Low" = "red")) +
  #   scale_color_manual(values = c("High" = "blue", "Low" = "red")) # Ensure point colors match fill

  ggsave(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/CD8T_cibersort_prop_perc_",percentile,".png"), 
          plot, width = 2500, height = 2500, units = "px")

  # Boxplot for Cibersort_LM22
  plot <- ggplot(filtered_data, aes(x = NMD_type, y = quanTIseq_lsei_TIL10*100, fill = NMD_type)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
    geom_point(
      aes(color = NMD_type), # Map points to the NMD_type group
      position = position_dodge(0.8), # Ensure points align with the boxplots
      size = 2, # Adjust point size
      alpha = 0.7 # Adjust transparency of points
    ) +
    facet_wrap(.~cancer_type_new, scale = "free_y") +
    labs(
      title = "quanTIseq_lsei_TIL10",
      x = "",
      y = "CD8 T %"
    ) +
    theme_classic() +
  #   coord_cartesian(ylim=c(0,0.2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      # aes(group = interaction(TIB_type, NMD_type)),
      method = "wilcox.test",
      # method.args = list(alternative = "less"),
      label = "p.format"
    )
  #   scale_fill_manual(values = c("High" = "blue", "Low" = "red")) +
  #   scale_color_manual(values = c("High" = "blue", "Low" = "red")) # Ensure point colors match fill

  ggsave(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/CD8T_prop_quantiseq_perc_",percentile,".png"), 
          plot, width = 2500, height = 2500, units = "px")

}

# 3) Plots with treated samples

TCGA_TCIA_treated <- merge(TCGA_TCIA,TCGA_cancer_metadata_drug_therapy_all, by.x = "sample", by.y = "sample")

TCGA_therapy_counts <- TCGA_TCIA_treated %>%
      dplyr::group_by(sample) %>%
      dplyr::summarize(immunotherapy_counts = sum(pharmaceutical_therapy_type == "Immunotherapy"),
              chemotherapy_counts = sum(pharmaceutical_therapy_type == "Chemotherapy"),
              other_therapy_counts = sum(!pharmaceutical_therapy_type %in% c("Immunotherapy","Chemotherapy")))

# Immunotherapy treated patients

# samples_to_keep <- as.character(TCGA_therapy_counts[TCGA_therapy_counts$immunotherapy_counts >= 1,"sample"]$sample)
# TCGA_df_filt <- TCGA_df %>%
#     filter(sample %in% samples_to_keep) %>%
#     filter(pharmaceutical_therapy_type %in% c("Immunotherapy")) #%>%# grep("ivolumab",TCGA_df$pharmaceutical_therapy_drug_name)

# Chemotherapy treated patients
samples_to_keep <- as.character(TCGA_therapy_counts[TCGA_therapy_counts$chemotherapy_counts >= 1 & TCGA_therapy_counts$immunotherapy_counts == 0,"sample"]$sample)
# TCGA_df_filt <- TCGA_df %>%
#     dplyr::filter(sample %in% samples_to_keep) %>%
#     dplyr::filter(pharmaceutical_therapy_type %in% c("Chemotherapy")) %>%
#     dplyr::filter(!grepl("[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab|[Cc]emiplimab|[Aa]tezolizumab|[Aa]velumab|[Dd]urvalumab|[Dd]ostarlimab|[Rr]etifanlimab",pharmaceutical_therapy_drug_name))


# ipilimumab == anti CTLA-4 --> metastatic melanoma
# pembrolizumab and nivolumab (anti PD-1) --> melanoma and KIRC, LUSC

# table(TCGA_TCIA_treated[TCGA_TCIA_treated$cancer_type == "SKCM","pharmaceutical_therapy_type"])
# df <- TCGA_TCIA_treated[TCGA_TCIA_treated$cancer_type == "SKCM" & TCGA_TCIA_treated$pharmaceutical_therapy_type == "Targeted Molecular therapy",]
# sort(unique(df$pharmaceutical_therapy_drug_name))
# filter <- grep("[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab|[Cc]emiplimab|[Aa]tezolizumab|[Aa]velumab|[Dd]urvalumab|[Dd]ostarlimab|[Rr]etifanlimab",df$pharmaceutical_therapy_drug_name)
# sort(unique(df[-filter,"pharmaceutical_therapy_drug_name"]))

# add to immuno

TCGA_TCIA_treated <- TCGA_TCIA_treated %>%
  mutate(
    therapy_category = case_when(
      # Immunotherapy condition
      pharmaceutical_therapy_type %in% c("Immunotherapy") ~ "immunotherapy",
      # Chemotherapy condition only for samples in samples_to_keep
      sample %in% samples_to_keep &
        pharmaceutical_therapy_type %in% c("Chemotherapy") &
        !grepl(
          "[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab|[Cc]emiplimab|[Aa]tezolizumab|[Aa]velumab|[Dd]urvalumab|[Dd]ostarlimab|[Rr]etifanlimab",
          pharmaceutical_therapy_drug_name
        ) ~ "chemotherapy",
      # Default case
      TRUE ~ NA_character_
    )
  )
TCGA_TCIA_treated <- TCGA_TCIA_treated %>%
  mutate(
    therapy_category = case_when(
      # Immunotherapy condition
      pharmaceutical_therapy_type %in% c("Immunotherapy","Chemotherapy") ~ "chemo_immuno",
      # Default case
      TRUE ~ NA_character_
    )
  )

table(TCGA_TCIA_treated$therapy_category,TCGA_TCIA_treated$pharmaceutical_therapy_type)

# Each sample can be treated with >1 immunotherapy drug
# Let's go one by one and choose the best treatment

# Define treatment hierarchy
response_hierarchy <- c("Complete Response", "Partial Response", "Stable Disease", 
                        "[Not Applicable]", "[Not Available]", "[Unknown]")

# Create a new column with response ranks
TCGA_TCIA_treated <- TCGA_TCIA_treated %>%
  mutate(
    response_rank = match(treatment_best_response, response_hierarchy)
  )

# Select the best treatment row per sample for each treatment class
TCGA_TCIA_best_treated <- TCGA_TCIA_treated %>%
  group_by(sample, therapy_category) %>% # Group by sample and therapy type
  arrange(response_rank) %>%                       # Sort by response rank (lower is better)
  slice_head(n = 1) %>%                            # Take the top row for each group
  ungroup()                                        # Ungroup for further processing

# ipilimumab == anti CTLA-4 --> metastatic melanoma
# pembrolizumab and nivolumab (anti PD-1) --> melanoma and KIRC, LUSC

table(TCGA_TCIA_best_treated$treatment_best_response)
TCGA_TCIA_best_treated$treatment_response <- TCGA_TCIA_best_treated$treatment_best_response 
TCGA_TCIA_best_treated$treatment_response <- ifelse(TCGA_TCIA_best_treated$treatment_response == "Clinical Progressive Disease" | TCGA_TCIA_best_treated$treatment_response == "Stable Disease","Non-Responders",TCGA_TCIA_best_treated$treatment_response)
TCGA_TCIA_best_treated$treatment_response <- ifelse(TCGA_TCIA_best_treated$treatment_response == "Partial Response" | TCGA_TCIA_best_treated$treatment_response == "Complete Response","Responders",TCGA_TCIA_best_treated$treatment_response)
TCGA_TCIA_best_treated$treatment_response <- ifelse(TCGA_TCIA_best_treated$treatment_response == "[Not Applicable]" | TCGA_TCIA_best_treated$treatment_response == "[Not Available]" | TCGA_TCIA_best_treated$treatment_response == "[Unknown]" |TCGA_TCIA_best_treated$treatment_response == "[Discrepancy]",
                                    NA,TCGA_TCIA_best_treated$treatment_response)
TCGA_TCIA_best_treated[is.na(TCGA_TCIA_best_treated$therapy_category),"therapy_category"] <- "other"
table(TCGA_TCIA_best_treated$treatment_response,TCGA_TCIA_best_treated$therapy_category)

# Merge non-treated samples
tmp <- TCGA_TCIA[!TCGA_TCIA$sample %in% unique(TCGA_TCIA_best_treated$sample),] 
tmp[,c("pharmaceutical_therapy_drug_name","pharmaceutical_therapy_type",
  "treatment_best_response","therapy_category","response_rank","treatment_response")] <- NA
tmp$treatment_response <- "Non-treated"
tmp$therapy_category <- "Non-treated"
TCGA_TCIA_best_treated_final <- rbind(TCGA_TCIA_best_treated,tmp)
table(TCGA_TCIA_best_treated_final$treatment_response,TCGA_TCIA_best_treated_final$therapy_category)

TCGA_TCIA_best_treated_final$treatment_response <- factor(TCGA_TCIA_best_treated_final$treatment_response,
                                                levels = c("Non-treated","Non-Responders","Responders"))
# TCGA_TCIA_best_treated_final$therapy_category <- factor(TCGA_TCIA_best_treated_final$therapy_category,
#                                                 levels = c("chemotherapy","immunotherapy","other","Non-treated"))
TCGA_TCIA_best_treated_final$therapy_category <- factor(TCGA_TCIA_best_treated_final$therapy_category,
                                                levels = c("chemo_immuno","other","Non-treated"))

TCGA_TCIA_best_treated_final <- TCGA_TCIA_best_treated_final %>%
  filter(therapy_category %in% c("chemo_immuno"))

TCGA_TCIA_best_treated_final <- TCGA_TCIA_best_treated_final %>%
  mutate(cancer_type_new = case_when(
    cancer_type %in% c("KIRP","KIRC") ~ "Pan-kidney",
    # cancer_type %in% "KIRC" ~ "KIRC",
    # cancer_type %in% "KICH" ~ "KICH",
    # cancer_type %in% "KIRP" ~ "KIRP",
    cancer_type %in% "BRCA" ~ "BRCA",
    cancer_type %in% "SKCM" ~ "SKCM",
    TRUE ~ "Other" 
  ))

# Define percentiles
percentiles <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Create a combined dataframe with NMD_type and TIB_type for each percentile
combined_data <- bind_rows(
  lapply(percentiles, function(percentile) {
    TCGA_TCIA_best_treated_final %>%
      group_by(cancer_type) %>%
      mutate(
        NMD_type = case_when(
          endogenous_NMD_Consensus >= quantile(endogenous_NMD_Consensus, 1 - percentile, na.rm = TRUE) ~ "High",
          endogenous_NMD_Consensus <= quantile(endogenous_NMD_Consensus, percentile, na.rm = TRUE) ~ "Low",
          TRUE ~ NA_character_
        ),
        TIB_type = case_when(
          TIB >= quantile(TMB, 0.5, na.rm = TRUE) ~ "TIB-high",
          TIB <= quantile(TMB, 0.5, na.rm = TRUE) ~ "TIB-low",
          TRUE ~ NA_character_
        ),
        percentile = percentile # Add the percentile as a column
      )
  })
)

# df <- TCGA_TCIA_best_treated_final[TCGA_TCIA_best_treated_final$cancer_type=="BRCA",]
# summary(df$TIB)

# Filter the data to ensure no missing values
filtered_data <- combined_data %>%
  filter(
    !is.na(Cibersort_LM22), 
    !is.na(quanTIseq_lsei_TIL10), 
    !is.na(NMD_type)
    # Cibersort_LM22 != 0 | quanTIseq_lsei_TIL10 != 0
  ) %>%
  filter(percentile == 0.3) %>%
  mutate(NMD_type = factor(NMD_type, levels = c("High", "Low"))) %>%
  filter(cancer_type_new != "Other")

library(rstatix)

# Calculate Wilcoxon test for each facet
facet_pvalues <- filtered_data %>%
  group_by(cancer_type_new) %>% # Group by facets
  rstatix::wilcox_test(
    Cibersort_LM22 ~ NMD_type,
    alternative = "less" # Always test High < Low
  ) %>%
  mutate(p.format = format(p, digits = 2, nsmall = 2)) %>% # Add formatted p-value column
  mutate(NMD_type = "High")

# df <- filtered_data[filtered_data$cancer_type_new %in% c("Other"),]
# df$NMD_type <- factor(df$NMD_type, levels = c("Low","High"))
# df$NMD_type <- factor(df$NMD_type, levels = c("High","Low"))
# test_result <- oneway_test(quanTIseq_lsei_TIL10 ~ NMD_type, data = df, distribution = approximate(B = 10000))
# test_result

# Plot Cibersort_LM22 with percentiles on the x-axis
plot_cibersort <- ggplot(filtered_data, aes(x = as.factor(NMD_type), y = Cibersort_LM22 * 100, fill = NMD_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
  geom_point(
    aes(color = NMD_type), 
    position = position_dodge(0.8), 
    size = 2, 
    alpha = 0.7
  ) +
  facet_wrap(.~cancer_type_new, scale = "free_y") +
  labs(
    title = "Cibersort_LM22",
    x = "",
    y = "CD8 T %"
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 15)) +
  theme(axis.text.x = element_blank()) +
  # stat_compare_means(
  #   aes(group = NMD_type,
  #     # label = sprintf("p=%5.4f", as.numeric(..p.format..))
  #     ),
  #   method = "wilcox.test",
  #   label = "p.format",
  #   label.y = 18,
  #   method.args = list(alternative = "less")
  # )
  stat_pvalue_manual(
    facet_pvalues, 
    label = "p.format", # Use formatted p-values
    y.position = 12.5
  )

# Calculate Wilcoxon test for each facet
facet_pvalues <- filtered_data %>%
  group_by(cancer_type_new) %>% # Group by facets
  rstatix::wilcox_test(
    quanTIseq_lsei_TIL10 ~ NMD_type,
    alternative = "less" # Always test High < Low
  ) %>%
  mutate(p.format = format(p, digits = 2, nsmall = 2)) %>% # Add formatted p-value column
  mutate(NMD_type = "High")

# Plot quanTIseq_lsei_TIL10 with percentiles on the x-axis
plot_quantiseq <- ggplot(filtered_data, aes(x = as.factor(NMD_type), y = quanTIseq_lsei_TIL10 * 100, fill = NMD_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
  geom_point(
    aes(color = NMD_type), 
    position = position_dodge(0.8), 
    size = 2, 
    alpha = 0.7
  ) +
  facet_wrap(.~cancer_type_new, scale = "free_y") +
  labs(
    title = "quanTIseq_lsei_TIL10",
    x = "",
    y = "CD8 T %"
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 10)) +
  theme(axis.text.x = element_blank()) +
  stat_pvalue_manual(
    facet_pvalues, 
    label = "p.format", # Use formatted p-values
    y.position = 7
  )

# Save plots
ggsave("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/CD8T_cibersort_percentile_plot.png", 
       plot_cibersort, width = 2000, height = 1000, units = "px")

ggsave("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/CD8T_quantiseq_percentile_plot.png", 
       plot_quantiseq, width = 2000, height = 1000, units = "px")

saveRDS(filtered_data, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/panel_B.RData")

############################
######### HARTWIG ##########
############################

# 1.1) CD8T proportions from UCD deconvolution
infiltration_Hartwig <- read.csv(file = "/g/strcombio/fsupek_cancer1/jlanillos/UCDeconvolve/Hartwig/ucd_HMF_AdjTPM/ucd_HMF_AdjTPM_preds.csv",
        header = TRUE, sep = "\t")
infiltration_Hartwig <- infiltration_Hartwig[,grep("X|cd8",colnames(infiltration_Hartwig))]
infiltration_Hartwig <- infiltration_Hartwig[,c("X","cd8.positive..alpha.beta.cytotoxic.t.cell",
        "cd8.positive..alpha.beta.t.cell","activated.cd8.positive..alpha.beta.t.cell",
        "cd8.positive..alpha.beta.memory.t.cell","cd8.positive..alpha.beta.cytokine.secreting.effector.t.cell"), drop = FALSE]
colnames(infiltration_Hartwig) <- c("sample_ID","CD8T_cytotoxic","CD8T_t_cell","CD8T_positive_t_cell","CD8T_memory","CD8T_cytokine")
# 1.1) CD8T proportions from 
# infiltration_Hartwig <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/Hartwig/RNA_seq_immune_infiltration_Martinez_Jimenez_et_al_2023.xlsx", 
#       sheet = "Immune infiltration (RNA deconv")
# infiltration_Hartwig <- data.frame(infiltration_Hartwig)
# # Convert each column to numeric safely
# infiltration_Hartwig[, c(5:10)] <- lapply(infiltration_Hartwig[, c(5:10)], function(col) {
#   as.numeric(as.character(col))  # Convert factors or characters to numeric
# })
# 1.2) proxy iNMDeff estimates
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Hartwig_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# table(unique(Immuno_NMDeff_all$hmfSampleId) %in% infiltration_Hartwig$sample_id)
# table(infiltration_Hartwig$sample_id %in% Immuno_NMDeff_all$hmfSampleId)
# Missing samples are from PCAWG
infiltration_Hartwig_iNMDeff <- merge(Immuno_NMDeff_all,infiltration_Hartwig, by.x = "Row.names", by.y = "sample_ID")
# infiltration_Hartwig_iNMDeff <- merge(Immuno_NMDeff_all,infiltration_Hartwig, by.x = "hmfSampleId", by.y = "sample_id")

# Metadata 2
# Hartwig_metadata_2 <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/Hartwig/Hartwig_metadata_2.xlsx")
# Hartwig_metadata_2 <- data.frame(Hartwig_metadata_2)
table(infiltration_Hartwig_iNMDeff$consolidatedTreatmentType)

# 2) Analysis in treated samples
percentiles <- c(0.1, 0.2, 0.3, 0.4, 0.5)

combined_data <- bind_rows(
  lapply(percentiles, function(percentile) {
    infiltration_Hartwig_iNMDeff %>%
      # filter(treatmentGiven == "Yes") %>%
      filter(dataset == "GTEx") %>%
      filter(NMD_method == "endogenous_NMD_Consensus") %>%
      filter(consolidatedTreatmentType %in% c("Chemotherapy","Immunotherapy")) %>%
      group_by(primaryTumorLocation) %>%
      mutate(
        NMD_type = case_when(
          NMDeff >= quantile(NMDeff, 1 - percentile, na.rm = TRUE) ~ "High",
          NMDeff <= quantile(NMDeff, percentile, na.rm = TRUE) ~ "Low",
          TRUE ~ NA_character_
        ),
        percentile = percentile # Add the percentile as a column
      )
  })
)

filtered_data <- combined_data %>%
  filter(
    !is.na(CD8T_t_cell), 
    # !is.na(CD8T_parent), 
    # !is.na(cd8_davoli),
    !is.na(NMD_type)
    # Cibersort_LM22 != 0 | quanTIseq_lsei_TIL10 != 0
  ) %>%
  filter(percentile == 0.3) %>%
  mutate(NMD_type = factor(NMD_type, levels = c("High", "Low"))) %>%
  filter(primaryTumorLocation %in% c("Skin","Breast","Lung","Kidney"))
  # filter(!primaryTumorLocation %in% c("Adrenal gland","Appendix","Fallopian tube","Gallbladder",
  # "Gastroesophageal","Hepatobiliary system","Penis","Stomach","Small intestine","Testis","Thymus",
  # "Thyroid gland","Unknown","Vagina","Vulva","Anus","Bile duct"))

immune_cells <- c("CD8T_cytotoxic","CD8T_t_cell","CD8T_positive_t_cell","CD8T_memory","CD8T_cytokine")

# Reshape the dataframe from wide to long format
filtered_data <- filtered_data %>%
  pivot_longer(
    cols = all_of(immune_cells),  # Columns to stack
    names_to = "immune_cell_type",  # New column for immune cell types
    values_to = "value"  # New column for values
  )

head(data.frame(filtered_data))

# for (imm_cell in c("cd4_davoli","cd8_davoli","nk","infiltration_davoli","ifn.gamma","t_cell_grasso")) {

# print(imm_cell)
# Plot quanTIseq_lsei_TIL10 with percentiles on the x-axis
plot <- ggplot(filtered_data, aes(x = factor(NMD_type), y = value*100, fill = NMD_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
  geom_point(
    aes(color = NMD_type), 
    position = position_dodge(0.8), 
    size = 2, 
    alpha = 0.7
  ) +
  facet_grid(immune_cell_type ~ primaryTumorLocation, scale = "free_y") +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_classic() #+
  # coord_cartesian(ylim = c(0, 0.05))
  # theme(axis.text.x = element_blank()) +
  # stat_compare_means(
  #   # aes(group = interaction(TIB_type, NMD_type)),
  #   # label.y = 0.06,
  #   method = "wilcox.test",
  #   # method.args = list(alternative = "less"),
  #   label = "p.format"
  # )
  # stat_pvalue_manual(
  #   facet_pvalues, 
  #   label = "p.format", # Use formatted p-values
  #   y.position = 7
  # )

# Save plots
# ggsave(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/Hartwig_Martinez_Jimenez_et_al_2023_percentile_plot.png"), 
ggsave(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/Hartwig_UCD_CD8T.png"), 
      plot, width = 3000, height = 3000, units = "px")

# }



# Define specific y-axis ranges for each cell type
cell_type_ranges <- list(
  "CD8T_cytotoxic" = c(0, 0.025),
  "CD8T_t_cell" = c(0, 20),
  "CD8T_positive_t_cell" = c(0, 0.02),
  "CD8T_memory" = c(0, 0.01),
  "CD8T_cytokine" = c(0, 20)
)

# Create a named vector for renaming cell types
cell_type_labels <- c(
  "CD8T_cytotoxic" = "CD8T Cytotoxic",
  "CD8T_t_cell" = "CD8T T-cell",
  "CD8T_positive_t_cell" = "CD8T Positive T-cell",
  "CD8T_memory" = "CD8T Memory",
  "CD8T_cytokine" = "CD8T Cytokine"
)

# Initialize an empty list to store the plots
plots <- list()
is_first_plot <- TRUE

# Loop over each immune cell type to create individual plots
for (cell_type in names(cell_type_ranges)) {
  # Subset the data for the current immune cell type
  subset_data <- filtered_data %>% 
    filter(immune_cell_type == cell_type)
  
  # Get the y-axis range for the current cell type
  y_range <- cell_type_ranges[[cell_type]]

  legend_position <- if (is_first_plot) "top" else "none"

  # Create the plot
  plot <- ggplot(subset_data, aes(x = factor(NMD_type), y = value * 100, fill = NMD_type)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
    geom_point(
      aes(color = NMD_type), 
      position = position_dodge(0.8), 
      size = 2, 
      alpha = 0.7
    ) +
    facet_grid(immune_cell_type~primaryTumorLocation, scale = "free_y",
            labeller = labeller(immune_cell_type = cell_type_labels)) +
    labs(
      title = "",
      x = "",
      y = "cell type (%)",
      fill = "ETG iNMDeff"
    ) +
    theme_classic() +
    coord_cartesian(ylim = y_range) +
    stat_compare_means(
      # aes(group = interaction(TIB_type, NMD_type)),
      label.y = y_range[2]-(y_range[2]*0.1),
      label.x = 1.5,
      size = 3.5,
      method = "wilcox.test",
      # method.args = list(alternative = "less"),
      label = "p.signif"
      # label = function(x) sprintf("%.3g", x$p)
    ) + guides(color = "none") + 
    theme(strip.text = element_text(size = 8),
            legend.position = legend_position
            )
  # After the first plot, set is_first_plot to FALSE
  is_first_plot <- FALSE
  
  # Save the plot in the list
  plots[[cell_type]] <- plot
}

figure <- plot_grid(plots[[1]], plots[[2]],plots[[3]],plots[[4]],plots[[5]], nrow = 5,
                       labels = c(""), label_size = 20, align = "v", rel_heights = c(0.23,0.2,0.2,0.2,0.2))

ggsave(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/immuno/Hartwig_UCD_CD8T.png"), 
      figure, width = 1500, height = 3500, units = "px")

# Save


saveRDS(filtered_data, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig28/panel_A.RData")
