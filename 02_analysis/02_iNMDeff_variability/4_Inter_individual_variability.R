# Libraries
library(dplyr)
library(tidyverse)
library(ggpubr)

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

coefficient_of_variation <- function(x) {
    sd(x) / mean(x,na.rm=TRUE)
}

RRMSE <- function(observed) {
  n <- length(observed)
  predicted <- rep(0,n)
  residuals <- predicted - observed
  rmse <- sqrt(mean(residuals^2))
  mean_obs <- mean(observed, na.rm = TRUE)
  rrmse <- rmse / mean_obs
  return(rrmse)
}

fix_df <- function(df) {
  df$intra_ind_var$type <- "Intra-individual"
  df$inter_ind_var$type <- "Inter-individual"
  df$intra_ind_var_rand$type <- "Intra-individual randomization"
  df$inter_ind_var_rand$type <- "Inter-individual randomization"
  colnames(df$intra_ind_var)[1] <- "ID"
  colnames(df$inter_ind_var)[1] <- "ID"
  colnames(df$intra_ind_var_rand)[1] <- "ID"
  colnames(df$inter_ind_var_rand)[1] <- "ID"
  NMDeff_PTC_all <- rbind( df$intra_ind_var, df$inter_ind_var, df$intra_ind_var_rand, df$inter_ind_var_rand)
  cols <- c("PTC_NMDeff_var")
  NMDeff_PTC_all_stacked <- stack(NMDeff_PTC_all[,cols])
  NMDeff_PTC_all_stacked$type <- rep(NMDeff_PTC_all$type,length(cols))
  NMDeff_PTC_all_stacked <- na.omit(NMDeff_PTC_all_stacked) 
  NMDeff_PTC_all_stacked$ind <- NULL
  colnames(NMDeff_PTC_all_stacked) <- c("variance","type")
  return(NMDeff_PTC_all_stacked)
}

PTCs_inter_intra_ind_corr <- function(samples_NMDeff_df, n) {

  # Obtain a correlation for each iteration with n = 1000
  # 2.1) Intra-individual PTC-NMDeff correlation
  # In each iteration we obtain PTC-NMDeff of two sampled PTCs for a given individual and we repeat this for each individual
  # Rows: Individuals. Column 1: Sampled PTC-1, Column 2: Sampled PTC-2. Values = PTC-NMDeff
  # samples_NMDeff_df_no_dup <- samples_NMDeff_df %>%
  #         mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
  # dup_to_keep <- unique(samples_NMDeff_df_no_dup[duplicated(samples_NMDeff_df_no_dup$PTC_ID),"PTC_ID"])
  # samples_NMDeff_df_no_dup <- samples_NMDeff_df_no_dup[-which(samples_NMDeff_df_no_dup$PTC_ID %in% dup_to_keep),]
  samples_NMDeff_df_no_dup <- samples_NMDeff_df
  # Count the number of PTCs for each individual
  individual_PTC_counts <- samples_NMDeff_df_no_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(individual) %>%
    summarise(PTCs_count = n())

  all_intra_corr <- c()
  for (i in 1:n) {
    print(i)
    Intra_NMDeff_PTC_pairs <- samples_NMDeff_df_no_dup %>%
      inner_join(individual_PTC_counts, by = "individual") %>%
      filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
      filter(PTCs_count >= 2) %>%
      group_by(individual) %>%
      summarise(PTCs = list(sample(NMD_efficiency_TPM,2)))
    Intra_NMDeff_PTC_pairs_df <- do.call(rbind,Intra_NMDeff_PTC_pairs$PTCs)
    # Correlation
    corr_res <- cor.test(Intra_NMDeff_PTC_pairs_df[,1], Intra_NMDeff_PTC_pairs_df[,2], method = "spearman", use = "pairwise.complete.obs")
    all_intra_corr <- c(all_intra_corr,as.numeric(corr_res$estimate))
  }

  # 2.2) Inter-individual PTC-NMDeff correlation
  # In each iteration we obtain PTC-NMDeff of two sampled individuals for a given shared PTC and we repeat this for each PTC
  # Rows: PTCs. Column 1: Sampled individual-1, Column 2: Sampled individual-2. Values = PTC-NMDeff
  # cols <- c('UTR5s_length', 'exons_length_prePTC', 'introns_length_prePTC', 'PTC_CDS_exon_length', 'introns_length_postPTC', 'exons_length_postPTC', 'UTR3s_length', 'PTC_CDS_pos')
  samples_NMDeff_df_dup <- samples_NMDeff_df %>%
          filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
          mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
          # mutate(PTC_ID = paste0(UTR5s_length,"_",exons_length_prePTC,"_",introns_length_prePTC,"_",PTC_CDS_exon_length,"_",introns_length_postPTC,"_",exons_length_postPTC,"_",UTR3s_length,"_",PTC_CDS_pos))
  dup_to_keep <- unique(samples_NMDeff_df_dup[duplicated(samples_NMDeff_df_dup$PTC_ID),"PTC_ID"])
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup[which(samples_NMDeff_df_dup$PTC_ID %in% dup_to_keep),]
  # Count the number of PTCs for each individual
  interindividual_PTC_counts <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(PTC_ID) %>%
    summarise(PTCs_count = n())

  all_inter_corr <- c()
  for (i in 1:n) {
    print(i)
    Inter_NMDeff_PTC_pairs <- samples_NMDeff_df_dup %>%
      inner_join(interindividual_PTC_counts, by = "PTC_ID") %>%
      filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
      filter(PTCs_count >= 2) %>%
      group_by(PTC_ID) %>%
      summarise(PTCs = list(sample(NMD_efficiency_TPM,2)))
    Inter_NMDeff_PTC_pairs_df <- do.call(rbind,Inter_NMDeff_PTC_pairs$PTCs)
    # Correlation
    corr_res <- cor.test(Inter_NMDeff_PTC_pairs_df[,1], Inter_NMDeff_PTC_pairs_df[,2], method = "spearman", use = "pairwise.complete.obs")
    all_inter_corr <- c(all_inter_corr,as.numeric(corr_res$estimate))
  }

  # 2.3) Randomization - Intra-individual PTC-NMDeff correlation
  samples_NMDeff_df_no_dup <- samples_NMDeff_df
  # Sampling of NMD efficiencies
  samples_NMDeff_df_no_dup <- samples_NMDeff_df_no_dup %>%
        group_by(tissue) %>%
        mutate(sample_size = n()) %>%
        filter(sample_size != 1) %>%
        mutate(NMD_efficiency_TPM = base::sample(NMD_efficiency_TPM))
  # samples_NMDeff_df_no_dup$NMD_efficiency_TPM <- sample(samples_NMDeff_df_no_dup$NMD_efficiency_TPM)
  # Count the number of PTCs for each individual
  individual_PTC_counts <- samples_NMDeff_df_no_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(individual) %>%
    summarise(PTCs_count = n())

  all_intra_corr_rand <- c()
  for (i in 1:n) {
    print(i)
    Intra_NMDeff_PTC_pairs <- samples_NMDeff_df_no_dup %>%
      inner_join(individual_PTC_counts, by = "individual") %>%
      filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
      filter(PTCs_count >= 2) %>%
      group_by(individual) %>%
      summarise(PTCs = list(sample(NMD_efficiency_TPM,2)))
    Intra_NMDeff_PTC_pairs_df <- do.call(rbind,Intra_NMDeff_PTC_pairs$PTCs)
    # Correlation
    corr_res <- cor.test(Intra_NMDeff_PTC_pairs_df[,1], Intra_NMDeff_PTC_pairs_df[,2], method = "spearman", use = "pairwise.complete.obs")
    all_intra_corr_rand <- c(all_intra_corr_rand,as.numeric(corr_res$estimate))
  }
  # # 2.4) Randomization - Inter-individual PTC-NMDeff correlation
  samples_NMDeff_df_dup <- samples_NMDeff_df %>%
          filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
          mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
  dup_to_keep <- unique(samples_NMDeff_df_dup[duplicated(samples_NMDeff_df_dup$PTC_ID),"PTC_ID"])
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup[which(samples_NMDeff_df_dup$PTC_ID %in% dup_to_keep),]
  # Sampling of NMD efficiencies
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup %>%
        group_by(tissue) %>%
        mutate(sample_size = n()) %>%
        filter(sample_size != 1) %>%
        mutate(NMD_efficiency_TPM = sample(NMD_efficiency_TPM))  # Count the number of PTCs for each individual
  interindividual_PTC_counts <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(PTC_ID) %>%
    summarise(PTCs_count = n())

  all_inter_corr_rand <- c()
  for (i in 1:n) {
    print(i)
    Inter_NMDeff_PTC_pairs <- samples_NMDeff_df_dup %>%
      inner_join(interindividual_PTC_counts, by = "PTC_ID") %>%
      filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
      filter(PTCs_count >= 2) %>%
      group_by(PTC_ID) %>%
      summarise(PTCs = list(sample(NMD_efficiency_TPM,2)))
    Inter_NMDeff_PTC_pairs_df <- do.call(rbind,Inter_NMDeff_PTC_pairs$PTCs)
    # Correlation
    corr_res <- cor.test(Inter_NMDeff_PTC_pairs_df[,1], Inter_NMDeff_PTC_pairs_df[,2], method = "spearman", use = "pairwise.complete.obs")
    all_inter_corr_rand <- c(all_inter_corr_rand,as.numeric(corr_res$estimate))
  }
  # Results
  # print("Intra individual correlations")
  # print(summary(all_intra_corr))
  # print("Inter individual correlations")
  # print(summary(all_inter_corr))
  # Save
  return(list(intra_ind_corr = all_intra_corr, inter_ind_corr = all_inter_corr,
              intra_ind_corr_rand = all_intra_corr_rand, inter_ind_corr_rand = all_inter_corr_rand))
}

PTCs_inter_intra_ind_variability <- function(samples_NMDeff_df) {

  # 3.1) Intra Individual PTC-NMDeff variation
  # samples_NMDeff_df_no_dup <- samples_NMDeff_df %>%
  #         mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
  # dup_to_keep <- unique(samples_NMDeff_df_no_dup[duplicated(samples_NMDeff_df_no_dup$PTC_ID),"PTC_ID"])
  # samples_NMDeff_df_no_dup <- samples_NMDeff_df_no_dup[-which(samples_NMDeff_df_no_dup$PTC_ID %in% dup_to_keep),]
  samples_NMDeff_df_no_dup <- samples_NMDeff_df
  individual_PTC_counts <- samples_NMDeff_df_no_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(individual) %>%
    summarise(PTCs_count = n())
  NMDeff_PTC_intraindividual <- samples_NMDeff_df_no_dup %>%
    inner_join(individual_PTC_counts, by = "individual") %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    filter(PTCs_count >= 3) %>%
    group_by(individual) %>%
    summarise(PTC_NMDeff_var = var(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_SD = sd(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_CV = coefficient_of_variation(NMD_efficiency_TPM),
              PTC_NMDeff_MAD= mad(NMD_efficiency_TPM),
              PTC_NMDeff_RRMSE = RRMSE(NMD_efficiency_TPM))
  # 3.2) Inter Individual PTC-NMDeff variation
  samples_NMDeff_df_dup <- samples_NMDeff_df %>%
          filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
          mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
  dup_to_keep <- unique(samples_NMDeff_df_dup[duplicated(samples_NMDeff_df_dup$PTC_ID),"PTC_ID"])
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup[which(samples_NMDeff_df_dup$PTC_ID %in% dup_to_keep),]
  # Count the number of PTCs for each individual
  interindividual_PTC_counts <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(PTC_ID) %>%
    summarise(PTCs_count = n())
  NMDeff_PTC_interindividual <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    inner_join(interindividual_PTC_counts, by = "PTC_ID") %>%
    filter(PTCs_count >= 3) %>%
    group_by(PTC_ID) %>%
    summarise(PTC_NMDeff_var = var(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_SD = sd(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_CV = coefficient_of_variation(NMD_efficiency_TPM),
              PTC_NMDeff_MAD= mad(NMD_efficiency_TPM),
              PTC_NMDeff_RRMSE = RRMSE(NMD_efficiency_TPM))
  # 3.3) Intra Individual PTC-NMDeff variation Randomization
  samples_NMDeff_df_no_dup <- samples_NMDeff_df
  # Randomize NMD efficiency
  samples_NMDeff_df_no_dup <- samples_NMDeff_df_no_dup %>%
        group_by(tissue) %>%
        mutate(sample_size = n()) %>%
        filter(sample_size != 1) %>%
        mutate(NMD_efficiency_TPM = sample(NMD_efficiency_TPM)) 
  # Count PTCs per sample
  individual_PTC_counts <- samples_NMDeff_df_no_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(individual) %>%
    summarise(PTCs_count = n())
  NMDeff_PTC_intraindividual_rand <- samples_NMDeff_df_no_dup %>%
    inner_join(individual_PTC_counts, by = "individual") %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    filter(PTCs_count >= 3) %>%
    group_by(individual) %>%
    summarise(PTC_NMDeff_var = var(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_SD = sd(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_CV = coefficient_of_variation(NMD_efficiency_TPM),
              PTC_NMDeff_MAD= mad(NMD_efficiency_TPM),
              PTC_NMDeff_RRMSE = RRMSE(NMD_efficiency_TPM))

  # 3.4) Inter Individual PTC-NMDeff variation Randomization
  samples_NMDeff_df_dup <- samples_NMDeff_df %>%
          filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
          mutate(PTC_ID = paste0(transcript_id,"_",PTC_CDS_pos,"_",start_pos))
  dup_to_keep <- unique(samples_NMDeff_df_dup[duplicated(samples_NMDeff_df_dup$PTC_ID),"PTC_ID"])
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup[which(samples_NMDeff_df_dup$PTC_ID %in% dup_to_keep),]
  # Randomize NMD efficiency
  samples_NMDeff_df_dup <- samples_NMDeff_df_dup %>%
        group_by(tissue) %>%
        mutate(sample_size = n()) %>%
        filter(sample_size != 1) %>%
        mutate(NMD_efficiency_TPM = sample(NMD_efficiency_TPM)) 
  # Count the number of PTCs for each individual
  interindividual_PTC_counts <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    group_by(PTC_ID) %>%
    summarise(PTCs_count = n())
  NMDeff_PTC_interindividual_rand <- samples_NMDeff_df_dup %>%
    filter(X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 200) %>%
    inner_join(interindividual_PTC_counts, by = "PTC_ID") %>%
    filter(PTCs_count >= 3) %>%
    group_by(PTC_ID) %>%
    summarise(PTC_NMDeff_var = var(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_SD = sd(NMD_efficiency_TPM,na.rm=TRUE),
              PTC_NMDeff_CV = coefficient_of_variation(NMD_efficiency_TPM),
              PTC_NMDeff_MAD= mad(NMD_efficiency_TPM),
              PTC_NMDeff_RRMSE = RRMSE(NMD_efficiency_TPM))
  # Save
  return(list(intra_ind_var = NMDeff_PTC_intraindividual, inter_ind_var = NMDeff_PTC_interindividual,
              intra_ind_var_rand = NMDeff_PTC_intraindividual_rand, inter_ind_var_rand = NMDeff_PTC_interindividual_rand))
}

# 1) Data
# 1.1) Samples NMD-eff
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
# TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = FALSE)

# # 1.2) PTC NMD efficiencies
# PTC
PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/germline_PTCs_all_TCGA_confident_seq.txt"
germline_PTCs_NMD_efficiencies_TCGA <- read.table(file = PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
somatic_PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/somatic_PTCs_all_TCGA_confident_seq.txt"
somatic_PTCs_NMD_efficiencies_TCGA <- read.table(file = somatic_PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
germline_PTCs_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/germline_PTCs_all_GTEx_confident_seq.txt"
germline_PTCs_NMD_efficiencies_GTEx <- read.table(file = germline_PTCs_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(germline_PTCs_NMD_efficiencies_GTEx)[colnames(germline_PTCs_NMD_efficiencies_GTEx) == "GTEx_sample"] <- "individual"
colnames(germline_PTCs_NMD_efficiencies_TCGA)[colnames(germline_PTCs_NMD_efficiencies_TCGA) == "TCGA_barcode"] <- "individual"
colnames(somatic_PTCs_NMD_efficiencies_TCGA)[colnames(somatic_PTCs_NMD_efficiencies_TCGA) == "TCGA_barcode"] <- "individual"
colnames(germline_PTCs_NMD_efficiencies_GTEx)[colnames(germline_PTCs_NMD_efficiencies_GTEx) == "GTEx_tissue"] <- "tissue"
colnames(germline_PTCs_NMD_efficiencies_TCGA)[colnames(germline_PTCs_NMD_efficiencies_TCGA) == "TCGA_cancer"] <- "tissue"
colnames(somatic_PTCs_NMD_efficiencies_TCGA)[colnames(somatic_PTCs_NMD_efficiencies_TCGA) == "TCGA_cancer"] <- "tissue"

# ASE
ASE_PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/germline_PTCs_ASE_all_TCGA_confident_seq.txt"
germline_ASE_PTCs_NMD_efficiencies_TCGA <- read.table(file = ASE_PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
somatic_ASE_PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/somatic_PTCs_ASE_all_TCGA_confident_seq.txt"
somatic_ASE_PTCs_NMD_efficiencies_TCGA <- read.table(file = somatic_ASE_PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
germline_ASE_PTCs_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/Marcell/germline_PTCs_ASE_all_GTEx_confident_seq.txt"
germline_ASE_PTCs_NMD_efficiencies_GTEx <- read.table(file = germline_ASE_PTCs_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx)[colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx) == "GTEx_sample"] <- "individual"
colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA) == "TCGA_barcode"] <- "individual"
colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA) == "TCGA_barcode"] <- "individual"
germline_ASE_PTCs_NMD_efficiencies_GTEx$NMD_efficiency_TPM <- NULL
germline_ASE_PTCs_NMD_efficiencies_TCGA$NMD_efficiency_TPM <- NULL
somatic_ASE_PTCs_NMD_efficiencies_TCGA$NMD_efficiency_TPM <- NULL
colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx)[colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx) == "ASE_NMD_efficiency_TPM"] <- "NMD_efficiency_TPM"
colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA) == "ASE_NMD_efficiency_TPM"] <- "NMD_efficiency_TPM"
colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA) == "ASE_NMD_efficiency_TPM"] <- "NMD_efficiency_TPM"
colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx)[colnames(germline_ASE_PTCs_NMD_efficiencies_GTEx) == "GTEx_tissue"] <- "tissue"
colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(germline_ASE_PTCs_NMD_efficiencies_TCGA) == "TCGA_cancer"] <- "tissue"
colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA)[colnames(somatic_ASE_PTCs_NMD_efficiencies_TCGA) == "TCGA_cancer"] <- "tissue"

# Filters
# MAF
#germline_PTCs_NMD_efficiencies_TCGA[which(germline_PTCs_NMD_efficiencies_TCGA$VAF >= 0.01),]
#germline_ASE_PTCs_NMD_efficiencies_TCGA[which(germline_PTCs_NMD_efficiencies_TCGA$VAF >= 0.01),]

# 2) Inter-individual NMDeff variability -- Correlations
######## PTCs ########
set.seed(42)
germline_TCGA_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = germline_PTCs_NMD_efficiencies_TCGA, n = 1000) 
germline_TCGA_corr_df <- data.frame(intra_ind_corr = germline_TCGA_corr_res$intra_ind_corr, 
                                    inter_ind_corr = germline_TCGA_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = germline_TCGA_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = germline_TCGA_corr_res$inter_ind_corr_rand)
germline_TCGA_corr_stack <- stack(germline_TCGA_corr_df)
germline_TCGA_corr_stack$dataset <- "germline_TCGA"
somatic_TCGA_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = somatic_PTCs_NMD_efficiencies_TCGA, n = 1000) 
somatic_TCGA_corr_df <- data.frame(intra_ind_corr = somatic_TCGA_corr_res$intra_ind_corr, 
                                    inter_ind_corr = somatic_TCGA_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = somatic_TCGA_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = somatic_TCGA_corr_res$inter_ind_corr_rand)
somatic_TCGA_corr_stack <- stack(somatic_TCGA_corr_df)
somatic_TCGA_corr_stack$dataset <- "somatic_TCGA"
germline_GTEx_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = germline_PTCs_NMD_efficiencies_GTEx, n = 1000) 
germline_GTEx_corr_df <- data.frame(intra_ind_corr = germline_GTEx_corr_res$intra_ind_corr, 
                                    inter_ind_corr = germline_GTEx_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = germline_GTEx_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = germline_GTEx_corr_res$inter_ind_corr_rand)
germline_GTEx_corr_stack <- stack(germline_GTEx_corr_df)
germline_GTEx_corr_stack$dataset <- "germline_GTEx"

# Plot
# Merge
all_PTC_corr_stack <- rbind(germline_TCGA_corr_stack,somatic_TCGA_corr_stack,germline_GTEx_corr_stack)
all_PTC_corr_stack$dataset <- factor(all_PTC_corr_stack$dataset)
colnames(all_PTC_corr_stack) <- c("correlation","type","dataset")
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "inter_ind_corr","Inter-individual",as.character(all_PTC_corr_stack$type))
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "intra_ind_corr","Intra-individual",all_PTC_corr_stack$type)
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "inter_ind_corr_rand","Inter-individual randomization",all_PTC_corr_stack$type)
all_PTC_corr_stack$type <- factor(ifelse(all_PTC_corr_stack$type == "intra_ind_corr_rand","Intra-individual randomization",all_PTC_corr_stack$type))
all_PTC_corr_stack$randomization <- "Observed"
all_PTC_corr_stack[grep("randomization",all_PTC_corr_stack$type),"randomization"] <- "Randomized"
all_PTC_corr_stack$type_var <- gsub(" randomization","",all_PTC_corr_stack$type)

plot_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_individual/intra_vs_inter_individual_PTC_NMDeff_correlations.png"
png(plot_path, width = 4500, height = 3500, res = 300)
p <- ggplot(data = all_PTC_corr_stack, aes(x = factor(type_var), y = correlation, fill = factor(randomization))) +
  geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(-0.1,0.4)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(aes(fill = factor(randomization)), 
              position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
              alpha = 0.05, size = 1) +
  facet_wrap(.~dataset) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("Spearman correlation") + xlab("") + ggtitle("PTC-NMDeff") +
  theme_bw(base_size = 30) + guides(fill = guide_legend(title = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=22, angle = 45, hjust = 1),
        axis.text.y = element_text(color="black", size=26),
        legend.position = "top",
        legend.text = element_text(size = 28)) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(nrow = 2)) +
  stat_compare_means(aes(group=type), size = 7, paired = TRUE,
                    label.y = c(0.75,1,1.25),
                    label = "p.format", method = "wilcox.test", hide.ns = TRUE)
print(p)
dev.off()

######## ASE ########
# germline_ASE_PTCs_NMD_efficiencies_TCGA <- germline_ASE_PTCs_NMD_efficiencies_TCGA[!germline_ASE_PTCs_NMD_efficiencies_TCGA$tissue %in% c("TCGA-ACC","TCGA-DLBC"),]
germline_TCGA_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = germline_ASE_PTCs_NMD_efficiencies_TCGA, n = 1000) 
germline_TCGA_corr_df <- data.frame(intra_ind_corr = germline_TCGA_corr_res$intra_ind_corr, 
                                    inter_ind_corr = germline_TCGA_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = germline_TCGA_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = germline_TCGA_corr_res$inter_ind_corr_rand)
germline_TCGA_corr_stack <- stack(germline_TCGA_corr_df)
germline_TCGA_corr_stack$dataset <- "germline_TCGA"
# somatic_ASE_PTCs_NMD_efficiencies_TCGA <- somatic_ASE_PTCs_NMD_efficiencies_TCGA[!somatic_ASE_PTCs_NMD_efficiencies_TCGA$tissue %in% c("TCGA-ACC","TCGA-SKCM","TCGA-TGCT","TCGA-OV","TCGA-THYM","TCGA-UCS","TCGA-KICH"),]
somatic_TCGA_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = somatic_ASE_PTCs_NMD_efficiencies_TCGA, n = 1000) 
somatic_TCGA_corr_df <- data.frame(intra_ind_corr = somatic_TCGA_corr_res$intra_ind_corr, 
                                    inter_ind_corr = somatic_TCGA_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = somatic_TCGA_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = somatic_TCGA_corr_res$inter_ind_corr_rand)
somatic_TCGA_corr_stack <- stack(somatic_TCGA_corr_df)
somatic_TCGA_corr_stack$dataset <- "somatic_TCGA"
germline_GTEx_corr_res <- PTCs_inter_intra_ind_corr(samples_NMDeff_df = germline_ASE_PTCs_NMD_efficiencies_GTEx, n = 1000) 
germline_GTEx_corr_df <- data.frame(intra_ind_corr = germline_GTEx_corr_res$intra_ind_corr, 
                                    inter_ind_corr = germline_GTEx_corr_res$inter_ind_corr,
                                    intra_ind_corr_rand = germline_GTEx_corr_res$intra_ind_corr_rand, 
                                    inter_ind_corr_rand = germline_GTEx_corr_res$inter_ind_corr_rand)
germline_GTEx_corr_stack <- stack(germline_GTEx_corr_df)
germline_GTEx_corr_stack$dataset <- "germline_GTEx"

# Plot
# Merge
all_PTC_corr_stack <- rbind(germline_TCGA_corr_stack,somatic_TCGA_corr_stack,germline_GTEx_corr_stack)
all_PTC_corr_stack$dataset <- factor(all_PTC_corr_stack$dataset)
colnames(all_PTC_corr_stack) <- c("correlation","type","dataset")
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "inter_ind_corr","Inter-individual",as.character(all_PTC_corr_stack$type))
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "intra_ind_corr","Intra-individual",all_PTC_corr_stack$type)
all_PTC_corr_stack$type <- ifelse(all_PTC_corr_stack$type == "inter_ind_corr_rand","Inter-individual randomization",all_PTC_corr_stack$type)
all_PTC_corr_stack$type <- factor(ifelse(all_PTC_corr_stack$type == "intra_ind_corr_rand","Intra-individual randomization",all_PTC_corr_stack$type))
all_PTC_corr_stack$randomization <- "Observed"
all_PTC_corr_stack[grep("randomization",all_PTC_corr_stack$type),"randomization"] <- "Randomized"
all_PTC_corr_stack$type_var <- gsub(" randomization","",all_PTC_corr_stack$type)

plot_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_individual/intra_vs_inter_individual_ASE_PTC_NMDeff_correlations.png"
png(plot_path, width = 4500, height = 3500, res = 300)
p <- ggplot(data = all_PTC_corr_stack, aes(x = factor(type_var), y = correlation, fill = factor(randomization))) +
  geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(-0.1,1)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(aes(fill = factor(randomization)), 
              position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
              alpha = 0.05, size = 1) +
  facet_wrap(.~dataset) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("Spearman correlation") + xlab("") + ggtitle("ASE PTC-NMDeff") +
  theme_bw(base_size = 30) + guides(fill = guide_legend(title = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=22, angle = 45, hjust = 1),
        axis.text.y = element_text(color="black", size=26),
        legend.position = "top",
        legend.text = element_text(size = 28)) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(nrow = 2)) +
  stat_compare_means(aes(group=type_var), size = 7, paired = TRUE,
                    label.y = c(0.75,1,1.25),
                    label = "p.format", method = "wilcox.test", hide.ns = TRUE)
print(p)
dev.off()

write.table(all_PTC_corr_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6C.txt",
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(all_PTC_corr_stack, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6C.RData")

# 3) Inter-individual NMDeff variability -- Variation metrics

######## PTC ########
germline_GTEx_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = germline_PTCs_NMD_efficiencies_GTEx)
germline_GTEx_var_df <- fix_df(germline_GTEx_var_res)
germline_GTEx_var_df$dataset <- "germline_GTEx"
germline_TCGA_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = germline_PTCs_NMD_efficiencies_TCGA)
germline_TCGA_var_df <- fix_df(germline_TCGA_var_res)
germline_TCGA_var_df$dataset <- "germline_TCGA"
somatic_TCGA_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = somatic_PTCs_NMD_efficiencies_TCGA)
somatic_TCGA_var_df <- fix_df(somatic_TCGA_var_res)
somatic_TCGA_var_df$dataset <- "somatic_TCGA"

# Plot
# Merge
all_PTC_var_stack <- rbind(germline_TCGA_var_df,somatic_TCGA_var_df,germline_GTEx_var_df)
all_PTC_var_stack$dataset <- factor(all_PTC_var_stack$dataset)
colnames(all_PTC_var_stack) <- c("variance","type","dataset")
all_PTC_var_stack$randomization <- "Observed"
all_PTC_var_stack[grep("randomization",all_PTC_var_stack$type),"randomization"] <- "Randomized"
all_PTC_var_stack$type_var <- gsub(" randomization","",all_PTC_var_stack$type)

plot_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_individual/intra_vs_inter_individual_PTC_NMDeff.png"
png(plot_path, width = 4500, height = 3500, res = 300)
p <- ggplot(data = all_PTC_var_stack, aes(x = factor(type_var), y = variance, fill = factor(randomization))) +
  geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(0,1.5)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(aes(fill = factor(randomization)), 
              position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
              alpha = 0.05, size = 1) +
  facet_wrap(.~dataset) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("Variance") + xlab("") + ggtitle("PTC-NMDeff") +
  theme_bw(base_size = 30) + guides(fill = guide_legend(title = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=22, angle = 45, hjust = 1),
        axis.text.y = element_text(color="black", size=26),
        legend.position = "top",
        legend.text = element_text(size = 28)) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(nrow = 2)) +
  stat_compare_means(aes(group=type_var), size = 7,
                    label.y = c(0.75,1,1.25), #paired = TRUE,
                    label = "p.format", method = "wilcox.test", hide.ns = TRUE)
print(p)
dev.off()

######## ASE ########
germline_GTEx_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = germline_ASE_PTCs_NMD_efficiencies_GTEx)
germline_GTEx_var_df <- fix_df(germline_GTEx_var_res)
germline_GTEx_var_df$dataset <- "germline_GTEx"
germline_TCGA_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = germline_ASE_PTCs_NMD_efficiencies_TCGA)
germline_TCGA_var_df <- fix_df(germline_TCGA_var_res)
germline_TCGA_var_df$dataset <- "germline_TCGA"
somatic_TCGA_var_res <- PTCs_inter_intra_ind_variability(samples_NMDeff_df = somatic_ASE_PTCs_NMD_efficiencies_TCGA)
somatic_TCGA_var_df <- fix_df(somatic_TCGA_var_res)
somatic_TCGA_var_df$dataset <- "somatic_TCGA"

# Plot
# Merge
all_PTC_var_stack <- rbind(germline_TCGA_var_df,somatic_TCGA_var_df,germline_GTEx_var_df)
all_PTC_var_stack$dataset <- factor(all_PTC_var_stack$dataset)
colnames(all_PTC_var_stack) <- c("variance","type","dataset")
all_PTC_var_stack$randomization <- "Observed"
all_PTC_var_stack[grep("randomization",all_PTC_var_stack$type),"randomization"] <- "Randomized"
all_PTC_var_stack$type_var <- gsub(" randomization","",all_PTC_var_stack$type)

plot_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_individual/intra_vs_inter_individual_ASE_PTC_NMDeff.png"
png(plot_path, width = 4500, height = 3500, res = 300)
p <- ggplot(data = all_PTC_var_stack, aes(x = factor(type_var), y = variance, fill = factor(randomization))) +
  geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(0,5)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(aes(fill = factor(randomization)), 
              position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
              alpha = 0.05, size = 1) +
  facet_wrap(.~dataset) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("Variance") + xlab("") + ggtitle("ASE-NMDeff") +
  theme_bw(base_size = 30) + guides(fill = guide_legend(title = "")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=22, angle = 45, hjust = 1),
        axis.text.y = element_text(color="black", size=26),
        legend.position = "top",
        legend.text = element_text(size = 28)) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(nrow = 2)) +
  stat_compare_means(aes(group=type_var), size = 7,
                    label.y = c(0.75,1,1.25), #paired = TRUE,
                    label = "p.format", method = "wilcox.test", hide.ns = TRUE)
print(p)
dev.off()

write.table(all_PTC_var_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6D.txt",
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(all_PTC_var_stack, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6D.RData")
