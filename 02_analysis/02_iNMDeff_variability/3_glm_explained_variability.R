# Libraries

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

# 1.2) somatic CNV PCA data - Individuals weights
TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))
TCGA_cancer <- "pancancer"
scale <- TRUE
center <- TRUE
alpha <- "3e-04"
num_PCs <- "100"

tryCatch({
    TCGA_CNV_PCA_ind <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer,"_sparse_PCA_ind_",alpha,"_robust_no_num_PCs_",num_PCs,".txt"),
                                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    rownames(TCGA_CNV_PCA_ind) <- gsub("\\.","-",rownames(TCGA_CNV_PCA_ind))
    # Remove PCs with 0
    cols <- colnames(TCGA_CNV_PCA_ind)[which( colSums(TCGA_CNV_PCA_ind) != 0 )]
    TCGA_CNV_PCA_ind <- TCGA_CNV_PCA_ind[,cols]
    print("PCA dimensions --> ")
    print(dim(TCGA_CNV_PCA_ind))
    nPCs <- ncol(TCGA_CNV_PCA_ind)
    # Scale PCs
    TCGA_CNV_PCA_ind <- data.frame(scale(TCGA_CNV_PCA_ind, scale = scale, center = center))
}, error = function(e){
    print(e)
    }
)
# Merge
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_PCA_ind, by.x = "sample", by.y = "row.names", all.x = TRUE)

# 2) Linear Model
NMD_methods <- c("ASE_PTC_NMD_triggering_0.2","endogenous_NMD_Consensus")
# 2.1) TCGA
#factor(batch_portion) + factor(batch_plate)

TCGA_variables_to_remove <- c("cancer_subtype","endogenous_purity","TMB","CNV_burden","sample_lib_size","POLE",
                "sex","age","vital_status","race","num_NMD_targets","CNV-PCs","all")
TCGA_continuous_var <- c("endogenous_purity","TMB","CNV_burden","sample_lib_size","age","num_NMD_targets")
TCGA_discrete_var <- c("POLE","cancer_subtype","sex","vital_status","race")
GTEx_variables_to_remove <- c("tissue","death_group","sample_lib_size","sex","age","num_NMD_targets","all")
GTEx_continuous_var <- c("sample_lib_size","num_NMD_targets")
GTEx_discrete_var <- c("age","tissue","sex","death_group")

NMDeff_variability_lm <- function(NMD_method, dataset) {

    if (dataset == "TCGA") {
        NMDeff_df <- sample_NMD_efficiencies_TCGA
        variables_to_remove <- TCGA_variables_to_remove
        continuous_var <- TCGA_continuous_var
        discrete_var <- TCGA_discrete_var
    } else if (dataset == "GTEx") {
        NMDeff_df <- sample_NMD_efficiencies_GTEx
        variables_to_remove <- GTEx_variables_to_remove
        continuous_var <- GTEx_continuous_var
        discrete_var <- GTEx_discrete_var
    }
    if (NMD_method == "endogenous_NMD_Consensus") {
        num_targets_var <- "endogenous_num_NMD_targets"
    } else {
        num_targets_var <- "ASE_num_PTCs_0.2"
    }
    variables_to_remove[variables_to_remove %in% "num_NMD_targets"] <- num_targets_var
    continuous_var[continuous_var %in% "num_NMD_targets"] <- num_targets_var

    final_res_df <- data.frame(variable = variables_to_remove, R2 = NA)
    for (remove_var in variables_to_remove) {
        # Remove variable if necessary
        if (remove_var != "all") {
            continuous_var_tmp <- continuous_var[!continuous_var %in% remove_var]
            discrete_var_tmp <- discrete_var[!discrete_var %in% remove_var]
        } else {
            continuous_var_tmp <- continuous_var
            discrete_var_tmp <- discrete_var
        }
        # Add continous variables
        continuous_var_char <- paste0(continuous_var_tmp,collapse=" + ")
        # Add discrete variables
        discrete_var_char <- ""
        for(var in discrete_var_tmp) {discrete_var_char <- paste0(discrete_var_char," + factor(",var,")")}
        # Merge
        all_var_char <- paste0(continuous_var_char,discrete_var_char)
        # Final Model
        if (remove_var == "CNV-PCs" | dataset == "TCGA") {
            lm_char <- paste0("lm(",NMD_method," ~ ",all_var_char,", data = NMDeff_df)")
        } else if (dataset == "TCGA") {
            lm_char <- paste0("lm(",NMD_method," ~ ",PCs_char," + ", all_var_char,", data = NMDeff_df)")
        } else if (dataset == "GTEx") {
            lm_char <- paste0("lm(",NMD_method," ~ ",all_var_char,", data = NMDeff_df)")
        }
        # Save results
        lm_res <- eval(parse(text=lm_char))
        lm_res_df <- summary(lm_res)
        row <- final_res_df$variable %in% remove_var
        R2 <- as.numeric(lm_res_df$adj.r.squared)
        final_res_df[row,"R2"] <- round(R2,3)
    }
    final_res_df$R2_diff <- final_res_df$R2[nrow(final_res_df)] - final_res_df$R2 

    return(final_res_df)

}

NMDeff_variability_lm(NMD_method = "endogenous_NMD_Consensus", dataset = "TCGA")
NMDeff_variability_lm(NMD_method = "ASE_PTC_NMD_triggering_0.2", dataset = "TCGA")
NMDeff_variability_lm(NMD_method = "endogenous_NMD_Consensus", dataset = "GTEx")
NMDeff_variability_lm(NMD_method = "ASE_PTC_NMD_triggering_0.2", dataset = "GTEx")
