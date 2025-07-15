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

library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggh4x)

# 1) Data
endogenous_NMD_genesets <-  c("NMD Colombo","NMD Karousis","NMD Tani","NMD Courtney","NMD Ensembl",
                      "NMD All","NMD Consensus","NMD SMG6","NMD SMG7",
                      "RandomGenes without NMD features","RandomGenes with NMD features")
ASE_NMD_genesets <- c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01","Synonymous 0.01",
                      "PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2")

# 1.1) sample NMD efficiencies TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)
# Change some columns
filter <- colnames(sample_NMD_efficiencies_TCGA) %in% c("sample","sex")
colnames(sample_NMD_efficiencies_TCGA)[filter] <- c("sample_short","sex")

# 1.2) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)
# Change some columns
filter <- colnames(sample_NMD_efficiencies_GTEx) %in% c("sample","sex")
colnames(sample_NMD_efficiencies_GTEx)[filter] <- c("sample_short","sex")

# 1.3) GTEx Rare germline variants
rare_germline_variants_name <- "PTV_Missense_CADD15_0.1perc"
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_",rare_germline_variants_name,".txt")
GTEx_variants_gnomad_allinfo <- read.csv(file = input_path, head=T,sep ="\t",stringsAsFactors = F)
# Change columns
GTEx_variants_gnomad_allinfo$Gene.refGene <- GTEx_variants_gnomad_allinfo$gene_name
filter <- colnames(GTEx_variants_gnomad_allinfo) %in% c("GTEx_sample")
colnames(GTEx_variants_gnomad_allinfo)[filter] <- c("sample_short")

# 1.4) TCGA Rare germline variants
rare_germline_variants_name <- "PTV_Missense_CADD15_0.1perc"
input_path <- paste('/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_',rare_germline_variants_name,'.txt',sep='')
TCGA_variants_gnomad_allinfo <- read.csv(file = input_path, head=T,sep ="\t",stringsAsFactors = F)

# 1.5) GTEx - PCA from common variants and other covariates
GTEx_tissues <- as.character(read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt")$V1)
GTEX_PCA_commonVariants <- c()
for (GTEx_tissue in GTEx_tissues) {
    if (GTEx_tissue == "Brain_Spinal_cord_cervical_c1") {
        GTEx_tissue <- "Brain_Spinal_cord_cervical_c-1"
    } else if (GTEx_tissue == "Cells_EBVtransformed_lymphocytes") {
        GTEx_tissue <- "Cells_EBV-transformed_lymphocytes"
    }
    print(GTEx_tissue)
    input_data <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/GTEx_Analysis_v8_eQTL_covariates/",GTEx_tissue,".v8.covariates.txt")
    tryCatch({
        error <- FALSE
        tissue_covs <- read.table( file = input_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        rownames(tissue_covs) <- tissue_covs$ID 
        tissue_covs$ID <- NULL
        tissue_covs <- data.frame(t(tissue_covs))
        tissue_covs$sample_short <- gsub("\\.","-",rownames(tissue_covs))
        tissue_covs$tissue <- GTEx_tissue
        tissue_covs <- tissue_covs[,-grep("InferredCov.*",colnames(tissue_covs))]
    },error = function(e) {
        print("Tissue with no covariates file...")
        error <- TRUE
    })
    if (length(GTEX_PCA_commonVariants) == 0) {
        GTEX_PCA_commonVariants <- tissue_covs
    } else {
        GTEX_PCA_commonVariants <- rbind(GTEX_PCA_commonVariants, tissue_covs)
    }
}
# Only 6  tissues do not have covariates (we will not use them)
colnames(GTEX_PCA_commonVariants)[1:5] <- paste0("Dim.",1:5)
rownames(GTEX_PCA_commonVariants) <- NULL
GTEX_PCA_commonVariants <- GTEX_PCA_commonVariants[!duplicated(GTEX_PCA_commonVariants),]
print(head(GTEX_PCA_commonVariants))

# Filter non-European samples and missing PCs samples
outlier_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/ancestry_and_admixture_and_PCA_outliers.txt")
all_outlier_samples <- read.table(file = outlier_path)$V1

# Remove samples
sample_NMD_efficiencies_GTEx <- sample_NMD_efficiencies_GTEx %>%
                    filter(!sample_short %in% all_outlier_samples)
GTEX_PCA_commonVariants <- GTEX_PCA_commonVariants %>%
                    select(-sex) %>%
                    filter(!sample_short %in% all_outlier_samples)
# Add PCs as covariates
sample_NMD_efficiencies_GTEx <- merge(sample_NMD_efficiencies_GTEx, GTEX_PCA_commonVariants, by = c("sample_short","tissue"), all.x = TRUE)

# 1.6) TCGA - PCA from common variants and other covariates
TCGA_PCA_commonVariants <- read.csv(file= "/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_5percent_common_ingnomAD_variants_PCA_samples_5perc_europeans.txt",
                            head=T,sep ="\t",stringsAsFactors = F) %>%
  mutate(cancer_type=sub('TCGA-','',project_id)) %>%
  mutate(sample_short=sub('_','-',sample_short)) %>%
  mutate(sample_short=sub('_','-',sample_short)) %>%
  dplyr::filter(sample_short %in% sample_NMD_efficiencies_TCGA$sample_short) %>%
  dplyr::select(sample_short,gender,age,Dim.1,Dim.2,Dim.3,Dim.4,Dim.5,Dim.6,Dim.7,Dim.8,Dim.9,Dim.10) 
TCGA_PCA_commonVariants$age[is.na(TCGA_PCA_commonVariants$age)] <- median(TCGA_PCA_commonVariants$age,na.rm=T) #replace na age with median (102 NA's)

# Add PCs as covariates
sample_NMD_efficiencies_TCGA$age <- NULL
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA, TCGA_PCA_commonVariants, by = c("sample_short"), all.x = TRUE)

# 2) Show associations
NMDeff_germ_mut_assoc <- function(NMD_method, sample_NMD_efficiencies, tissue, dataset) {
    # Fix
    NMD_method_char <- gsub("\\s","_",NMD_method)
    NMD_method_char <- gsub("\\-","_",NMD_method_char)
    colnames(sample_NMD_efficiencies)[colnames(sample_NMD_efficiencies) %in% NMD_method] <- NMD_method_char

    if (dataset == "TCGA") {
        variants_gnomad_allinfo <- TCGA_variants_gnomad_allinfo
        tissue_var <- "cancer_type"
    } else if (dataset == "GTEx") {
        variants_gnomad_allinfo <- GTEx_variants_gnomad_allinfo
        tissue_var <- "acronyms"
    }
    if ( !tissue %in% c("pancancer","pantissue")) {
        sample_NMD_efficiencies_filt <- sample_NMD_efficiencies[sample_NMD_efficiencies[,tissue_var] %in% tissue,]
    } else {
        sample_NMD_efficiencies_filt <- sample_NMD_efficiencies
    }
    
    # Create DF results
    df_res <- data.frame(gene = 1:length(genes), coefficient = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
    rownames(df_res) <- genes
    if (nrow(sample_NMD_efficiencies_filt) == 0) {return(df_res)}

    glm_char <- paste0(NMD_method_char," ~ ")

    for (gene in genes) {
        print(gene)
        # Filter variants in the gene
        rare_variants_filter <- variants_gnomad_allinfo %>% 
                        filter(Gene.refGene %in% gene)
        # Create variable of MUT vs WT samples
        sample_NMD_efficiencies_filt$germ_mut <- "WT"
        filter <- sample_NMD_efficiencies_filt$sample_short %in% rare_variants_filter$sample_short
        sample_NMD_efficiencies_filt[filter,"germ_mut"] <- "MUT"
        # At least 2 samples mutated
        n_mut <- sum(sample_NMD_efficiencies_filt$germ_mut == "MUT")
        # At least 1 male 1 female in each category
        n_sex <- sum(table(sample_NMD_efficiencies_filt$sex,sample_NMD_efficiencies_filt$germ_mut) == 0)
        # Have NMDeff values...
        n_NAs <- sum(is.na(sample_NMD_efficiencies_filt[sample_NMD_efficiencies_filt$germ_mut=="MUT",NMD_method_char])) 

        # if ( (n_mut < 2 | n_sex >= 1) | (n_mut-n_NAs < 2) ) { next}
        if ( (n_mut < 2) | (n_mut-n_NAs < 2) ) { next}
        
        if (dataset == "TCGA") {
            if (tissue %in% c("OV","BRCA","PRAD","UCEC","CESC","UCS","TGCT")) {
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + endogenous_purity + age + sample_lib_size + as.factor(cancer_subtype) + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + Dim.6, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            } else if (tissue %in% c("LAML", "THYM","DLBC")) {
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + as.factor(cancer_subtype) + age + sample_lib_size + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + Dim.6, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            } else {
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + endogenous_purity + as.factor(cancer_subtype) + as.factor(sex) + age +  sample_lib_size + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + Dim.6, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            }
        } else if (dataset == "GTEx") {
            # No sex covariate
            # if (tissue %in% c("Cervix_Ectocervix","Cervix_Endocervix","Fallopian_Tube","Kidney_Medulla","Ovary","Testis","Uterus","Vagina","Prostate")) {
            if (tissue %in% c("CVXECT","CVXEND","FLLPNT","KDNMDL","OVARY","TESTIS","UTERUS","VAGINA","PRSTTE")) {
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + factor(age) + sample_lib_size + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            # } else if (tissue %in% c("Cells_EBVtransformed_lymphocytes","Brain_Spinal_cord_cervical_c1","Bladder") ) { # Don't have PCs
            } else if (tissue %in% c("LCL","BRNSPC","BLDDER") ) { # Don't have PCs
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + factor(sex) + factor(age) + sample_lib_size, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            } else { 
                glm_model <- paste0("glm(",glm_char," + relevel(as.factor(germ_mut), ref = \"WT\") + factor(sex) + factor(age) + sample_lib_size + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5, data = sample_NMD_efficiencies_filt, family = \"gaussian\", na.action = na.exclude)")
            }        
        }

        df <- sample_NMD_efficiencies_filt[sample_NMD_efficiencies_filt$germ_mut == "MUT",]
        table(df$Dim.1)
        tryCatch({
          # Test regression
          glm_res <- eval(parse(text=glm_model))
          CI <- confint(glm_res, parm = "relevel(as.factor(germ_mut), ref = \"WT\")MUT" ,level = 0.95)
          glm_res <- summary(glm_res)
          hasMUTreg <- grep("germ_mut",rownames(glm_res$coefficients))
          if (length(hasMUTreg) != 0) {
              df_res[gene,"coefficient"] <- glm_res$coefficients[hasMUTreg,"Estimate"]
              df_res[gene,"p_value"] <- glm_res$coefficients[hasMUTreg,"Pr(>|t|)"]
              df_res[gene,"CI_2.5"] <- as.numeric(CI[1])
              df_res[gene,"CI_97.5"] <- as.numeric(CI[2])
              
          } 
        }, error = function(e){
            print(e)
            }
        )
    }
    return(df_res)
}

# 2.1) GTEx -- NMDeff vs germline mutation association
#GTEx_tissues <- c(GTEx_tissues,"pantissue")
GTEx_tissues <- unique(sample_NMD_efficiencies_GTEx$acronyms)
GTEx_glm_res_df_final <- data.frame()
genes <- c("NUP153","KDM6B","TRAP1","CRTC1","PXDN","PDIA2","LAMC1","FIG4","COL19A1","CCDC114")

for (GTEx_tissue in GTEx_tissues) {
  print(GTEx_tissue)
  if (GTEx_tissue %in% c("Cells_Leukemia_cell_line_CML")) {next}

  glm_res_df_end <- NMDeff_germ_mut_assoc(NMD_method = "NMD Consensus", sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx,
                        tissue = GTEx_tissue, dataset = "GTEx")
  glm_res_df_ASE <- NMDeff_germ_mut_assoc(NMD_method = "PTC NMD-triggering 0.2", sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx,
                        tissue = GTEx_tissue, dataset = "GTEx")
  colnames(glm_res_df_ASE) <- paste0("ASE_",colnames(glm_res_df_ASE))
  colnames(glm_res_df_end) <- paste0("END_",colnames(glm_res_df_end))
   # Save all associations
  glm_res_tissue <- merge(glm_res_df_ASE,glm_res_df_end, by = "row.names")
  glm_res_tissue$tissue <- GTEx_tissue
  if (length(GTEx_glm_res_df_final) == 0){
    GTEx_glm_res_df_final <- glm_res_tissue
  } else {
    GTEx_glm_res_df_final <- rbind(GTEx_glm_res_df_final,glm_res_tissue)
  }
}
# Fix
GTEx_glm_res_df_final$gene <- GTEx_glm_res_df_final$Row.names
GTEx_glm_res_df_final$Row.names <- NULL
GTEx_glm_res_df_final$dataset <- "GTEx"
# Significant genes -- FDR correction
FDR <- 0.1
GTEx_glm_res_df_final <- GTEx_glm_res_df_final %>%
        group_by(tissue, gene) %>%
        mutate(ASE_pvalue_FDR_adjusted = p.adjust(ASE_p_value, method = "fdr")) %>%
        mutate(ASE_significant = ifelse(ASE_pvalue_FDR_adjusted < FDR, "yes","no")) %>%
        mutate(END_pvalue_FDR_adjusted = p.adjust(END_p_value, method = "fdr")) %>%
        mutate(END_significant = ifelse(END_pvalue_FDR_adjusted < FDR, "yes","no"))

# 2.2) TCGA -- NMDeff vs germline mutation association
TCGA_tissues <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))
TCGA_glm_res_df_final <- data.frame()
genes <- c("NUP153","KDM6B","TRAP1","CRTC1","PXDN","PDIA2","LAMC1","FIG4","COL19A1","CCDC114")

for (TCGA_tissue in TCGA_tissues) {
  print(TCGA_tissue)
  if (TCGA_tissue %in% c("")) {next}
  glm_res_df_end <- NMDeff_germ_mut_assoc(NMD_method = "NMD Consensus", sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA,
                        tissue = TCGA_tissue, dataset = "TCGA")
  glm_res_df_ASE <- NMDeff_germ_mut_assoc(NMD_method = "PTC NMD-triggering 0.2", sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA,
                        tissue = TCGA_tissue, dataset = "TCGA")
  colnames(glm_res_df_ASE) <- paste0("ASE_",colnames(glm_res_df_ASE))
  colnames(glm_res_df_end) <- paste0("END_",colnames(glm_res_df_end))
   # Save all associations
  glm_res_tissue <- merge(glm_res_df_ASE,glm_res_df_end, by = "row.names")
  glm_res_tissue$tissue <- TCGA_tissue
  if (length(TCGA_glm_res_df_final) == 0){
    TCGA_glm_res_df_final <- glm_res_tissue
  } else {
    TCGA_glm_res_df_final <- rbind(TCGA_glm_res_df_final,glm_res_tissue)
  }
}
# Fix
TCGA_glm_res_df_final$gene <- TCGA_glm_res_df_final$Row.names
TCGA_glm_res_df_final$Row.names <- NULL
TCGA_glm_res_df_final$dataset <- "TCGA"
# Significant genes -- FDR correction
FDR <- 0.1
TCGA_glm_res_df_final <- TCGA_glm_res_df_final %>%
        group_by(tissue, gene) %>%
        mutate(ASE_pvalue_FDR_adjusted = p.adjust(ASE_p_value, method = "fdr")) %>%
        mutate(ASE_significant = ifelse(ASE_pvalue_FDR_adjusted < FDR, "yes","no")) %>%
        mutate(END_pvalue_FDR_adjusted = p.adjust(END_p_value, method = "fdr")) %>%
        mutate(END_significant = ifelse(END_pvalue_FDR_adjusted < FDR, "yes","no"))

# 2.3) Merge
glm_res_df_final <- rbind(GTEx_glm_res_df_final,TCGA_glm_res_df_final)

# Save
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/RGWAS_glm_effect_size_assocation.txt"
# write.table(glm_res_df_final, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# saveRDS(glm_res_df_final, output_path)
glm_res_df_final <- readRDS(output_path)

###############################

df <- glm_res_df_final[glm_res_df_final$gene == "COL19A1",]
df %>% group_by(dataset) %>%
  # filter(ASE_pvalue_FDR_adjusted < 0.1) %>%
  summarise(a = median(ASE_coefficient, na.rm = TRUE))

df <- glm_res_df_final %>% group_by(dataset,tissue, gene) %>%
          filter(ASE_pvalue_FDR_adjusted < 0.35 | END_pvalue_FDR_adjusted < 0.35) %>%
          summarise(ETG_median = mean(END_coefficient, na.rm = TRUE),
                    ASE_median = mean(ASE_coefficient, na.rm = TRUE)) %>%
          ungroup() %>%
          group_by(dataset,gene) %>%
          summarise(ETG_tissues_median = mean(ETG_median, na.rm = TRUE),
                    ASE_tissues_median = mean(ASE_median, na.rm = TRUE)) %>%
          arrange(gene)
df

###############################33

# 2.4) Forest plot
glm_res_df_final$tissue <- factor(glm_res_df_final$tissue, levels= unique(glm_res_df_final$tissue))
glm_res_df_final_stack1 <- stack(glm_res_df_final[,c("ASE_coefficient","END_coefficient")])
glm_res_df_final_stack2 <- stack(glm_res_df_final[,c("ASE_CI_2.5","END_CI_2.5")])
glm_res_df_final_stack3 <- stack(glm_res_df_final[,c("ASE_CI_97.5","END_CI_97.5")])
glm_res_df_final_stack4 <- stack(glm_res_df_final[,c("ASE_pvalue_FDR_adjusted","END_pvalue_FDR_adjusted")])
glm_res_df_final_stack <- cbind(glm_res_df_final_stack1,glm_res_df_final_stack2,glm_res_df_final_stack3,glm_res_df_final_stack4)
colnames(glm_res_df_final_stack) <- c("coefficient","NMD_method","CI_2.5_values","CI_2.5","CI_97.5_values","CI_97.5","p_value_FDR_adjusted")
glm_res_df_final_stack$tissue <- rep(glm_res_df_final$tissue,2)
glm_res_df_final_stack$gene <- rep(glm_res_df_final$gene,2)
glm_res_df_final_stack$dataset <- rep(glm_res_df_final$dataset,2)
glm_res_df_final_stack$NMD_method <- gsub("_coefficient","",glm_res_df_final_stack$NMD_method)
glm_res_df_final_stack$NMD_method <- gsub("END","ETG",glm_res_df_final_stack$NMD_method)
# sig_genes <- c("KDM6B","PXDN")
# sig_genes <- c("NUP153","CRTC1")
sig_genes <- genes
# Add SKAT-O p_values adjusted by FDR only in the discovery databases
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/RGWAS_res_filtered.txt"
RGWAS_res <- read.table(file = output_path, header = TRUE, sep = "\t")
RGWAS_res_filt <- RGWAS_res %>%
                filter(randomization == "no") %>%
                group_by(database, tissue) %>%
                mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
                filter(dataset == "PTV_Missense_CADD15_0.1perc" & variant %in% sig_genes)
RGWAS_res_filt$NMD_method <- ifelse(RGWAS_res_filt$NMD_method == "Endogenous","ETG","ASE")
RGWAS_res_filt <- RGWAS_res_filt[,c("variant","NMD_method","database","tissue","pValue","pvalue_FDR_adjusted","rho_SKAT_O")]
RGWAS_res_filt <- data.frame(RGWAS_res_filt)
colnames(RGWAS_res_filt) <- c("gene","NMD_method","database","tissue","SKATO_pValue","SKATO_p_value_FDR_adjusted","rho_SKAT_O")

# 2.4.1) TCGA
glm_res_TCGA_stack <- glm_res_df_final_stack[glm_res_df_final_stack$dataset == "TCGA",]
# Add SKAT-O pvalues
glm_res_TCGA_stack <- merge(glm_res_TCGA_stack,RGWAS_res_filt, by.x = c("gene","NMD_method","dataset","tissue"),
          by.y = c("gene","NMD_method","database","tissue"), all.x = TRUE)
# Remove outliers
# glm_res_df_final_stack <- glm_res_df_final_stack[!glm_res_df_final_stack$cancer %in% c("ACC","KICH"),]
# Order
gene_char <- "KDM6B"
tmp <- TCGA_glm_res_df_final[TCGA_glm_res_df_final$gene == gene_char,]
tmp$mean_coefficient <- as.numeric(rowMeans(tmp[,c("END_coefficient","ASE_coefficient")], na.rm = TRUE))
tissue_order <- as.character(tmp[order(tmp$mean_coefficient, decreasing = TRUE),"tissue"]$tissue)
glm_res_TCGA_stack$tissue <- factor(glm_res_TCGA_stack$tissue, levels = tissue_order)
glm_res_TCGA_stack <- glm_res_TCGA_stack[glm_res_TCGA_stack$gene %in% sig_genes,]
# Filter cancers with all NAs
tmp_vec <- table(is.na(glm_res_TCGA_stack$coefficient),glm_res_TCGA_stack$tissue)[2,] 
# cancers_to_remove <- names(tmp_vec[tmp_vec == 4])
cancers_to_remove <- names(tmp_vec[tmp_vec == 20])
glm_res_TCGA_stack_filt <- glm_res_TCGA_stack[!glm_res_TCGA_stack$tissue %in% cancers_to_remove,]

plot_TCGA <- glm_res_TCGA_stack_filt %>%
        # dplyr::filter(gene %in% sig_genes) %>%
          ggplot(aes(y = tissue, x = coefficient, color = factor(gene))) +
            geom_point(size = 10, shape = 16) +
            facet_nested( ~ gene + NMD_method, scales = "free_x") +
            geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 2, height = 0.2) +
            geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
            labs(x = "Association effect size", y = "") + scale_color_brewer(palette = "Paired") +
            theme_classic(base_size = 25) + 
            theme(plot.title = element_text(hjust = 0.5, size = 55),
                    legend.position = "top",
                    axis.text.y = element_text(size = 30),
                    axis.title.x = element_text(size = 45),
                    axis.text.x = element_text(size= 28),
                    strip.text = element_text(size = 45),
                    panel.spacing = grid::unit(3, "lines"),
                    legend.title = element_text(size = 45),
                    legend.text = element_text(size = 40)) +
            geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted < 0.1, "*", "")), 
                                    position = position_dodge(width = 1), size = 15, hjust = 0, color = "black", vjust = 0) +
            #coord_cartesian(xlim=c(-0.5,0.5)) +
            #ggtitle(paste0("Significant gene at FDR < ",FDR*100,"%")) +
            guides(color = guide_legend(title = "PCs", override.aes = list(size = 16)), size = "none") #+
            #scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.1)) #+
            #scale_x_continuous(breaks = c(-0.2,0,0.2))

glm_res_TCGA_stack_final <- glm_res_TCGA_stack_filt[glm_res_TCGA_stack_filt$gene == "KDM6B",]
# Save
write.table(glm_res_TCGA_stack_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D_1.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_TCGA_stack_final, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D_1.RData")

# Save
write.table(glm_res_TCGA_stack_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/SuppFig17A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_TCGA_stack_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/SuppFig17A.RData")

# 2.4.2) GTEx
RGWAS_res_filt <- merge(RGWAS_res_filt,unique(sample_NMD_efficiencies_GTEx[,c("tissue","acronyms")]),
                       by.x = "tissue", by.y ="tissue", all.x = TRUE)
glm_res_GTEx_stack <- glm_res_df_final_stack[glm_res_df_final_stack$dataset == "GTEx",]
glm_res_GTEx_stack <- merge(glm_res_GTEx_stack,RGWAS_res_filt, by.x = c("gene","NMD_method","dataset","tissue"),
          by.y = c("gene","NMD_method","database","acronyms"), all.x = TRUE)

# Remove outliers
# glm_res_df_final_stack <- glm_res_df_final_stack[!glm_res_df_final_stack$cancer %in% c("ACC","KICH"),]
# Order
gene_char <- "KDM6B"
tmp <- GTEx_glm_res_df_final[GTEx_glm_res_df_final$gene == gene_char,]
tmp$mean_coefficient <- as.numeric(rowMeans(tmp[,c("END_coefficient","ASE_coefficient")], na.rm = TRUE))
tissue_order <- as.character(tmp[order(tmp$mean_coefficient, decreasing = TRUE),"tissue"]$tissue)
glm_res_GTEx_stack$tissue <- factor(glm_res_GTEx_stack$tissue, levels = tissue_order)

glm_res_GTEx_stack <- glm_res_GTEx_stack[glm_res_GTEx_stack$gene %in% sig_genes,]

# Filter cancers with all NAs
tmp_vec <- table(is.na(glm_res_GTEx_stack$coefficient),glm_res_GTEx_stack$tissue)[2,] 
cancers_to_remove <- names(tmp_vec[tmp_vec == 20])
# cancers_to_remove <- names(tmp_vec[tmp_vec == 4])
glm_res_GTEx_stack_filt <- glm_res_GTEx_stack[!glm_res_GTEx_stack$tissue %in% cancers_to_remove,]

plot_GTEx <- glm_res_GTEx_stack_filt %>%
        # dplyr::filter(gene %in% sig_genes) %>%
          ggplot(aes(y = tissue, x = coefficient, color = factor(gene))) +
            geom_point(size = 10, shape = 16) +
            facet_nested( ~ gene + NMD_method, scales = "free_x") +
            geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 2, height = 0.2) +
            geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
            labs(x = "Association effect size", y = "") + scale_color_brewer(palette = "Paired") +
            theme_classic(base_size = 25) + 
            theme(plot.title = element_text(hjust = 0.5, size = 55),
                    legend.position = "top",
                    axis.text.y = element_text(size = 30),
                    axis.title.x = element_text(size = 45),
                    axis.text.x = element_text(size= 28),
                    strip.text = element_text(size = 45),
                    panel.spacing = grid::unit(3, "lines"),
                    legend.title = element_text(size = 45),
                    legend.text = element_text(size = 40)) +
            geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted < 0.1, "*", "")), 
                                    position = position_dodge(width = 1), size = 15, hjust = 0, color = "black", vjust = 0) +
            #coord_cartesian(xlim=c(-0.5,0.5)) +
            #ggtitle(paste0("Significant gene at FDR < ",FDR*100,"%")) +
            guides(color = guide_legend(title = "PCs", override.aes = list(size = 16)), size = "none") #+
            #scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.1)) #+
            #scale_x_continuous(breaks = c(-0.2,0,0.2))
glm_res_GTEx_stack_final <- glm_res_GTEx_stack_filt[glm_res_GTEx_stack_filt$gene == "KDM6B",]
# Save
write.table(glm_res_GTEx_stack_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D_2.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_GTEx_stack_final, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D_2.RData")
# Save
write.table(glm_res_GTEx_stack_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/SuppFig17B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_GTEx_stack_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/SuppFig17B.RData")

forest_plot <- plot_grid(plot_TCGA, plot_GTEx, nrow = 1, ncol = 2, rel_widths = c(0.4,0.6))
final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/nogenelist/effect_size/TCGA_GTEx_hits_assoc_effect_size.png")
png(final_figure_path, width = 14000, height = 6500, res = 300)
print(forest_plot)
dev.off()








