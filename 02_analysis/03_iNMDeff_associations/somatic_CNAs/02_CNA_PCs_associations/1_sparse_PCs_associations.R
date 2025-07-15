replicated_hits <- function(discovery, validation, FDR) {
  # Split dataframe
  glm_res_df_final$
  # Discovery set --> ASE
  pval_col <- grep("p_value",colnames(discovery))
  discovery$pvalue_adjust <- p.adjust(discovery[,pval_col],method = "fdr")
  PCs_sig <- rownames(discovery[which(discovery$pvalue_adjust < FDR),])
  if (length(PCs_sig) == 0 ){
    print("No significant PCs")
    return(NA)
  } else {
    print(paste0("Significant PCs in discovery --> "))
    print(PCs_sig)
  }
  colnames(discovery)[-1] <- paste0(colnames(discovery)[-1],"_discovery")
  # Validation set --> Endogenous
  PCs_validated <- which(rownames(validation) %in% PCs_sig)
  validation <- validation[PCs_validated,]
  pval_col <- grep("p_value",colnames(validation))
  validation$pvalue_adjust <- p.adjust(validation[,pval_col],method = "fdr")
  colnames(validation)[-1] <- paste0(colnames(validation)[-1],"_validation")
  validation <- validation[which(validation$pvalue_adjust_validation < FDR),]
  if (nrow(validation) == 0 ){
    print("No PCs validated")
    return(NA)
  }
  # Merge discovery and validation
  results <- merge(validation,discovery, by = "row.names", all.x = TRUE)
  results <- results[order(results$pvalue_adjust_discovery),]
  rownames(results) <- results$Row.names
  results$Row.names <- NULL
  results$PC_names <- rownames(results)
  return(results)
}

create_CNV_PCs_glm_char <- function(n) {
  PCs <- c()
  for (i in 1:as.numeric(n)) {
      PC <- paste0("Dim.",i)
      PCs <- c(PCs,PC)
  }
  CNV_PCs_glm_char <- paste(PCs,collapse=" + ")
  return(CNV_PCs_glm_char)
}

NMDeff_CNV_PCs_assoc <- function(NMD_geneset, sample_NMD_efficiencies_TCGA, TCGA_cancer) {
  PCs_names <- colnames(sample_NMD_efficiencies_TCGA)[grep("Dim.",colnames(sample_NMD_efficiencies_TCGA))]
  if ( TCGA_cancer != "pancancer") {
    sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA[sample_NMD_efficiencies_TCGA$cancer_type %in% TCGA_cancer,]
  } else {
    sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA
  }
  # Create DF results
  df_res <- data.frame(CNV_PC = 1:length(PCs_names), coefficient = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
  rownames(df_res) <- PCs_names
  glm_char <- paste0(NMD_geneset," ~ ")
  for (CNV_PC in PCs_names) {
    print(CNV_PC)
    if (TCGA_cancer %in% c("OV","BRCA","PRAD","UCEC","CESC","UCS","TGCT")) {
      glm_model <- paste0("glm(",glm_char," endogenous_purity + as.factor(cancer_subtype) + age +  sample_lib_size + ",CNV_PC,", data = sample_NMD_efficiencies_TCGA_filt, family = \"gaussian\", na.action = na.exclude)")
    } else if (TCGA_cancer %in% c("LAML", "THYM","DLBC")) {
      glm_model <- paste0("glm(",glm_char," as.factor(cancer_subtype) + age +  sample_lib_size + ",CNV_PC,", data = sample_NMD_efficiencies_TCGA_filt, family = \"gaussian\", na.action = na.exclude)")    
    } else {
      glm_model <- paste0("glm(",glm_char," endogenous_purity + as.factor(cancer_subtype) + as.factor(sex) + age +  sample_lib_size + ",CNV_PC,", data = sample_NMD_efficiencies_TCGA_filt, family = \"gaussian\", na.action = na.exclude)")
    }
    # 5) Test regression
    glm_res <- eval(parse(text=glm_model))
    CI <- confint(glm_res, parm = CNV_PC ,level = 0.95)
    glm_res <- summary(glm_res)
    hasCNVreg <- grep("Dim.",rownames(glm_res$coefficients))
    if (length(hasCNVreg) != 0) {
      df_res[CNV_PC,"coefficient"] <- glm_res$coefficients[hasCNVreg,"Estimate"]
      df_res[CNV_PC,"p_value"] <- glm_res$coefficients[hasCNVreg,"Pr(>|t|)"]
      df_res[CNV_PC,"CI_2.5"] <- as.numeric(CI[1])
      df_res[CNV_PC,"CI_97.5"] <- as.numeric(CI[2])
      
    } 
  }
  return(df_res)
}

qqplot <- function(glm_res_df, NMD_method_VAF, TCGA_cancer) {
  num_PCs <- nrow(glm_res_df)
  qqplot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/qqplots/",TCGA_cancer,"_",NMD_method_VAF,"_qqplot_PCs_",num_PCs,".png")
  # 1) QQ-plots
  png(qqplot, width = 3500, height = 2500, res = 300)
  p_values <- glm_res_df[,paste0(NMD_method_VAF,"_p_value")]
  # Lambda calculation
  chisq <- qchisq(1-p_values,1)
  lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
  p <- qqPlot(p_values, truncate = FALSE, ylim = NULL, thinThreshold = NULL, 
              ci=TRUE, main = paste0(TCGA_cancer," -- QQ-plot for CNV-PCs associations. Lambda = ",lambda)) #col=as.character(df_ass_final$color_qqplot))
  # legend(x = "topleft", legend=c("cancer_gene","NMD_factor","NMD_related"),
  #             fill = 1:3, cex=1.2)        
  print(p)
  dev.off()
  return(lambda)
}

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

library("GWASTools")
library("ggrepel")
library("dplyr")
library("ggplot2")
#library("Karyoplotter")
library("qvalue")
library("ggbreak")
library("scales")
library("patchwork")
library("readxl")
library("viridis")

# 1) Data
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")

# 1.1) NMD efficiencies
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# 1.2) somatic CNV PCA data
TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
#TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
TCGA_cancers <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))
# 1.2.1) Individuals weights
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

# 1.3) TCGA/PCAWG pancancer CNV signatures (Nature 2022 Geoff)

input_data <- "/g/strcombio/fsupek_cancer1/gpalou/CNV_signatures/EMS150036-supplement-Supp_Tables_15_22.xlsx"
pancancer_TCGA_CNV_signatures <- read_excel(input_data, sheet = "ST_18_TCGA_Activities_raw", 
                            col_names = TRUE, col_types = NULL, na = "", skip = 0)
pancancer_TCGA_CNV_signatures <- data.frame(pancancer_TCGA_CNV_signatures)
colnames(pancancer_TCGA_CNV_signatures)[1] <- "sample"
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA, pancancer_TCGA_CNV_signatures, by.x = "sample", by.y = "sample", all.x = TRUE)

#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################

# 2) NMDeff vs sparse CNV-PCs associations
TCGA_cancers <- c(TCGA_cancers,"pancancer")
glm_res_df_final <- data.frame()
for (TCGA_cancer in TCGA_cancers) {
  print(TCGA_cancer)
  if (TCGA_cancer == "LAML") {next}
  glm_res_df_ASE <- NMDeff_CNV_PCs_assoc(NMD_geneset = "ASE_PTC_NMD_triggering_0.2", sample_NMD_efficiencies_TCGA = sample_NMD_efficiencies_TCGA, TCGA_cancer = TCGA_cancer)
  glm_res_df_end <- NMDeff_CNV_PCs_assoc(NMD_geneset = "endogenous_NMD_Consensus", sample_NMD_efficiencies_TCGA = sample_NMD_efficiencies_TCGA, TCGA_cancer = TCGA_cancer)
  colnames(glm_res_df_ASE) <- paste0("ASE_",colnames(glm_res_df_ASE))
  colnames(glm_res_df_end) <- paste0("END_",colnames(glm_res_df_end))
  # QQ-plots
  lambda_ASE <- qqplot(glm_res_df = glm_res_df_ASE, NMD_method_VAF = "ASE", TCGA_cancer = TCGA_cancer)
  lambda_END <- qqplot(glm_res_df = glm_res_df_end, NMD_method_VAF = "END", TCGA_cancer = TCGA_cancer)
  # Save all associations
  glm_res_cancer <- merge(glm_res_df_ASE,glm_res_df_end, by = "row.names")
  glm_res_cancer$cancer <- TCGA_cancer
  if (length(glm_res_df_final) == 0){
    glm_res_df_final <- glm_res_cancer
  } else {
    glm_res_df_final <- rbind(glm_res_df_final,glm_res_cancer)
  }
}
# Fix
glm_res_df_final$CNV_PCs_full <- glm_res_df_final$Row.names
glm_res_df_final$Row.names <- NULL
glm_res_df_final$CNV_PCs <- glm_res_df_final$ASE_CNV_PC

# Save dataframe
file_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/NMDeff_CNV_PCs_associations.txt"
# write.table(glm_res_df_final, file = file_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
glm_res_df_final <- read.table(file = file_path, header = TRUE, sep = "\t")

# 3) Replicated hits across cancers and pancancer
# ASE as discovery and Endogenous as validation
FDR <- 0.1
replicated_hits <- glm_res_df_final %>%
        group_by(cancer) %>%
        mutate(ASE_pvalue_FDR_adjusted = p.adjust(ASE_p_value, method = "fdr")) %>%
        mutate(END_pvalue_FDR_adjusted = p.adjust(END_p_value, method = "fdr")) %>%
        filter(ASE_pvalue_FDR_adjusted < FDR) %>%
        mutate(END_pvalue_FDR_adjusted_validation = p.adjust(END_p_value, method = "fdr")) %>%
        filter(END_pvalue_FDR_adjusted < FDR)
data.frame(replicated_hits)
file_path2 <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/NMDeff_CNV_PCs_replicated_associations_FDR_",FDR,".txt")
# write.table(replicated_hits, file = file_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
replicated_hits <- read.table(file = file_path2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

########### cancer type specific hits with no replication #######

FDR <- 0.1
replicated_hits <- glm_res_df_final %>%
        group_by(cancer) %>%
        mutate(ASE_pvalue_FDR_adjusted = p.adjust(ASE_p_value, method = "fdr")) %>%
        mutate(END_pvalue_FDR_adjusted = p.adjust(END_p_value, method = "fdr")) %>%
        filter(ASE_pvalue_FDR_adjusted < FDR) %>%
        mutate(END_pvalue_FDR_adjusted_validation = p.adjust(END_p_value, method = "fdr")) 
data.frame(replicated_hits[,c("CNV_PCs","cancer","ASE_coefficient","END_coefficient")])

#######################################

# 3.1) Volcanot plot with all CNA-PCs
replicated_CNV_PCs <- replicated_hits[replicated_hits$CNV_PCs & !replicated_hits$cancer %in% c("PCPG","UCEC"),]

glm_res_df_final$replicated_hits <- "Not replicated"
glm_res_df_final[glm_res_df_final$CNV_PCs %in% replicated_CNV_PCs$CNV_PCs,"replicated_hits"] <- 
          glm_res_df_final[glm_res_df_final$CNV_PCs %in% replicated_CNV_PCs$CNV_PCs,"CNV_PCs"]
glm_res_df_final$replicated_hits_labels <- NA
filter <- glm_res_df_final$CNV_PCs %in% replicated_CNV_PCs$CNV_PCs & glm_res_df_final$cancer %in% "pancancer"
glm_res_df_final[filter,"replicated_hits_labels"] <- 
          paste0("pancancer_",glm_res_df_final[filter,"CNV_PCs"])
# Colors
# Extract colors from the Brewer palette
brewer_colors <- scales::brewer_pal(palette = "Set1")(length(unique(glm_res_df_final$replicated_hits)))
# Modify the color for the "other" group
names(brewer_colors) <- unique(glm_res_df_final$replicated_hits)
brewer_colors[names(brewer_colors) == "Not replicated"] <- "#787B80"
# Manual change
brewer_colors[names(brewer_colors) == "71"] <- "#BDCB14"

NMD_method <- "END"
volcano <- ggplot(glm_res_df_final, aes(x = eval(parse(text=paste0(NMD_method,"_coefficient"))), 
              y = -log10(eval(parse(text=paste0(NMD_method,"_p_value")))), color = replicated_hits )) +
                  xlim(-0.5,0.5) +
                  geom_point(data = glm_res_df_final %>% filter(replicated_hits == "Not replicated"),
                      size = 2, alpha = 0.3, color = "grey") + 
                  geom_point(data = glm_res_df_final %>% filter(replicated_hits != "Not replicated"),
                      size = 2, alpha = 0.75) +
                  geom_label_repel(aes(label=replicated_hits_labels),hjust=0.5, vjust=0.5, size = 2, 
                      max.overlaps = nrow(glm_res_df_final), alpha = 1, key_glyph = draw_key_blank) +
                  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
                  geom_vline(xintercept=0, linetype="dashed", color = "red") +
                  labs(x = "Association coefficients", y = "-log10(p-value)", 
                        title = NMD_method, color = "Replicated CNA-PCs") +
                  scale_color_manual(values = brewer_colors) + # Use the modified colors
                  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
                  theme_classic() +
                  theme(legend.position="top") +
                  guides(color = guide_legend(override.aes = list(size=4)))
volcanot_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/CNV_PCs_associations_volcano_plot.png")
png(volcanot_plot, width = 1500, height = 1500, res = 300)
print(volcano)
dev.off()

# Save
write.table(glm_res_df_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3A_1.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_df_final, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3A_1.RData")

# 4) Significant CNV-PCs forest plot
#glm_res_df_final <- glm_res_df_final[order(glm_res_df_final$PC_names, glm_res_df_final$ASE_coefficient, decreasing = FALSE),]

# Swap direction of coefficients for Dim.3
numeric_cols <- c("ASE_coefficient","ASE_CI_2.5","ASE_CI_97.5","END_coefficient","END_CI_2.5","END_CI_97.5")
glm_res_df_final[glm_res_df_final$ASE_CNV_PC %in% 3,numeric_cols] <- -glm_res_df_final[glm_res_df_final$ASE_CNV_PC %in% 3,numeric_cols]

glm_res_df_final$cancer <- factor(glm_res_df_final$cancer, levels= unique(glm_res_df_final$cancer))
glm_res_df_final_stack1 <- stack(glm_res_df_final[,c("ASE_coefficient","END_coefficient")])
glm_res_df_final_stack2 <- stack(glm_res_df_final[,c("ASE_CI_2.5","END_CI_2.5")])
glm_res_df_final_stack3 <- stack(glm_res_df_final[,c("ASE_CI_97.5","END_CI_97.5")])
glm_res_df_final_stack <- cbind(glm_res_df_final_stack1,glm_res_df_final_stack2,glm_res_df_final_stack3)
colnames(glm_res_df_final_stack) <- c("coefficient","NMD_method","CI_2.5_values","CI_2.5","CI_97.5_values","CI_97.5")
glm_res_df_final_stack$cancer <- rep(glm_res_df_final$cancer,2)
glm_res_df_final_stack$CNV_PCs <- rep(glm_res_df_final$CNV_PCs,2)
glm_res_df_final_stack$NMD_method <- gsub("_coefficient","",glm_res_df_final_stack$NMD_method)
glm_res_df_final_stack$NMD_method <- gsub("END","ETG",glm_res_df_final_stack$NMD_method)

sig_CNV_PCs <- unique(replicated_hits[replicated_hits$cancer == "pancancer","CNV_PCs"])
# Remove PCs not interesed in
remove_PCs <- c(11,27,71)
sig_CNV_PCs <- sig_CNV_PCs[!sig_CNV_PCs %in% remove_PCs]
# Remove outliers
glm_res_df_final_stack <- glm_res_df_final_stack[!glm_res_df_final_stack$cancer %in% c("ACC","KICH"),]
# Order
tmp <- glm_res_df_final[glm_res_df_final$CNV_PCs == 3,]
# tmp$Mean_coefficient <- rowMeans(tmp[,c("END_coefficient","ASE_coefficient")])
cancer_order <- as.character(tmp[order(tmp$END_coefficient, decreasing = TRUE),"cancer"])
glm_res_df_final_stack$cancer <- factor(glm_res_df_final_stack$cancer, levels = cancer_order)

# 4.1) Selected CNA-PCs
forest <- glm_res_df_final_stack %>%
        filter(CNV_PCs %in% sig_CNV_PCs) %>%
          ggplot(aes(y = cancer, x = coefficient, color = factor(CNV_PCs))) +
          geom_point(size = 10, shape = 16) +
          facet_wrap( ~ CNV_PCs + NMD_method, scales = "free_x", nrow = 1) +
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
          #coord_cartesian(xlim=c(-0.5,0.5)) +
          ggtitle(paste0("Significant pan-cancer CNA-PCs at FDR < ",FDR*100,"%")) +
          guides(color = guide_legend(title = "PCs", override.aes = list(size = 16)), size = "none") +
          scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.1)) #+
          #scale_x_continuous(breaks = c(-0.2,0,0.2))

forest_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/CNV_PCs_associations_forest_pancancer_sig_PCs_FDR_",FDR,".png")
png(forest_plot, width = 8000, height = 5500, res = 300)
print(forest)
dev.off()

# Save
write.table(glm_res_df_final_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig10/SuppFig10A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_df_final_stack, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig10/SuppFig10A.RData")

# 4.2) CNA-PC3 only
sig_CNV_PCs <- unique(replicated_hits[replicated_hits$cancer == "pancancer","CNV_PCs"])
# Remove PCs not interesed in
keep_PCs <- c(3)
sig_CNV_PCs <- sig_CNV_PCs[sig_CNV_PCs %in% keep_PCs]
# Remove outliers
glm_res_df_final_stack <- glm_res_df_final_stack[!glm_res_df_final_stack$cancer %in% c("ACC","KICH"),]
# Order
tmp <- glm_res_df_final[glm_res_df_final$CNV_PCs == 3,]
# tmp$Mean_coefficient <- rowMeans(tmp[,c("END_coefficient","ASE_coefficient")])
cancer_order <- as.character(tmp[order(tmp$END_coefficient, decreasing = TRUE),"cancer"])
glm_res_df_final_stack$cancer <- factor(glm_res_df_final_stack$cancer, levels = cancer_order)
glm_res_df_final_stack$NMD_method <- factor(glm_res_df_final_stack$NMD_method, levels = c("ETG","ASE"))

forest <- glm_res_df_final_stack %>%
        filter(CNV_PCs %in% sig_CNV_PCs) %>%
          ggplot(aes(y = cancer, x = coefficient, color = factor(CNV_PCs))) +
          geom_point(size = 10, shape = 16) +
          facet_wrap( ~ NMD_method, scales = "free_x", nrow = 1) +
          geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 2, height = 0.2) +
          geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
          labs(x = "Effect size", y = "") + scale_color_brewer(palette = "Paired") +
          theme_classic(base_size = 25) + 
          theme(plot.title = element_text(hjust = 0.5, size = 63),
                legend.position = "none",
                axis.text.y = element_text(size = 40),
                axis.title.x = element_text(size = 65),
                axis.text.x = element_text(size= 47),
                strip.text = element_text(size = 55),
                panel.spacing = grid::unit(3, "lines")) +
          ggtitle(paste0("CNA-PC3 associations to iNMDeff")) +
          scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.1))
          #scale_x_continuous(breaks = c(-0.1,0,0.1))

forest_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/significant_PCs/CNA_PC3_associations_forest_pancancer_sig_PCs_FDR_",FDR,".png")
png(forest_plot, width = 5000, height = 5500, res = 300)
print(forest)
dev.off()

# Save
write.table(glm_res_df_final_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(glm_res_df_final_stack, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3A.RData")

