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

sample_NMDeff_corr <- function(sample_NMD_efficiencies, tissue = NULL, ASE_VAF, dataset, cov_corr = NULL) {

  output_list <- list()
  if (!is.null(tissue)) {
    char <- "tissues/"
    if (dataset == "TCGA") {
      sample_NMD_efficiencies <- sample_NMD_efficiencies[which(sample_NMD_efficiencies$cancer_type_strat %in% tissue),]
      #sample_NMD_efficiencies$tissue_selection <- factor(sample_NMD_efficiencies$MSI_status)
      sample_NMD_efficiencies$tissue_selection <- "All"
    } else if (dataset == "GTEx") {
      sample_NMD_efficiencies <- sample_NMD_efficiencies[which(sample_NMD_efficiencies$acronyms %in% tissue),]
      #sample_NMD_efficiencies$tissue_selection <- factor(tissue)
      sample_NMD_efficiencies$tissue_selection <- "All"
    }
    if ( length(levels(sample_NMD_efficiencies$tissue_selection)) == 0 ) {
      sample_NMD_efficiencies$tissue_selection <- "All"
    }
  } else {
    tissue_selection <- NULL
    if (dataset == "TCGA") {
      tissue <- "Pan-cancer"
    } else if (dataset == "GTEx") {
      tissue <- "Pan-tissue"
    }
    char <- "/"
  }

  # NMDeff scatterplot correlations

  formula <- y ~ x
  method1 <- paste0("ASE_PTC_NMD_triggering_0.2")
  method2 <- paste0("endogenous_NMD_Consensus")
  # Check sample size
  m1_ss <- sum(!is.na(sample_NMD_efficiencies[,method1]))
  m2_ss <- sum(!is.na(sample_NMD_efficiencies[,method2]))
  if (m1_ss >= m2_ss) {
    sample_size <- m2_ss
  } else {
    sample_size <- m1_ss
  }
  # Plot
  p <- ggplot(data = sample_NMD_efficiencies, mapping = aes(x = eval(parse(text = method1)), y = eval(parse(text = method2)), color = tissue_selection)) +
    geom_point(alpha = 0.3, color = "black") + ggtitle(paste0(tissue,", n = ",sample_size)) +
    #stat_density_2d(geom = "polygon", alpha=0.3,size = 0.10, aes(alpha = ..level.., fill = tissue_selection), bins = 20, color = "blue") +
    geom_smooth(method = "lm", formula = formula, se = TRUE, size = 1, color = "#69b3a2") +
    xlab(paste0("ASE NMD efficiency")) + ylab(paste0("Endogenous NMD efficiency")) + coord_cartesian(xlim=c(-4,4),ylim=c(-4,4))+
    theme_classic() + scale_colour_discrete(na.translate = F) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=18, face="bold"),
          axis.text.x = element_text(color="black", size=24),
          axis.text.y = element_text(color="black", size=24),
          axis.title.y = element_text(color="black", size=18, face="bold"),
          legend.position='top')  + #+ xlim(c(0,3.5)) + ylim(c(0,3.5)) +
    stat_cor(p.accuracy = NULL, r.accuracy = 0.001, size = 5, label.x = -3.75, label.y = 3)
  #> `geom_smooth()` using formula 'y ~ x'
  png(paste0(corr_output_path,"/sample_level/",char,dataset,"_",tissue,"_methods_correlations_ASE_VAF_",ASE_VAF,".png"), 
        width = 2500, height = 1750, res = 300)
  print(p)  
  dev.off()

  # NMDeff + covariates matrix of correlation 
  if (tissue %in% c("Pan-cancer","Pan-tissue")) {
    if (dataset == "TCGA") {
      cols <- c("0.2|.*_Consensus.*|endogenous|.*non_NMD.*|num_NMD_targets|TMB|TIB|TNB|CNV_burden|sample_lib_size|MSI_score|days|age|sex|purity$")
      sample_NMD_efficiencies$sex <- ifelse(sample_NMD_efficiencies$sex == "male",0,1)
    } else if (dataset == "GTEx") {
      cols <- c("0.2|.*_Consensus.*|endogenous|.*non_NMD.*|num_NMD_targets|sample_lib_size|age|death_group|sex|purity")
      sample_NMD_efficiencies$age <- as.numeric(substr(sample_NMD_efficiencies$age,1,2))
    }
    # Cols and rename
    num_cols <- grep(cols,colnames(sample_NMD_efficiencies))
    sample_NMD_efficiencies_filt <- sample_NMD_efficiencies[,num_cols]
    rem_cols <- grep("randomized",colnames(sample_NMD_efficiencies_filt))
    sample_NMD_efficiencies_filt <- sample_NMD_efficiencies_filt[,-rem_cols]
    colnames(sample_NMD_efficiencies_filt)[colnames(sample_NMD_efficiencies_filt) %in% c("endogenous_LF","endogenous_purity")] <- c("Leukocyte Fraction","Purity")
    colnames(sample_NMD_efficiencies_filt) <- gsub("_|0.2"," ",colnames(sample_NMD_efficiencies_filt))
    # Correlation matrix
    corr <- cor(sample_NMD_efficiencies_filt, use = "pairwise.complete.obs", method = "pearson")
    #corr <- corr[rownames(corr) %in% c("endogenous NMD Consensus","ASE PTC NMD triggering  "),]
    if (!is.null(cov_corr)) {
      return(corr)
    }
    png(paste0(corr_output_path,"/sample_level/",dataset,"_matrix_correlation_all.png"), width = 5000, height = 3500, res = 300)
    # p <- corrplot(as.matrix(corr), is.corr = FALSE)
    p <- corrplot(as.matrix(corr), type = "upper", method = "pie",
              #p.mat = as.matrix(pvas_mat),
              pch = 10,
              title = "",
              tl.cex = 0.9,
              tl.srt= 30,
              pch.cex = 15,
              #sig.level = 0.05
              insig = "label_sig",
              #mar=c(0,0,0,0),tl.offset = 1,
              diag = FALSE,
              order = 'alphabet')
    print(p)
    dev.off()
  } else {
    # Return correlations
    cor_res <- cor.test(sample_NMD_efficiencies$ASE_PTC_NMD_triggering_0.2,sample_NMD_efficiencies$endogenous_NMD_Consensus, 
              na.rm = TRUE, method = "pearson")
    cor_res_df <- data.frame(tissue = tissue, R = as.numeric(cor_res$estimate), p_value = cor_res$p.value, 
                sample_size = as.numeric(cor_res$parameter), R_conf_int_low = cor_res$conf.int[1], R_conf_int_high = cor_res$conf.int[2])
    return(cor_res_df)
  }
}

PTC_NMDeff_corr <- function(PTCs_ASE_NMD_efficiencies_TCGA, tissue = NULL) {

  if(!is.null(tissue)) {
    PTCs_ASE_NMD_efficiencies_TCGA <- PTCs_ASE_NMD_efficiencies_TCGA[which(PTCs_ASE_NMD_efficiencies_TCGA$cancer_type_strat %in% tissue),]
    PTCs_ASE_NMD_efficiencies_TCGA$tissue_selection <- factor(PTCs_ASE_NMD_efficiencies_TCGA$MSI_status)
    if ( length(levels(PTCs_ASE_NMD_efficiencies_TCGA$tissue_selection)) == 0 ) {
      PTCs_ASE_NMD_efficiencies_TCGA$tissue_selection <- "All"
    }
  } else {
    PTCs_ASE_NMD_efficiencies_TCGA$tissue_selection <- "All"
    tissue <- "all_cancers"
  }

  formula <- y ~ x
  png(paste0(corr_output_path,"/PTC_level/",tissue,"_methods_correlations.png"), width = 3500, height = 2000, res = 300)
  p <- ggplot(data = PTCs_ASE_NMD_efficiencies_TCGA, mapping = aes(x = NMD_efficiency_TPM, y = ASE_NMD_efficiency_TPM, color = tissue_selection)) +
    geom_point(alpha = 0.5) + ggtitle(tissue) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
    xlab(paste0("PTCs_stopgain_NMD_triggering")) + ylab(paste0("ASE_stopgain")) +
    theme_classic() + scale_colour_discrete(na.translate = F) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          legend.position='right') #+ xlim(c(0,3.5)) + ylim(c(0,3.5))
  print(p + stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5))
  dev.off()

}

tissue_ranking_corr_matrix <- function(ASE_VAF, TCGA_GTEx_tissues, sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx) {

  # 1) Sample level medians
  NMD_genesets <- c("ASE_PTC_NMD_triggering","endogenous_NMD_Consensus")
  # TCGA
  samples_NMDeff_median_TCGA_samples <- data.frame(cancer_type_strat = levels(sample_NMD_efficiencies_TCGA$cancer_type_strat))
  for (NMD_method in NMD_genesets) {
    if (NMD_method == "ASE_PTC_NMD_triggering") {
      NMD_geneset <- paste0("ASE_PTC_NMD_triggering_",ASE_VAF)
    } else {
      NMD_geneset <- NMD_method
    }
    samples_NMDeff_median <- aggregate(eval(parse(text = NMD_geneset)) ~ cancer_type_strat, data=sample_NMD_efficiencies_TCGA,median)
    samples_NMDeff_median <- samples_NMDeff_median[order(samples_NMDeff_median$cancer_type_strat, decreasing = FALSE),]
    colnames(samples_NMDeff_median)[2] <- NMD_method
    samples_NMDeff_median_TCGA_samples <- merge(samples_NMDeff_median_TCGA_samples,samples_NMDeff_median, by = c("cancer_type_strat"), all.x = TRUE)
  }
  #colnames(samples_NMDeff_median_TCGA_samples) <- c("cancer_type_strat","TCGA_NMDeff_ASE","TCGA_NMDeff_endogenous","TCGA_NMDeff_PTCs","TCGA_NMDeff_mean","TCGA_NMDeff_PCA")
  colnames(samples_NMDeff_median_TCGA_samples) <- c("cancer_type_strat","TCGA_NMDeff_ASE","TCGA_NMDeff_endogenous")
  # GTEx
  samples_NMDeff_median_GTEx_samples <- data.frame(tissue = unique(sample_NMD_efficiencies_GTEx$acronyms))
  for (NMD_method in NMD_genesets) {
    if (NMD_method == "ASE_PTC_NMD_triggering") {
      NMD_geneset <- paste0("ASE_PTC_NMD_triggering_",ASE_VAF)
    } else {
      NMD_geneset <- NMD_method
    }
    samples_NMDeff_median <- aggregate(eval(parse(text = NMD_geneset)) ~ acronyms, data=sample_NMD_efficiencies_GTEx,median)
    samples_NMDeff_median <- samples_NMDeff_median[order(samples_NMDeff_median$acronyms, decreasing = FALSE),]
    colnames(samples_NMDeff_median)[2] <- NMD_method
    samples_NMDeff_median_GTEx_samples <- merge(samples_NMDeff_median_GTEx_samples,samples_NMDeff_median, by.x = "tissue", by.y = "acronyms", all.x = TRUE)
  }
  #colnames(samples_NMDeff_median_GTEx_samples) <- c("tissue","GTEx_NMDeff_ASE","GTEx_NMDeff_endogenous","GTEx_NMDeff_PTCs","GTEx_NMDeff_mean","GTEx_NMDeff_PCA")
  colnames(samples_NMDeff_median_GTEx_samples) <- c("tissue","GTEx_NMDeff_ASE","GTEx_NMDeff_endogenous")

  GTEx_vs_TCGA_median_tissue_correlations <- function(TCGA_GTEx_tissues, NMDeff_median_GTEx, NMDeff_median_TCGA) {
    # Structure pairs of tissues into a list
    TCGA_GTEx_rank_list <- list()
    for (i in 1:nrow(TCGA_GTEx_tissues)) {
      tissue <- TCGA_GTEx_tissues[i,"tissues"]
      TCGA_cancer <- TCGA_GTEx_tissues[i,"TCGA_cancers"]
      GTEx_tissue <- TCGA_GTEx_tissues[i,"GTEx_tissues"]
      TCGA_GTEx_rank_list[[tissue]] <- list(GTEx = NMDeff_median_GTEx[grep(GTEx_tissue, NMDeff_median_GTEx$tissue),], TCGA = NMDeff_median_TCGA[grep(TCGA_cancer, NMDeff_median_TCGA$cancer_type_strat),])
    }
    # Remove TCGA UCEC POLE
    TCGA_GTEx_rank_list$Uterus$TCGA <- TCGA_GTEx_rank_list$Uterus$TCGA[TCGA_GTEx_rank_list$Uterus$TCGA$cancer_type_strat != "TCGA-UCEC_MSS_SBS10ab",]
    # Mean of the tissue values
    TCGA_GTEx_rank_list_means <- lapply(TCGA_GTEx_rank_list, function(list) {
      if (length(unlist(list[["GTEx"]])) > 1) {
        list[["GTEx"]] <- mean(unlist(list[["GTEx"]][,-1]))
      } else {
        list[["GTEx"]] <- as.numeric(unlist(list[["GTEx"]][,-1]))
      }
      if (length(unlist(list[["TCGA"]])) > 1){
        list[["TCGA"]] <- mean(unlist(list[["TCGA"]][,-1]))
      } else {
        list[["TCGA"]] <- as.numeric(unlist(list[["TCGA"]][,-1]))
      }
      list
    })
    GTEx_rank_list_means <- lapply(TCGA_GTEx_rank_list_means, function(list) {list[["GTEx"]]})
    GTEx_ranking_values <- unlist(GTEx_rank_list_means)
    TCGA_rank_list_means <- lapply(TCGA_GTEx_rank_list_means, function(list) {list[["TCGA"]]})
    TCGA_ranking_values <- unlist(TCGA_rank_list_means)
    df <- data.frame(GTEx_ranking_values,TCGA_ranking_values)
    correlation <- cor.test(df$GTEx_ranking_values,df$TCGA_ranking_values)
    return(correlation)
  }

  # 2) Correlations
  n <- length(NMD_genesets)*2
  # corr matrix
  col_names <- c(colnames(samples_NMDeff_median_TCGA_samples)[-1],colnames(samples_NMDeff_median_GTEx_samples)[-1])
  cor_mat <- matrix(1, nrow = n, ncol = n)
  colnames(cor_mat) <- col_names
  rownames(cor_mat) <- col_names
  cor_mat <- data.frame(cor_mat)
  # P-val matrix
  pvas_mat <- matrix(1, nrow = n, ncol = n)
  colnames(pvas_mat) <- col_names
  rownames(pvas_mat) <- col_names
  pvas_mat <- data.frame(pvas_mat)
  # TCGA corr
  corr <- cor(samples_NMDeff_median_TCGA_samples[,-1], use = "pairwise.complete.obs", method = "spearman")
  # TCGA P-values matrix
  pvals <- corr.test(samples_NMDeff_median_TCGA_samples[,-1])$p    # Apply corr.test function
  n2 <- ncol(samples_NMDeff_median_TCGA_samples)-1
  # Add
  cor_mat[1:n2,1:n2] <- corr
  pvas_mat[1:n2,1:n2] <- pvals
  # GTEx corr
  corr <- cor(samples_NMDeff_median_GTEx_samples[,-1], use = "pairwise.complete.obs", method = "spearman")
  # GTEx P-values matrix
  pvals <- corr.test(samples_NMDeff_median_GTEx_samples[,-1])$p    # Apply corr.test function
  # Add
  n3 <- ((n2+1):(n2*2))
  cor_mat[n3,n3] <- corr
  pvas_mat[n3,n3] <- pvals
  # 3) Add GTEx vs TCGA correlations
  # for (TCGA_NMD_geneset in colnames(samples_NMDeff_median_TCGA_samples)[-1]) {
  #   for (GTEx_NMD_geneset in colnames(samples_NMDeff_median_GTEx_samples)[-1]) {
  #     cor_test <- GTEx_vs_TCGA_median_tissue_correlations(TCGA_GTEx_tissues = TCGA_GTEx_tissues, 
  #                                             NMDeff_median_GTEx = samples_NMDeff_median_GTEx_samples[,c("tissue",GTEx_NMD_geneset)],
  #                                             NMDeff_median_TCGA = samples_NMDeff_median_TCGA_samples[,c("cancer_type_strat",TCGA_NMD_geneset)])
  #     cor_mat[TCGA_NMD_geneset,GTEx_NMD_geneset] <- as.numeric(cor_test$estimate)
  #     cor_mat[GTEx_NMD_geneset,TCGA_NMD_geneset] <- as.numeric(cor_test$estimate)
  #     pvas_mat[TCGA_NMD_geneset,GTEx_NMD_geneset] <- as.numeric(cor_test$p.value)
  #     pvas_mat[GTEx_NMD_geneset,TCGA_NMD_geneset] <- as.numeric(cor_test$p.value)
  #   }
  # }

  png(paste0(corr_output_path,"/tissue_ranking/all_methods_correlations_ASE_VAF_",ASE_VAF,".png"), width = 3500, height = 2000, res = 300)
  tryCatch({
    p <- corrplot(as.matrix(cor_mat), type = "upper", method = "pie",
              p.mat = as.matrix(pvas_mat),
              pch = 4,
              title = "Ranking Spearman Correlations of tissue medians between methods",
              tl.cex = 0.7,
              tl.srt= 10,
              pch.cex = 3,
              #sig.level = 0.05
              insig = "label_sig",mar=c(0,0,5,0),tl.offset = 1,
              diag = FALSE,
              order = 'hclust')
    },error = function(e) {
    p <- corrplot(as.matrix(cor_mat), type = "upper", method = "pie",
              p.mat = as.matrix(pvas_mat),
              pch = 4,
              title = "Ranking Spearman Correlations of tissue medians between methods",
              tl.cex = 0.7,
              tl.srt= 10,
              pch.cex = 3,
              #sig.level = 0.05
              insig = "p-value", mar=c(0,0,5,0),tl.offset = 1,
              diag = FALSE,
              order = 'hclust')
      })
  print(p)
  dev.off()

  cols <- grep("ASE|endogenous",colnames(cor_mat))
  rows <- grep("ASE|endogenous",rownames(cor_mat))
  cor_mat_filt <- cor_mat[rows,cols]
  pvas_mat_filt <- pvas_mat[rows,cols]

  png(paste0(corr_output_path,"/tissue_ranking/all_methods_correlations_ASE_VAF_",ASE_VAF,"_subset.png"), width = 3500, height = 2000, res = 300)
  tryCatch({
    p <- corrplot(as.matrix(cor_mat_filt), type = "upper", method = "pie",
              p.mat = as.matrix(pvas_mat_filt),
              pch = 4,
              title = "Ranking Spearman Correlations of tissue medians between methods",
              tl.cex = 0.7,
              tl.srt= 10,
              pch.cex = 3,
              #sig.level = 0.05
              insig = "label_sig",mar=c(0,0,5,0),tl.offset = 1,
              diag = FALSE,
              order = 'hclust')
    },error = function(e) {
    p <- corrplot(as.matrix(cor_mat), type = "upper", method = "pie",
              p.mat = as.matrix(pvas_mat),
              pch = 4,
              title = "Ranking Spearman Correlations of tissue medians between methods",
              tl.cex = 0.7,
              tl.srt= 10,
              pch.cex = 3,
              #sig.level = 0.05
              insig = "p-value", mar=c(0,0,5,0),tl.offset = 1,
              diag = FALSE,
              order = 'hclust')
      })
  print(p)
  dev.off()

  # Do the ranking of tissues by median iNMDeff
  GTEx_iNMDeff_median_ranking <- samples_NMDeff_median_GTEx_samples[order(samples_NMDeff_median_GTEx_samples[,3], decreasing = TRUE),]
  GTEx_iNMDeff_median_ranking[,paste0("GTEx_ETG_iNMDeff_method_ranking")] <- 1:length(GTEx_iNMDeff_median_ranking$tissue)
  GTEx_iNMDeff_median_ranking[,paste0("GTEx_ASE_iNMDeff_method_ranking")] <- rank(-GTEx_iNMDeff_median_ranking$GTEx_NMDeff_ASE)
  # Merge
  GTEx_iNMDeff_median_ranking <- na.omit(GTEx_iNMDeff_median_ranking)
  # Save
  # GTEx_ranking_all[[NMD_method_char]] <<- GTEx_ranking
  # GTEx_ranking[,2] <- NULL

  # Plot
  formula <- y ~ x
  png(paste0(corr_output_path,"/tissue_ranking/GTEx_scatterplot_ASE_VAF_",ASE_VAF,".png"), width = 3500, height = 2000, res = 300)
  p <- ggplot(data = GTEx_iNMDeff_median_ranking, 
            mapping = aes(x = GTEx_ASE_iNMDeff_method_ranking, 
                y = GTEx_ETG_iNMDeff_method_ranking)) +
      geom_point(alpha = 0.5) +
      geom_label_repel(aes(label=tissue, color = "black"), size=3, nudge_y=0.05, max.overlaps = nrow(GTEx_iNMDeff_median_ranking)) +
      geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
      xlab(paste0("ASE iNMDeff tissue ranking")) + 
      ylab(paste0("ETG iNMDeff tissue ranking")) +
      ggtitle("GTEx") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          legend.position='none')
  print(p + stat_cor(p.accuracy = 0.0000000001, r.accuracy = 0.01, size = 5))
  dev.off()

  # save R object
  write.table(GTEx_iNMDeff_median_ranking, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/panel_B.txt", 
              sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
  saveRDS(GTEx_iNMDeff_median_ranking, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/panel_B.RData")

  # Do the ranking of tissues by median iNMDeff
  TCGA_iNMDeff_median_ranking <- samples_NMDeff_median_TCGA_samples[order(samples_NMDeff_median_TCGA_samples[,3], decreasing = TRUE),]
  TCGA_iNMDeff_median_ranking[,paste0("TCGA_ETG_iNMDeff_method_ranking")] <- 1:length(TCGA_iNMDeff_median_ranking$cancer_type_strat)
  TCGA_iNMDeff_median_ranking[,paste0("TCGA_ASE_iNMDeff_method_ranking")] <- rank(-TCGA_iNMDeff_median_ranking$TCGA_NMDeff_ASE)
  # Merge
  TCGA_iNMDeff_median_ranking <- na.omit(TCGA_iNMDeff_median_ranking)
  # Save
  # GTEx_ranking_all[[NMD_method_char]] <<- GTEx_ranking
  # GTEx_ranking[,2] <- NULL

  # Plot
  png(paste0(corr_output_path,"/tissue_ranking/TCGA_scatterplot_ASE_VAF_",ASE_VAF,".png"), width = 3500, height = 2000, res = 300)
  p <- ggplot(data = TCGA_iNMDeff_median_ranking, 
            mapping = aes(x = TCGA_ASE_iNMDeff_method_ranking, 
                y = TCGA_ETG_iNMDeff_method_ranking)) +
      geom_point(alpha = 0.5) +
      geom_label_repel(aes(label=cancer_type_strat, color = "black"), size=3, nudge_y=0.05, max.overlaps = nrow(TCGA_iNMDeff_median_ranking)) +
      geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
      xlab(paste0("ASE iNMDeff tissue ranking")) + 
      ylab(paste0("ETG iNMDeff tissue ranking")) +
      ggtitle("TCGA") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          legend.position='none')
  print(p + stat_cor(p.accuracy = 0.000001, r.accuracy = 0.01, size = 5))
  dev.off()  

  # save R object
  write.table(TCGA_iNMDeff_median_ranking, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/panel_C.txt", 
              sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
  saveRDS(TCGA_iNMDeff_median_ranking, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/panel_C.RData")

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

# # 1.2) PTC NMD efficiencies TCGA
# # PTC
# PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_all_TCGA_confident.txt"
# PTCs_NMD_efficiencies_TCGA <- read.table(file = PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # ASE
# PTCs_ASE_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_ASE_all_TCGA_confident.txt"
# PTCs_ASE_NMD_efficiencies_TCGA <- read.table(file = PTCs_ASE_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # Add sample metadata
# PTCs_ASE_NMD_efficiencies_TCGA <- merge(PTCs_ASE_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)
# PTCs_NMD_efficiencies_TCGA <- merge(PTCs_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)

# 1.3) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)

# 1.4) GTEx - TCGA matched tissues
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_GTEx_match.txt")
TCGA_GTEx_tissues <- read.table(file = output_path, header = TRUE, sep = "\t")

# 2) Correlations between methods
corr_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods"

# 2.1) Sample or PTC-level correlations

# 2.1.1) Sample NMDeff
#### Pancancer ####
# sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA, ASE_VAF = 0.01, dataset = "TCGA")
sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA, ASE_VAF = 0.2, dataset = "TCGA")
# save R object
write.table(sample_NMD_efficiencies_TCGA, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3A.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(sample_NMD_efficiencies_TCGA, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3A.RData")
sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA, ASE_VAF = 0.2, dataset = "TCGA", tissue = NULL)
# sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx, ASE_VAF = 0.01, dataset = "GTEx")
sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx, ASE_VAF = 0.2, dataset = "GTEx")
# save R object
write.table(sample_NMD_efficiencies_GTEx, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3B.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(sample_NMD_efficiencies_GTEx, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3B.RData")
sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx, ASE_VAF = 0.2, dataset = "GTEx", tissue = NULL)

# TCGA Correlations between covariates and iNMDeff
cov_corr_df <- sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA, ASE_VAF = 0.2, dataset = "TCGA", cov_corr = "yes")
# save R object
write.table(cov_corr_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3E.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(cov_corr_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3E.RData")

# GTEx Correlations between covariates and iNMDeff
cov_corr_df <- sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx, ASE_VAF = 0.2, dataset = "GTEx", cov_corr = "yes")
# save R object
write.table(cov_corr_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3F.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(cov_corr_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3F.RData")

#### By cancer ####
# One by one
NMDeff_methods_corr_TCGA <- data.frame()
for (tissue in na.omit(unique(sample_NMD_efficiencies_TCGA$cancer_type_strat))) {
  print(tissue)
  NMDeff_corr <- sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_TCGA, tissue = tissue, ASE_VAF = 0.2, dataset = "TCGA")
  if (length(NMDeff_methods_corr_TCGA) == 0) {
    NMDeff_methods_corr_TCGA <- NMDeff_corr
  } else {
    NMDeff_methods_corr_TCGA <- rbind(NMDeff_methods_corr_TCGA,NMDeff_corr)
  }
}
# How many correlations are positive
sum(NMDeff_methods_corr_TCGA$R > 0)
nrow(NMDeff_methods_corr_TCGA)

NMDeff_methods_corr_GTEx <- data.frame()
for (tissue in unique(sample_NMD_efficiencies_GTEx$acronyms)) {
  print(tissue)
  NMDeff_corr <- sample_NMDeff_corr(sample_NMD_efficiencies = sample_NMD_efficiencies_GTEx, tissue = tissue, ASE_VAF = 0.2, dataset = "GTEx")
    if (length(NMDeff_methods_corr_GTEx) == 0) {
    NMDeff_methods_corr_GTEx <- NMDeff_corr
  } else {
    NMDeff_methods_corr_GTEx <- rbind(NMDeff_methods_corr_GTEx,NMDeff_corr)
  }
}
sum(NMDeff_methods_corr_GTEx$R > 0)
nrow(NMDeff_methods_corr_GTEx)

# Correlations barplot of all tissues/cancers
# Plot


NMDeff_methods_corr_TCGA$database <- "TCGA"
tissue_order <- order(NMDeff_methods_corr_TCGA$R)
NMDeff_methods_corr_TCGA$tissue <- factor(NMDeff_methods_corr_TCGA$tissue, levels = NMDeff_methods_corr_TCGA$tissue[tissue_order])

# save R object
write.table(NMDeff_methods_corr_TCGA, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3C.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(NMDeff_methods_corr_TCGA, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3C.RData")

p1 <- ggplot(data = NMDeff_methods_corr_TCGA, aes(x = tissue, y = R, fill = NULL)) +
    geom_bar(stat = "identity", linewidth = 0.7, fill = "skyblue", width = 0.6) + coord_flip(ylim = c(-0.5,0.5)) + #coord_cartesian(ylim = c(-5,5)) +
    guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
    facet_grid(database ~.) +
    ylab("R") + ggtitle("") + xlab("") +
    geom_errorbar(aes(ymin = R_conf_int_low, ymax = R_conf_int_high), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Pastel1") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=18),# angle = 90, hjust = 1),
        axis.text.y = element_text(color="black", size=14),
        strip.text = element_text(color="black", size = 20),
        legend.position='top', legend.text = element_text(size = 18)) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
                         position = position_dodge(width = 1), size = 15, hjust = -0.1, color = "black", vjust = 0.75)

NMDeff_methods_corr_GTEx$database <- "GTEx"
tissue_order <- order(NMDeff_methods_corr_GTEx$R)
NMDeff_methods_corr_GTEx$tissue <- factor(NMDeff_methods_corr_GTEx$tissue, levels = NMDeff_methods_corr_GTEx$tissue[tissue_order])

# save R object
write.table(NMDeff_methods_corr_GTEx, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3D.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(NMDeff_methods_corr_GTEx, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3D.RData")

p2 <- ggplot(data = NMDeff_methods_corr_GTEx, aes(x = tissue, y = R, fill = NULL)) +
    geom_bar(stat = "identity", linewidth = 0.7, fill = "skyblue", width = 0.6) + coord_flip(ylim = c(-0.5,0.5)) + #coord_cartesian(ylim = c(-5,5)) +
    guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
    facet_grid(database ~.) +
    ylab("R") + ggtitle("") + xlab("") +
    geom_errorbar(aes(ymin = R_conf_int_low, ymax = R_conf_int_high), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Pastel1") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=18),# angle = 90, hjust = 1),
        axis.text.y = element_text(color="black", size=14),
        strip.text = element_text(color="black", size = 20),
        legend.position='top', legend.text = element_text(size = 18)) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
                         position = position_dodge(width = 1), size = 15, hjust = -0.1, color = "black", vjust = 0.75)

p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)

png(paste0(corr_output_path,"/sample_level/correlations_all_tissues.png"), width = 5000, height = 3500, res = 300)
#ggsave(gsub(".png",".pdf",plot_path), p, width = 4500, height = 3500, units = "px")
print(p)
dev.off()


# 2.1.2) PTC-level correlations
# 2.1.2.1) Pancancer
#PTC_NMDeff_corr(PTCs_ASE_NMD_efficiencies_TCGA = PTCs_ASE_NMD_efficiencies_TCGA)

# 2.1.2.2) By cancer
# for (tissue in cancers) {
#   print(tissue)
#   PTC_NMDeff_corr(PTCs_ASE_NMD_efficiencies_TCGA = PTCs_ASE_NMD_efficiencies_TCGA, tissue = tissue)
# }

# 2.2) Tissue-level (median ranking) correlations

# All 5 methods together (either Sample or PTC NMDeff)
# Spearman Rank Correlation of median of tissues
tissue_ranking_corr_matrix(ASE_VAF = "0.01", TCGA_GTEx_tissues = TCGA_GTEx_tissues, 
                    sample_NMD_efficiencies_TCGA = sample_NMD_efficiencies_TCGA,
                     sample_NMD_efficiencies_GTEx = sample_NMD_efficiencies_GTEx)
tissue_ranking_corr_matrix(ASE_VAF = "0.2", TCGA_GTEx_tissues = TCGA_GTEx_tissues, 
                    sample_NMD_efficiencies_TCGA = sample_NMD_efficiencies_TCGA,
                     sample_NMD_efficiencies_GTEx = sample_NMD_efficiencies_GTEx)

# 2.2.1) Scatterplot for each method combination

# 3) Correlation between ASE and ETG for different thresholds of ASE

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

  return(sample_NMDeff)
}

# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)
# sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)

# Add a column to classify tissues as "Brain" or "Non-Brain"
# df <- data.frame(sort(table(sample_NMD_efficiencies_GTEx$acronyms)))
# df$Category <- ifelse(grepl("^BRN", df$Var1), "Brain", "Non-Brain")
# # Calculate mean sample count for each category
# mean_counts <- aggregate(Freq ~ Category, data = df, FUN = mean)

TCGA_all_corr_res_final <- c()
GTEx_all_corr_res_final <- c()

for (VAF in c("0.01","0.2")) {

  TCGA_all_corr_res <- c()
  TCGA_all_brain_res <- c()

  for (threshold in c(1,2,3,4,5,6,7,8,9,10)) {

    # 1) Observed correlations
    # Dataframe
    corr_res_df <- data.frame(corr = NA, pvalue = NA, threshold = NA, sample_size = NA, type = NA)
    TCGA_NMDeff_filt <- sample_NMD_efficiencies_TCGA
    total_samples <- nrow(TCGA_NMDeff_filt)
    # Filter samples for the PTC number threshold
    TCGA_NMDeff_filt <- TCGA_NMDeff_filt[which(TCGA_NMDeff_filt[,paste0("ASE_num_PTCs_",VAF)] >= threshold),]
    # print(cor.test(TCGA_NMDeff_filt$ASE_PTC_NMD_triggering_0.2,TCGA_NMDeff_filt$ASE_num_PTCs_0.2))
    # Check sample size
    sample_size <- sum(!is.na(TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]))
    # Perform Pearson correlation
    pearson_corr_res <- cor.test(TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)], 
          TCGA_NMDeff_filt$endogenous_NMD_Consensus, 
          method = "pearson")
    corr <- pearson_corr_res$estimate
    pvalue <- pearson_corr_res$p.value
    # Save
    corr_res_df$corr <- corr
    corr_res_df$pvalue <- pvalue
    corr_res_df$sample_size <- sample_size
    corr_res_df$threshold <- threshold
    corr_res_df$sample_size_perc <- round(sample_size/total_samples,2)*100
    corr_res_df$type <- "threshold"
    TCGA_all_corr_res <- rbind(TCGA_all_corr_res,corr_res_df)

    # 2) Brain - rest of cancers iNMDeff differences
    
    brain_res_df <- data.frame(mean_brain = NA, mean_rest = NA, p_value = NA, 
        NMD_method = NA, sample_size = NA, threshold = NA, type = NA)

    # Filter for brain cancer types (LGG, GBM) and the rest
    brain_cancer <- TCGA_NMDeff_filt$cancer_type %in% c("LGG", "GBM")
    rest_tissues <- !brain_cancer

    if (sum(brain_cancer) >= 1) {

      # Perform t-tests for the two columns
      t_test_endogenous_NMD <- t.test(
        TCGA_NMDeff_filt$endogenous_NMD_Consensus[brain_cancer],
        TCGA_NMDeff_filt$endogenous_NMD_Consensus[rest_tissues]
      )

      t_test_ASE_PTC <- t.test(
        TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)][brain_cancer],
        TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)][rest_tissues]
      )

      # Display t-test results
      brain_res_df[1,] <- c( round(as.numeric(t_test_endogenous_NMD$estimate[1]),2),
        round(as.numeric(t_test_endogenous_NMD$estimate[2]),2),t_test_endogenous_NMD$p.value,"ETG",
        sample_size, threshold, "threshold")
      brain_res_df[2,] <- c(round(as.numeric(t_test_ASE_PTC$estimate[1]),2),
        round(as.numeric(t_test_ASE_PTC$estimate[2]),2),t_test_ASE_PTC$p.value,"ASE",
        sample_size, threshold, "threshold")
    }
    TCGA_all_brain_res <- rbind(TCGA_all_brain_res,brain_res_df)

    # 3) Correlations with a randomized subset
    subset_all_corr_res <- c()
    for (i in 1:10) {
      # Dataframe
      corr_res_df <- data.frame(corr = NA, pvalue = NA, threshold = NA, sample_size = NA, type = NA)
      TCGA_NMDeff_filt <- sample_NMD_efficiencies_TCGA
      # Remove NAs
      TCGA_NMDeff_filt <- TCGA_NMDeff_filt[!is.na(TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]),]
      # Subset of samples with same sample size as before
      # set.seed(123)
      TCGA_NMDeff_filt <- TCGA_NMDeff_filt[sample(1:nrow(TCGA_NMDeff_filt), sample_size), ]
      # Sample size
      sample_size_subset <- sum(!is.na(TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]))
      # Perform Pearson correlation
      pearson_corr_res <- cor.test(TCGA_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)], 
            TCGA_NMDeff_filt$endogenous_NMD_Consensus, 
            method = "pearson")
      corr <- pearson_corr_res$estimate
      pvalue <- pearson_corr_res$p.value
      # Save
      corr_res_df$corr <- corr
      corr_res_df$pvalue <- pvalue
      corr_res_df$sample_size <- sample_size_subset
      corr_res_df$threshold <- threshold
      corr_res_df$sample_size_perc <- round(sample_size_subset/total_samples,2)*100
      corr_res_df$type <- "subset"
      subset_all_corr_res <- rbind(subset_all_corr_res,corr_res_df)
    }
      # Calculate the mean of corr and pvalue
    subset_all_corr_res <- data.frame(
      corr = mean(subset_all_corr_res$corr),
      pvalue = mean(subset_all_corr_res$pvalue),
      threshold = unique(subset_all_corr_res$threshold), # Assuming threshold is constant
      sample_size = unique(subset_all_corr_res$sample_size), # Assuming sample_size is constant
      sample_size_perc = unique(subset_all_corr_res$sample_size_perc), # Assuming sample_size is constant
      type = unique(subset_all_corr_res$type) # Assuming type is constant
    )
    TCGA_all_corr_res <- rbind(TCGA_all_corr_res,subset_all_corr_res)

  }
  TCGA_all_corr_res$VAF <- VAF
  TCGA_all_corr_res_final <- rbind(TCGA_all_corr_res_final,TCGA_all_corr_res)

  GTEx_all_corr_res <- c()
  GTEx_all_brain_res <- c()

  for (threshold in c(1,2,3,4,5,6,7,8,9,10)) {

    # 1) Observed correlations
    # Dataframe
    corr_res_df <- data.frame(corr = NA, pvalue = NA, threshold = NA, sample_size = NA, type = NA)
    GTEx_NMDeff_filt <- sample_NMD_efficiencies_GTEx
    total_samples <- nrow(GTEx_NMDeff_filt)
    # Filter samples for the PTC number threshold
    GTEx_NMDeff_filt <- GTEx_NMDeff_filt[which(GTEx_NMDeff_filt[,paste0("ASE_num_PTCs_",VAF)] >= threshold),]
    sort(table(GTEx_NMDeff_filt$acronyms))
    # print(cor.test(GTEx_NMDeff_filt$ASE_PTC_NMD_triggering_0.2,GTEx_NMDeff_filt$ASE_num_PTCs_0.2))
    # Check sample size
    sample_size <- sum(!is.na(GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]))
    # Perform Pearson correlation
    pearson_corr_res <- cor.test(GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)], 
          GTEx_NMDeff_filt$endogenous_NMD_Consensus, 
          method = "pearson")
    corr <- pearson_corr_res$estimate
    pvalue <- pearson_corr_res$p.value
    # Save
    corr_res_df$corr <- corr
    corr_res_df$pvalue <- pvalue
    corr_res_df$sample_size <- sample_size
    corr_res_df$threshold <- threshold
    corr_res_df$type <- "threshold"
    corr_res_df$sample_size_perc <- round(sample_size/total_samples,2)*100
    GTEx_all_corr_res <- rbind(GTEx_all_corr_res,corr_res_df)

    # 2) Brain - rest of cancers iNMDeff differences
    
    brain_res_df <- data.frame(mean_brain = NA, mean_rest = NA, p_value = NA, 
          NMD_method = NA, sample_size = NA, threshold = NA, type = NA)

    # Filter for brain cancer types (LGG, GBM) and the rest
    brain_tissues <- grepl("Brain",GTEx_NMDeff_filt$tissue)
    rest_tissues <- !brain_tissues

    # Perform t-tests for the two columns
    t_test_endogenous_NMD <- t.test(
      GTEx_NMDeff_filt$endogenous_NMD_Consensus[brain_tissues],
      GTEx_NMDeff_filt$endogenous_NMD_Consensus[rest_tissues]
    )

    t_test_ASE_PTC <- t.test(
      GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)][brain_tissues],
      GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)][rest_tissues]
    )

    # Display t-test results
    brain_res_df[1,] <- c( round(as.numeric(t_test_endogenous_NMD$estimate[1]),2),
      round(as.numeric(t_test_endogenous_NMD$estimate[2]),2),t_test_endogenous_NMD$p.value,"ETG",
      sample_size, threshold, "threshold")
    brain_res_df[2,] <- c(round(as.numeric(t_test_ASE_PTC$estimate[1]),2),
      round(as.numeric(t_test_ASE_PTC$estimate[2]),2),t_test_ASE_PTC$p.value,"ASE",
      sample_size, threshold, "threshold")

    GTEx_all_brain_res <- rbind(GTEx_all_brain_res,brain_res_df)

    # 3) Correlations with a randomized subset
    subset_all_corr_res <- c()
    for (i in 1:10) {
      # Dataframe
      corr_res_df <- data.frame(corr = NA, pvalue = NA, threshold = NA, sample_size = NA, type = NA)
      GTEx_NMDeff_filt <- sample_NMD_efficiencies_GTEx
      # Remove NAs
      GTEx_NMDeff_filt <- GTEx_NMDeff_filt[!is.na(GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]),]
      # Subset of samples with same sample size as before
      # set.seed(123)
      GTEx_NMDeff_filt <- GTEx_NMDeff_filt[sample(1:nrow(GTEx_NMDeff_filt), sample_size), ]
      # Sample size
      sample_size_subset <- sum(!is.na(GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)]))
      # Perform Pearson correlation
      pearson_corr_res <- cor.test(GTEx_NMDeff_filt[,paste0("ASE_PTC_NMD_triggering_",VAF)], 
            GTEx_NMDeff_filt$endogenous_NMD_Consensus, 
            method = "pearson")
      corr <- pearson_corr_res$estimate
      pvalue <- pearson_corr_res$p.value
      # Save
      corr_res_df$corr <- corr
      corr_res_df$pvalue <- pvalue
      corr_res_df$sample_size <- sample_size_subset
      corr_res_df$sample_size_perc <- round(sample_size_subset/total_samples,2)*100
      corr_res_df$threshold <- threshold
      corr_res_df$type <- "subset"
      subset_all_corr_res <- rbind(subset_all_corr_res,corr_res_df)
    }
      # Calculate the mean of corr and pvalue
    subset_all_corr_res <- data.frame(
      corr = mean(subset_all_corr_res$corr),
      pvalue = mean(subset_all_corr_res$pvalue),
      threshold = unique(subset_all_corr_res$threshold), # Assuming threshold is constant
      sample_size = unique(subset_all_corr_res$sample_size), # Assuming sample_size is constant
      sample_size_perc = unique(subset_all_corr_res$sample_size_perc), # Assuming sample_size is constant
      type = unique(subset_all_corr_res$type) # Assuming type is constant
    )
    GTEx_all_corr_res <- rbind(GTEx_all_corr_res,subset_all_corr_res)
  }

  GTEx_all_corr_res$VAF <- VAF
  GTEx_all_corr_res_final <- rbind(GTEx_all_corr_res_final,GTEx_all_corr_res)

}

GTEx_all_corr_res_final$dataset <- "GTEx"
TCGA_all_corr_res_final$dataset <- "TCGA"

# Merge
all_corr_res <- rbind(GTEx_all_corr_res_final,TCGA_all_corr_res_final)

# Add significance levels for p-values
all_corr_res$significance <- cut(all_corr_res$pvalue,
                                 breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                 labels = c("***", "**", "*", ".", "ns"))

# Add a new column to highlight threshold == 3
all_corr_res$highlight <- ifelse(all_corr_res$threshold == 3, "highlight", "normal")
all_corr_res$VAF <- paste0("AF <= ",all_corr_res$VAF)
all_corr_res$VAF <- factor(all_corr_res$VAF, levels = c("AF <= 0.2","AF <= 0.01"))

# Plot
plot <- ggplot(all_corr_res, aes(x = factor(threshold), y = corr, fill = factor(type), 
        color = factor(highlight), linetype = factor(highlight)) ) +
  geom_bar(stat = "identity", position = position_dodge(), size = 0.5) +
  geom_text(aes(label = sample_size_perc), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
  geom_text(aes(label = significance), position = position_dodge(width = 0.9), vjust = -2, size = 2) +
  labs(x = "Number of PTCs per sample (ASE iNMDeff method)", y = "ASE vs ETG iNMDeff Pearson Correlation (R)", fill = "Samples used") +
  facet_wrap(dataset ~ VAF, scale = "free") +
  scale_fill_manual(
    name = "Samples used",
    labels = c("Subset", "Threshold"),
    values = c(
    "subset" = "#91afe2",
    "threshold" = "#3cdaa0"
  )) +
  scale_color_manual(
    values = c("normal" = "black", "highlight" = "red")
  ) +
  scale_linetype_manual(
    values = c("normal" = "solid", "highlight" = "dashed")
  ) +
  guides(color = "none", linetype = "none") +  # Remove the fill legend
  theme_classic() +
  theme(
      axis.text.x = element_text(),
      legend.position = "top",
      strip.background = element_blank(),  # Remove panel background
      strip.text = element_text(size = 12) # Keep titles
        )

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/sample_level/corr_by_PTC_ASE_thresholds_and_sample_size.png")
ggsave(final_figure_path, plot, width = 180, height = 150, units = "mm") 

# save R object
write.table(all_corr_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5A.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(all_corr_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5A.RData")

# 4)  BRAIN DIFFERENCES BY THRESHOLD
GTEx_all_brain_res$dataset <- "GTEx"
TCGA_all_brain_res$dataset <- "TCGA"
all_brain_res <- rbind(GTEx_all_brain_res,TCGA_all_brain_res)

all_brain_res$mean_brain <- as.numeric(all_brain_res$mean_brain)
all_brain_res$mean_rest <- as.numeric(all_brain_res$mean_rest)
all_brain_res$threshold <- factor(all_brain_res$threshold, levels = 1:10)

# Create a line plot for mean_brain and mean_rest across thresholds
plot <- ggplot(all_brain_res, aes(x = threshold)) +
  geom_line(aes(y = mean_brain, color = "Mean Brain", linetype = NMD_method), size = 1) +
  geom_line(aes(y = mean_rest, color = "Mean Rest", linetype = NMD_method), size = 1) +
  geom_point(aes(y = mean_brain, color = "Mean Brain", shape = NMD_method), size = 3) +
  geom_point(aes(y = mean_rest, color = "Mean Rest", shape = NMD_method), size = 3) +
  coord_cartesian(ylim = c(-1.5,0.5)) +
  facet_wrap(dataset ~.)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) + # Add horizontal dashed line
  scale_color_manual(
    name = "Metric",
    values = c("Mean Brain" = "blue", "Mean Rest" = "red")
  ) +
  labs(
    x = "Threshold --> # Number PTCs per sample for ASE",
    y = "Mean tissue iNMDeff",
    title = "Mean Brain vs. Rest Across Thresholds by NMD Method"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 14)
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/Brain_t_test_diff_using_PTC_ASE_thresholds_VAF_0.01.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

# 4) Correlation between ASE and ETG for different thresholds of ASE, done by each tissue

library(purrr)
library(glue)

# Function to process data for a single threshold
corr_by_tissue <- function(threshold, samples_iNMDeff, dataset, VAF) {
  # Filter the data based on the current threshold
  samples_iNMDeff_tmp <- samples_iNMDeff %>%
    filter(.data[[glue("ASE_num_PTCs_{VAF}")]] >= threshold)
  if (nrow(samples_iNMDeff_tmp)==0){next}
  if (dataset == "TCGA") {
    grouping_var <- "cancer_type_strat"
  } else if (dataset == "GTex") {
    grouping_var <- "acronyms"
  }
  
  # Perform group-wise correlation test
  cor_results <- samples_iNMDeff_tmp %>%
    group_by(across(all_of(grouping_var))) %>%  # Dynamically group by the specified variable
    filter(n() > 2) %>%  # Skip groups with 2 or fewer rows
    # filter(sum(!is.na(ASE_PTC_NMD_triggering_0.01) & !is.na(endogenous_NMD_Consensus)) > 1) %>%  # Ensure enough valid pairs
    summarise(
      cor_test = list(
        tryCatch(
          cor.test(
            .data[[glue("ASE_PTC_NMD_triggering_{VAF}")]],
            endogenous_NMD_Consensus, 
            method = "pearson", 
            use = "pairwise.complete.obs"
          ), 
          error = function(e) NULL  # Handle groups where correlation can't be computed
        )
      ), .groups = "drop"
    ) %>%
    mutate(
      estimate = map_dbl(cor_test, ~ if (!is.null(.x)) .x$estimate else NA_real_),  # Extract correlation coefficient
      p_value = map_dbl(cor_test, ~ if (!is.null(.x)) .x$p.value else NA_real_)    # Extract p-value
    )
  
  # Calculate summary statistics for the current threshold
  cor_results %>%
    summarise(
      threshold = threshold,  # Add the threshold for reference
      mean_correlation = round(mean(estimate, na.rm = TRUE),2),
      median_correlation = round(median(estimate, na.rm = TRUE),2),
      positive_correlation_count = sum(estimate > 0, na.rm = TRUE),
      Q1 = round(quantile(estimate, 0.25, na.rm = TRUE), 2),
      Q3 = round(quantile(estimate, 0.75, na.rm = TRUE), 2),
      SD = round(sd(estimate, na.rm = TRUE), 2),
      total_tissues = n(),
      positive_correlation_percentage = round((positive_correlation_count / total_tissues) * 100,2)
    )
}

# Apply the function to all thresholds (1:10)

all_final_results <- c()

for (VAF in c("0.01","0.2")) {

  thresholds <- 1:10
  GTEx_final_results <- map_dfr(
    thresholds,
    ~ corr_by_tissue(
      threshold = .x,
      dataset = "GTex",
      VAF = VAF,
      samples_iNMDeff = sample_NMD_efficiencies_GTEx
    )
  )
  GTEx_final_results$dataset <- "GTex"

  TCGA_final_results <- map_dfr(
    thresholds,
    ~ corr_by_tissue(
      threshold = .x,
      dataset = "TCGA",
      VAF = VAF,
      samples_iNMDeff = sample_NMD_efficiencies_TCGA
    )
  )

  TCGA_final_results$dataset <- "TCGA"

  final_results <- rbind(TCGA_final_results,GTEx_final_results)
  final_results$VAF <- VAF
  all_final_results <- rbind(all_final_results, final_results)
}

# View results
print(data.frame(all_final_results))


# Create the plot
plot <- ggplot(all_final_results, aes(x = threshold)) +
  geom_line(aes(y = mean_correlation, color = "Mean Correlation"), size = 1) +
  geom_line(aes(y = median_correlation, color = "Median Correlation"), size = 1, linetype = "dashed") +
  # geom_line(aes(y = positive_correlation_percentage / 100, color = "% of tissues with positive correlation"), size = 1, linetype = "dotted") +
  geom_bar(aes(y = positive_correlation_percentage / 100, fill = "% of tissues with positive correlation"), 
           stat = "identity", position = "dodge", alpha = 0.6) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_x_continuous(
    breaks = 1:10  # Explicitly set X-axis breaks from 1 to 10
  ) +
  scale_y_continuous(
    breaks = c(0.25,0.5,0.75,1)  # Explicitly set X-axis breaks from 1 to 10
  ) +  scale_color_manual(
    values = c("Mean Correlation" = "blue", 
               "Median Correlation" = "green", 
               "% of tissues with positive correlation" = "red")
  ) +
  labs(
    title = "AF <= 0.2",
    fill = "",
    x = "Number of PTCs per sample (ASE iNMDeff method)",
    y = "ASE vs ETG iNMDeff Pearson Correlation (R)",
    color = "Metric"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/sample_level/by_tissue_corr_by_PTC_ASE_thresholds_and_sample_size.png"
ggsave(final_figure_path, plot, width = 185, height = 120, units = "mm") 

# save R object
write.table(all_final_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5B.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
saveRDS(all_final_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5B.RData")

# 5) Histogram of number of targets (ETG) and PTCs (ASE)
#ETG
plot <- ggplot(sample_NMD_efficiencies_TCGA, aes(x = endogenous_num_NMD_targets)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
  labs(
    title = "ETG iNMDeff method in TCGA",
    x = "Number of NMD Targets",
    y = "Frequency"
  ) +
  theme_classic() + xlim(0,130)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/TCGA_ETG_histogram_num_NMD_targets.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

plot <- ggplot(sample_NMD_efficiencies_GTEx, aes(x = endogenous_num_NMD_targets)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(
    title = "ETG iNMDeff method in GTEx",
    x = "Number of NMD Targets",
    y = "Frequency"
  ) +
  theme_classic()  + xlim(0,130)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/GTEx_ETG_histogram_num_NMD_targets.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

#ASE 
VAF <- "0.01"
plot <- ggplot(sample_NMD_efficiencies_TCGA, aes(x = ASE_num_PTCs_0.01)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(
    title = "ASE iNMDeff method in TCGA",
    x = "Number of PTCs per sample",
    y = "Frequency"
  ) +
  theme_classic()

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/TCGA_ASE_histogram_num_PTCs_VAF_",VAF,".png")
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

VAF <- "0.01"

plot <- ggplot(sample_NMD_efficiencies_GTEx, aes(x = ASE_num_PTCs_0.01)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(
    title = "ASE iNMDeff method in GTEx",
    x = "Number of PTCs per sample",
    y = "Frequency"
  ) +
  theme_classic()

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/GTEx_ASE_histogram_num_PTCs_VAF_",VAF,".png")
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

# 6) Correlation between ASE and ETG for different thresholds of AF (ASE)

