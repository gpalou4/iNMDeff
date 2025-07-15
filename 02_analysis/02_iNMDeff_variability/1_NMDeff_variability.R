rm(list=ls())

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

# Libraries
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("stringr")
library("ggpubr")

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################


pantissue_boxplot_NMDeff <- function(dataset, method, NMD_genesets, VAF = NULL) {

  if (dataset == "TCGA") {
    NMD_efficiencies_orig <- sample_NMD_efficiencies_TCGA
    if (method == "endogenous") {
      plot_path <- paste0(TCGA_boxplot_endogenous_pancancer_path,paths[paths$folder_or_object=="NMDeff_endogenous_pancancer_boxplot","path_or_filename"])
    } else if (method == "ASE") {
      plot_path <- paste0(TCGA_boxplot_ASE_pancancer_path,paths[paths$folder_or_object=="NMDeff_ASE_pancancer_boxplot","path_or_filename"])
    } else if (method == "PTCs") {
      plot_path <- paste0(TCGA_boxplot_PTCs_pancancer_path,paths[paths$folder_or_object=="NMDeff_PTCs_pancancer_boxplot","path_or_filename"])
    }
  } else if (dataset == "GTEx") {
    NMD_efficiencies_orig <- sample_NMD_efficiencies_GTEx
    if (method == "endogenous") {
      plot_path <- paste0(GTEx_boxplot_endogenous_pantissue_path,paths[paths$folder_or_object=="NMDeff_endogenous_pantissue_boxplot","path_or_filename"])
    } else if (method == "ASE") {
      plot_path <- paste0(GTEx_boxplot_ASE_pantissue_path,paths[paths$folder_or_object=="NMDeff_ASE_pantissue_boxplot","path_or_filename"])
    } else if (method == "PTCs") {
      plot_path <- paste0(GTEx_boxplot_PTCs_pantissue_path,paths[paths$folder_or_object=="NMDeff_PTCs_pantissue_boxplot","path_or_filename"])
    }
  }
  sample_NMD_efficiencies_filt <- NMD_efficiencies_orig[,NMD_genesets]
  nb_coeff_res_filt <- stack(sample_NMD_efficiencies_filt)
  nb_coeff_res_filt$NMD_geneset <- NA
  nb_coeff_res_filt[grep("NMD", nb_coeff_res_filt$ind),"NMD_geneset"] <- "NMD"
  nb_coeff_res_filt[grep("RandomGenes", nb_coeff_res_filt$ind),"NMD_geneset"] <- "RandomGenes"
  # Y label text 
  if (method == "endogenous") {
    nb_coeff_res_filt$ind <- str_wrap(nb_coeff_res_filt$ind, width = 15)  # Adjust the width as necessary
  } else {
    nb_coeff_res_filt$ind <- gsub(" 0.2| 0.01","",nb_coeff_res_filt$ind)
  }

  if ( method %in% c("ASE") ) {
    nb_coeff_res_filt$NMD_geneset <- nb_coeff_res_filt$ind
    legend <- "none"
    factor_levels <- paste0(c("PTC NMD-triggering","PTC NMD-evading","Synonymous"))
    nb_coeff_res_filt$ind <- factor(nb_coeff_res_filt$ind, levels = factor_levels)
    plot_path <- gsub(paste0("NMDeff_",method),paste0("NMDeff_",method,"_",VAF),plot_path)
    width <- 4500
    height <- 4000
    ylim <- c(-3,3)
    pval_ylabel <- c(1.25,1.75,2.25)
    boxplot_width <- 0.15
    annotate_hjust <- 1
    combinations <- combn(names(table(nb_coeff_res_filt$ind)), 2, simplify = FALSE)
  } else if (method == "endogenous") {
    legend <- "none"
    factor_levels <- c("RandomGenes","NMD")
    #factor_levels_2 <- c("RandomGenes without NMD features","RandomGenes with NMD features","NMD All","NMD Consensus")
    width <- 7000
    height <- 6000
    ylim <- c(-2,7)
    annotate_hjust <- -1
    pval_ylabel <- c(5,6)
    boxplot_width <- 0.3
    combinations <- list(c("RandomGenes\nwith NMD\nfeatures","NMD Consensus"),c("RandomGenes\nwithout NMD\nfeatures","NMD Consensus"))
  }
  # NAs
  tmp_df <- nb_coeff_res_filt[!is.na(nb_coeff_res_filt$values),]
  
  # All NMD genesets
  png(plot_path, width = width-500, height = height-500, res = 300)
  p <- ggplot(data = nb_coeff_res_filt, aes(x = ind, y = values, fill = factor(NMD_geneset, levels = factor_levels))) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = ylim) +
    geom_boxplot(width=boxplot_width, color="black", alpha=0.2) +
    geom_point(alpha=0.25) + xlab("") +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    ylab("NMD efficiency") +
    theme_bw(base_size = 15) + #guides(color = guide_legend(title = "NMD genesets")) +
    theme(plot.title = element_text(hjust = 0.5, size = 35),
          axis.title.x = element_text(color="black", size=40),
          axis.title.y = element_text(color="black", size=35),
          axis.text.x = element_text(color="black", size=40),
          axis.text.y = element_text(color="black", size=35),
          legend.position = legend, axis.text=element_text(size=18)) +
    scale_fill_brewer(palette = "Dark2") +
    stat_compare_means(comparisons = combinations, size = 10,
                  label.y = pval_ylabel,
                  label = "p.format", method = "wilcox.test", hide.ns = TRUE) +  
    annotate("text",
            x = 1:length(table(tmp_df$ind)),
            y = aggregate( values ~ ind, tmp_df, median)[ , 2],
            label = table(tmp_df$ind),
            col = "black",
            hjust = annotate_hjust,
            vjust = -1,
            size = 10)
  ggsave(gsub(".png",".pdf",plot_path), p, width = width, height = height, units = "px")
  print(p)
  dev.off()

  if (method == "endogenous") {
    # Subset of NMD genesets
    #NMD_genesets <- c("RandomGenes\nwith NMD\nfeatures","RandomGenes\nwithout NMD\nfeatures","NMD Consensus", "NMD All")
    #nb_coeff_res_filt2 <- nb_coeff_res_filt[nb_coeff_res_filt$ind %in% NMD_genesets,]
    nb_coeff_res_filt2 <- nb_coeff_res_filt
    # NAs
    tmp_df <- nb_coeff_res_filt2[!is.na(nb_coeff_res_filt2$values),]

    plot_path_2 <- gsub("pancancer_boxplot","pancancer_boxplot_subset",plot_path)
    png(plot_path_2, width = (width-3000), height = (height-3000), res = 300)
    p <- ggplot(data = nb_coeff_res_filt2, aes(x = factor(ind), y = values, fill = factor(NMD_geneset, levels = factor_levels))) +
      geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = ylim) +
      geom_boxplot(width=0.3, color="black", alpha=0.2) +
      geom_point(alpha=0.25) + xlab("") +
      geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
      ylab("NMD efficiency") +
      theme_bw(base_size = 15) + #guides(color = guide_legend(title = "NMD genesets")) +
      theme(plot.title = element_text(hjust = 0.5, size = 35),
            axis.title.x = element_text(color="black", size=35),
            axis.title.y = element_text(color="black", size=30),
            axis.text.x = element_text(color="black", size=35),
            axis.text.y = element_text(color="black", size=30),
            legend.position = legend, axis.text=element_text(size=18)) +
      scale_fill_brewer(palette = "Dark2") +
      stat_compare_means(comparisons = combinations, size = 10,
                    label.y = c(5,6),
                    label = "p.format", method = "wilcox.test", hide.ns = TRUE) +  
      annotate("text",
              x = 1:length(table(tmp_df$ind)),
              y = aggregate( values ~ ind, tmp_df, median)[ , 2],
              label = table(tmp_df$ind),
              col = "black",
              hjust = - 1,
              vjust = - 1,
              size = 8)
    #ggsave(gsub(".png",".pdf",plot_path), p, width = width, height = height, units = "px")
    print(p)
    dev.off()
    return(nb_coeff_res_filt2)
  }
  return(nb_coeff_res_filt)
}

randomization_test <- function(NMDeff_rand, NMD_geneset) {
  # Reorder NB coeffs randomly (artificial random tissues)
  NMDeff_rand$values <- sample(NMDeff_rand$values)
  # Calculate medians per tissue
  NMDeff_rand_medians <- aggregate(.~ind,NMDeff_rand[,c(1:2)],median,na.rm = TRUE)
  # Calculate SD per tissue
  NMDeff_rand_sd <- aggregate(.~ind,NMDeff_rand[,c(1:2)],sd)
  # Difference between each tissue median and median across tissues
  rand_scores <- as.numeric(abs(NMDeff_rand_medians$values) - abs(median(NMDeff_rand$values,na.rm = TRUE)))
  rand_scores_sorted <- sort(abs(rand_scores),decreasing = TRUE)
  # Randomized SD of medians
  rand_sd_of_medians <- sd(NMDeff_rand_medians$values)
  return(list(rand_scores_sorted = rand_scores_sorted, rand_sd_of_medians = rand_sd_of_medians, rand_tissue_sd = NMDeff_rand_sd, 
              rand_tissue_median = NMDeff_rand_medians, NMDeff_rand = NMDeff_rand))
}

boxplots_by_NMDgeneset <- function(NMD_geneset, dataset, method, NMD_efficiencies_df, n_rand) {
  
  print(NMD_geneset)
  NMD_geneset_name <- gsub("\\.","_", NMD_geneset)
  
  ylim_SD <- c(-0.4,0.4)
  if ( (dataset == "TCGA") && (method == "endogenous") )  {
    plot_path_tmp <- paste0(TCGA_boxplots_endogenous_NMD_genesets_path)
  } else if ( (dataset == "TCGA") && (method == "PTCs") ) {
    plot_path_tmp <- paste0(TCGA_boxplots_PTCs_NMD_genesets_path)
  } else if ( (dataset == "TCGA") && (method == "ASE") ) {
    plot_path_tmp <- paste0(TCGA_boxplots_ASE_NMD_genesets_path)
  } else if ( (dataset == "TCGA") && (method == "consensus") ) {
    plot_path_tmp <- paste0(TCGA_boxplots_consensus_NMD_genesets_path)
  } else if ( (dataset == "GTEx") && (method == "endogenous") ) {
    plot_path_tmp <- paste0(GTEx_boxplots_endogenous_NMD_genesets_path)
  } else if ( (dataset == "GTEx") && (method == "ASE") ) {
    plot_path_tmp <- paste0(GTEx_boxplots_ASE_NMD_genesets_path)
  } else if ( (dataset == "GTEx") && (method == "PTCs") ) {
    plot_path_tmp <- paste0(GTEx_boxplots_PTCs_NMD_genesets_path)
  } else if ( (dataset == "GTEx") && (method == "consensus") ) {
    plot_path_tmp <- paste0(GTEx_boxplots_consensus_NMD_genesets_path)
  }
  if (dataset == "GTEx") {
    tissue <- "acronyms"
     x.size <- 15
  } else if (dataset == "TCGA") {
    tissue <- "cancer_type_strat"
    x.size <- 18
  }

  NMD_efficiencies_df_filt <- NMD_efficiencies_df[,c(tissue,NMD_geneset)]
  # Add sample size
  NAs <- which(is.na(NMD_efficiencies_df_filt[,NMD_geneset]))
  if (length(NAs) != 0) {
    tmp <- NMD_efficiencies_df_filt[-NAs,]
  } else {
    tmp <- NMD_efficiencies_df_filt
  }
  sample_size_df <- data.frame(table(tmp[,tissue]))
  colnames(sample_size_df)[2] <- "sample_size"
  NMD_efficiencies_df_original <- NMD_efficiencies_df_filt
  # Scale
  NMD_efficiencies_df_filt[,NMD_geneset] <- scale(NMD_efficiencies_df_filt[,NMD_geneset])

  add_info <- function(NMD_eff_df, medians = NULL) {
    NMD_eff_df$id <- 1:nrow(NMD_eff_df)
    NMD_eff_df <- merge(NMD_eff_df,sample_size_df, by.x = tissue, by.y = "Var1", all.x = TRUE, sort = FALSE)
    NMD_eff_df$tissue_ss <- paste0(NMD_eff_df[,tissue]," (",NMD_eff_df$sample_size,")")
    NMD_eff_df <- NMD_eff_df[order(NMD_eff_df$id),]  
    # Calculate percentile NMDeff
    NMDeff_medians_df <- aggregate(NMD_eff_df[,2] ~ eval(parse(text = tissue )),NMD_eff_df, median)
    colnames(NMDeff_medians_df)[1] <- tissue
    NMDeff_medians_df <- NMDeff_medians_df[order(NMDeff_medians_df[,2]),]
    NMDeff_medians <- NMDeff_medians_df[,2]
    names(NMDeff_medians) <- NMDeff_medians_df[,tissue]  
    if (!is.null(medians)) {
      return(list(NMDeff_medians = NMDeff_medians,NMDeff_medians_df = NMDeff_medians_df))
    }
    # Order labels by NMD eff median
    NMD_eff_df[,tissue] <- factor(NMD_eff_df[,tissue], levels= names(NMDeff_medians) )
    sample_size_df_order <- sample_size_df[match(NMDeff_medians_df[,tissue],sample_size_df$Var1),]
    NMD_eff_df$tissue_ss <- factor(NMD_eff_df$tissue_ss, levels=c(paste0(NMDeff_medians_df[,tissue], " (",sample_size_df_order$sample_size,")")) )
    NMD_eff_df$NMD_geneset <- NMD_geneset
    NMD_eff_df <- na.omit(NMD_eff_df)
    colnames(NMD_eff_df) <- c("ind","values","id","sample_size","tissue_ss",NMD_geneset)
    return(NMD_eff_df)
  }
  NMDeff_medians <- add_info(NMD_eff_df = NMD_efficiencies_df_filt, medians = "yes")
  NMD_efficiencies_df_filt <- add_info(NMD_eff_df = NMD_efficiencies_df_filt)
  NMD_efficiencies_df_original <- add_info(NMD_eff_df = NMD_efficiencies_df_original)
  
  plot_name1 <- ""

  ##  Randomization test ##
  NMDeff_rand <- NMD_efficiencies_df_filt
  # Calculate observed score
  # Difference between each tissue median and median across tissues
  obs_scores_tissues <- abs(NMDeff_medians$NMDeff_medians) - abs(median(NMDeff_medians$NMDeff_medians_df[,2]))
  obs_scores <- as.numeric(obs_scores_tissues)
  obs_scores_sorted <- sort(abs(obs_scores),decreasing = TRUE)
  # Observed SD of medians
  obs_sd_of_medians <- sd(NMDeff_medians$NMDeff_medians)
  # Observed tissue SD
  obs_tissue_sd <- aggregate(.~ind,NMDeff_rand[,c(1:2)],sd,na.rm = TRUE)
  colnames(obs_tissue_sd)[2] <- "values"

  # Peform N randomizations and obtain a distribution of random scores (our null distribution)
  randomizations <- n_rand
  rand_scores_list <- list()
  
  for (iteration in 1:randomizations) {
    rand_scores <- randomization_test(NMDeff_rand = NMDeff_rand)
    rand_scores_list[[iteration]] <- rand_scores
  }
  # Non-scaled
  rand_scores_list_non_scaled <- list()
  for (iteration in 1:randomizations) {
    rand_scores <- randomization_test(NMDeff_rand = NMD_efficiencies_df_original)
    rand_scores_list_non_scaled[[iteration]] <- rand_scores
  }

  # 1) Median difference between best X tissues
  # Calculate p-value for randomization test
  # How many times the random score is higher than the observed score divided by number of randomizations performed (+1 is for avoiding pval = 0)
  obs_scores_tissues <- abs(obs_scores_tissues)
  for (tissue_iteration in 1:length(unique(NMD_efficiencies_df_filt$ind))) {
    # Random scores
    rand_scores_iteration_list  <- lapply(rand_scores_list, function(X){X[["rand_scores_sorted"]][tissue_iteration]})
    rand_scores_iteration <- unlist(rand_scores_iteration_list)
    obs_score_value <- obs_scores_sorted[tissue_iteration]
    rand_pval <- (sum(rand_scores_iteration > obs_score_value)+1)/randomizations
    # Tissue type
    obs_score_tissue <- names(obs_scores_tissues[round(obs_scores_tissues,8) == round(obs_score_value,8)])
    randomization_list[[dataset]][[method]]$randomization_pvals[paste0(obs_score_tissue,"_",tissue_iteration),NMD_geneset] <<- rand_pval
  }
  # 2) SD of medians across tissues (Inter-tissue randomization test)
  # Randomization SD of medians
  rand_sd_of_medians_iteration_list  <- lapply(rand_scores_list, function(X){X[["rand_sd_of_medians"]]})
  rand_sd_of_medians <- unlist(rand_sd_of_medians_iteration_list)
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["randomization_sd_of_medians"]][NMD_geneset] <<- rand_sd_of_medians
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["obs_sd_of_medians"]][NMD_geneset] <<- obs_sd_of_medians
  # Take the median, 5th and 95th percentiles of the SD of the medians of iNMDeff
  rand_median_sd_of_medians <- median(rand_sd_of_medians)
  rand_5perc_sd_of_medians <- as.numeric(quantile(rand_sd_of_medians, probs = 0.05, na.rm = TRUE))
  rand_95perc_sd_of_medians <- as.numeric(quantile(rand_sd_of_medians, probs = 0.95, na.rm = TRUE))
  # Difference between Observed SD of medians and (5th, 50th, 95th percentiles) Randomized SD of medians
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["sd_of_medians_diff"]][NMD_geneset,"median_SD_of_medians_diff"] <<- obs_sd_of_medians - rand_median_sd_of_medians
  print(obs_sd_of_medians - rand_median_sd_of_medians)
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["sd_of_medians_diff"]][NMD_geneset,"perc5_SD_of_medians_diff"] <<- obs_sd_of_medians - rand_5perc_sd_of_medians
  print(obs_sd_of_medians - rand_5perc_sd_of_medians)
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["sd_of_medians_diff"]][NMD_geneset,"perc95_SD_of_medians_diff"] <<- obs_sd_of_medians - rand_95perc_sd_of_medians
  print(obs_sd_of_medians - rand_95perc_sd_of_medians)
  # P-value
  rand_pval <- (sum(rand_sd_of_medians > obs_sd_of_medians)+1)/randomizations
  randomization_list[[dataset]][[method]]$randomization_sd_of_medians[["sd_of_medians_diff"]][NMD_geneset,"p_value"] <<- rand_pval
  
  # 3) Tissue SD and median
  # 3.1) Observed tissue SD vs Randomized tissue SD
  # Random tissue SD
  rand_tissue_sd_iteration_list  <- lapply(rand_scores_list_non_scaled, function(X){X[["rand_tissue_sd"]]})
  # Median tissue SD

  tissues_rand_SD <- do.call(rbind,rand_tissue_sd_iteration_list)
  tissues_obs_SD <- NMD_efficiencies_df_original %>%
                    group_by(ind) %>%
                    summarise(SD = sd(values)) %>%
                    arrange(desc(SD))  
  tissues_SD_all <- merge(tissues_rand_SD, tissues_obs_SD, all.x = TRUE)
  colnames(tissues_SD_all) <- c("tissues","rand_SD","obs_SD")
  tissues_SD_all$diff_SD <- tissues_SD_all$obs_SD - tissues_SD_all$rand_SD
  # P-value
  tissues_SD_all <- tissues_SD_all %>%
                    group_by(tissues) %>%
                    mutate(p_value_above_SD = (sum(rand_SD > obs_SD)+1)/randomizations) %>%
                    mutate(p_value_below_SD = (sum(rand_SD < obs_SD)+1)/randomizations) %>%
                    mutate(color = ifelse(p_value_below_SD < 0.05 | p_value_above_SD < 0.05,TRUE,FALSE))            
  # Order
  tissues_SD_all_medians <- tissues_SD_all %>% 
      group_by(tissues) %>%
      summarise(median = median(diff_SD)) %>%
      arrange(desc(median))
  tissues_SD_all$tissues <- factor(tissues_SD_all$tissues, levels = rev(as.character(tissues_SD_all_medians$tissues)))
    
  # 3.2) Observed tissue median vs Randomized tissue median
  # Random tissue NMDeff
  rand_tissue_median_iteration_list  <- lapply(rand_scores_list_non_scaled, function(X){X[["rand_tissue_median"]]})
  tissues_rand_NMDeff_median <- do.call(rbind,rand_tissue_median_iteration_list)
  tissues_obs_NMDeff_median <- NMD_efficiencies_df_original %>%
                    group_by(ind) %>%
                    summarise(median = median(values)) %>%
                    arrange(desc(median))  
  tissues_NMDeff_median_all <- merge(tissues_rand_NMDeff_median, tissues_obs_NMDeff_median, all.x = TRUE)
  colnames(tissues_NMDeff_median_all) <- c("tissues","rand_NMDeff_median","obs_NMDeff_median")
  tissues_NMDeff_median_all$diff_NMDeff_median <- tissues_NMDeff_median_all$obs_NMDeff_median - tissues_NMDeff_median_all$rand_NMDeff_median
  # P-value
  tissues_NMDeff_median_all <- tissues_NMDeff_median_all %>%
                    group_by(tissues) %>%
                    mutate(p_value_above_NMDeff_median = (sum(rand_NMDeff_median > obs_NMDeff_median)+1)/randomizations) %>%
                    mutate(p_value_below_NMDeff_median = (sum(rand_NMDeff_median < obs_NMDeff_median)+1)/randomizations) %>%
                    mutate(color = ifelse(p_value_below_NMDeff_median < 0.05 | p_value_above_NMDeff_median < 0.05,TRUE,FALSE))            
  # Order
  tissues_NMDeff_median_all_medians <- tissues_NMDeff_median_all %>% 
      group_by(tissues) %>%
      summarise(median = median(diff_NMDeff_median)) %>%
      arrange(desc(median))
  tissues_NMDeff_median_all$tissues <- factor(tissues_NMDeff_median_all$tissues, levels = rev(as.character(tissues_NMDeff_median_all_medians$tissues)))
    
  # 3.3) Violin plot of differences between SD (left panel) and medians (right panel) for each tissue
  
  # tissues_NMDeff_median_all$type <- "NMDeff_median"
  # tissues_SD_all$type <- "SD"
  tissues_NMDeff_median_and_SD <- merge(tissues_NMDeff_median_all, tissues_SD_all, by = c("tissues"), all.x = TRUE)
  tissues_NMDeff_median_and_SD <- cbind(tissues_NMDeff_median_all, tissues_SD_all)
  colnames(tissues_NMDeff_median_and_SD)[1] <- "tissues"
  # Stack df
  tissues_NMDeff_median_and_SD_df <- tissues_NMDeff_median_and_SD[,c("diff_NMDeff_median","diff_SD")]
  tissues_NMDeff_median_and_SD_df_stack <- stack(tissues_NMDeff_median_and_SD_df)
  tissues_NMDeff_median_and_SD_df_stack[1:(nrow(tissues_NMDeff_median_and_SD_df_stack) / 2),"p_value_below_NMDeff_median"] <- tissues_NMDeff_median_and_SD$p_value_below_NMDeff_median
  tissues_NMDeff_median_and_SD_df_stack[1:(nrow(tissues_NMDeff_median_and_SD_df_stack) / 2),"p_value_above_NMDeff_median"] <- tissues_NMDeff_median_and_SD$p_value_above_NMDeff_median
  tissues_NMDeff_median_and_SD_df_stack[( ((nrow(tissues_NMDeff_median_and_SD_df_stack) / 2) + 1) : nrow(tissues_NMDeff_median_and_SD_df_stack) ),"p_value_below_SD"] <- tissues_NMDeff_median_and_SD$p_value_below_SD
  tissues_NMDeff_median_and_SD_df_stack[( ((nrow(tissues_NMDeff_median_and_SD_df_stack) / 2) + 1) : nrow(tissues_NMDeff_median_and_SD_df_stack)),"p_value_above_SD"] <- tissues_NMDeff_median_and_SD$p_value_above_SD
  tissues_NMDeff_median_and_SD_df_stack$tissues <- rep(tissues_NMDeff_median_and_SD$tissues,2)
  # Color
  tissues_NMDeff_median_and_SD_df_stack <- tissues_NMDeff_median_and_SD_df_stack %>%
        mutate(color_NMDeff_median = ifelse(p_value_below_NMDeff_median < 0.05 | p_value_above_NMDeff_median < 0.05,TRUE,FALSE)) %>%
        mutate(color_SD = ifelse(p_value_below_SD < 0.05 | p_value_above_SD < 0.05,TRUE,FALSE)) %>%
        mutate(ind = ifelse(ind == "diff_NMDeff_median", "NMDeff median","SD"))
        
  tissues_NMDeff_median_and_SD_df_stack$color <- c(tissues_NMDeff_median_and_SD_df_stack[1:(nrow(tissues_NMDeff_median_and_SD_df_stack) / 2),"color_NMDeff_median"],
                                                    tissues_NMDeff_median_and_SD_df_stack[( ((nrow(tissues_NMDeff_median_and_SD_df_stack) / 2) + 1) : nrow(tissues_NMDeff_median_and_SD_df_stack)),"color_SD"])

  # Plot
  plot_path <- paste0(paste0(gsub("NMD_genesets","tissues_SD",plot_path_tmp)) , NMD_geneset_name,"_NMDeff_median_and_SD_by_tissue.png")  
  png(plot_path, width = 4250, height = 4500, res = 300)
  p <- ggplot(data = tissues_NMDeff_median_and_SD_df_stack, aes(x = tissues, y = values, color = color)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = ylim_SD) +
    geom_point(alpha=0.25) + xlab("") +
    facet_wrap(. ~ ind) +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    ylab("Observed - Randomized Difference") +
    theme_bw(base_size = 15) + guides(color = guide_legend(title = "P-value < 0.05")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=26, face="bold"),
          axis.text.x = element_text(color="black", size=26),
          axis.text.y = element_text(color="black", size=22),
          strip.text = element_text(size = 28),
          panel.spacing = unit(2, "cm"),
          legend.position = "top", 
          legend.text = element_text(size = 24),
          legend.title=element_text(size=26)) +
    scale_color_brewer(palette = "Accent", direction = -1)
  print(p)
  dev.off()

  # 4) Plot Observed NMDeff + Randomized NMDeff 
  rand_NMDeff_iteration_list  <- lapply(rand_scores_list, function(X){X[["NMDeff_rand"]]})
  NMD_efficiencies_df_filt$NMD_eff_rand <- NA
  for (sample in 1:nrow(rand_NMDeff_iteration_list[[1]])) {
    #print(sample)
    sample_NMDeff_list <- lapply(rand_NMDeff_iteration_list, function(NMDeff_rand_it){
      NMDeff_rand_it[sample,"values"]
      })
    sample_NMDeff <- unlist(sample_NMDeff_list)
    # Take one NMD eff random
    NMD_efficiencies_df_filt[sample,"NMD_eff_rand"] <- sample(na.omit(sample_NMDeff))[1]
  }

  NMD_efficiencies_df_filt <- na.omit(NMD_efficiencies_df_filt)
  NMD_efficiencies_df_merged <- NMD_efficiencies_df_filt[,!colnames(NMD_efficiencies_df_filt)%in%NMD_geneset]
  colnames(NMD_efficiencies_df_merged) <- c("cancer","NMD_efficiency","id","sample_size","tissue_ss","NMD_efficiency_rand")
  # Sort Randomized by median
  NMD_efficiencies_df_merged_rand <- NMD_efficiencies_df_merged[,c("tissue_ss","NMD_efficiency_rand")]
  NMDeff_medians_df_rand <- aggregate(NMD_efficiency_rand ~ tissue_ss,NMD_efficiencies_df_merged_rand, median)
  NMDeff_medians_df_rand <- NMDeff_medians_df_rand[order(NMDeff_medians_df_rand[,2]),]
  NMD_efficiencies_df_merged_rand_sort <- data.frame()
  for (i in 1:nrow(NMDeff_medians_df_rand)) {
    cancer <- NMDeff_medians_df_rand$tissue_ss[i]
    cancer_obs <- levels(NMD_efficiencies_df_filt$tissue_ss)[i]
    NMD_efficiencies_df_merged_rand_filt <- NMD_efficiencies_df_merged_rand[NMD_efficiencies_df_merged_rand$tissue_ss %in% cancer,]
    NMD_efficiencies_df_merged_rand_filt$tissue_ss <- cancer_obs
    if (nrow(NMD_efficiencies_df_merged_rand_sort) == 0) {
      NMD_efficiencies_df_merged_rand_sort <- NMD_efficiencies_df_merged_rand_filt
    } else {
      NMD_efficiencies_df_merged_rand_sort <- rbind(NMD_efficiencies_df_merged_rand_sort,NMD_efficiencies_df_merged_rand_filt)
    }
  }
  # Merge with previous NMDeff
  NMD_efficiencies_df_final <- NMD_efficiencies_df_merged[,c("tissue_ss","NMD_efficiency")]
  colnames(NMD_efficiencies_df_final) <- c("tissue_ss","NMD_efficiency_rand")
  NMD_efficiencies_df_final_sort <- rbind(NMD_efficiencies_df_final,NMD_efficiencies_df_merged_rand_sort)
  colnames(NMD_efficiencies_df_final_sort) <- c("tissue_ss","NMD_efficiency")
  NMD_efficiencies_df_final_sort$type <- NA
  NMD_efficiencies_df_final_sort[1:nrow(NMD_efficiencies_df_merged),"type"] <- "NMD_efficiency"
  NMD_efficiencies_df_final_sort[(nrow(NMD_efficiencies_df_merged)+1):(nrow(NMD_efficiencies_df_final_sort)),"type"] <- "Randomization"
  #NMD_efficiencies_df_final_sort$tissue_ss <- factor(NMD_efficiencies_df_final_sort_tmp$tissue_ss, levels = c(levels(NMD_efficiencies_df_filt$tissue_ss)))
  
  plot_path1 <- gsub("\\[X\\]",NMD_geneset_name,paste0(plot_path_tmp,paths[paths$folder_or_object==paste0(dataset,"_",method,"_boxplots_NMD_genesets"),"path_or_filename"]))
  plot_path2 <- gsub("\\[X2\\]",plot_name1,plot_path1)
  png(plot_path2, width = 4500, height = 3500, res = 300)
  if (method == "consensus") {
    p <- ggplot(data = NMD_efficiencies_df_final_sort, aes(x=tissue_ss, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","NMD_efficiency")))) +
      geom_boxplot() + coord_flip(ylim = ylim) +
      geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
      ylab("NMD efficiency") + xlab("tissues") + ggtitle(paste0("NMD efficiency --> ", NMD_geneset)) +
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text.x = element_text(color="black", size=x.size),
            axis.text.y = element_text(color="black", size=18),
            legend.position='top')
  } else {
    title <- gsub("endogenous ","",paste0(gsub("_"," ",NMD_geneset)))
    p <- ggplot(data = NMD_efficiencies_df_final_sort, aes(x=tissue_ss, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","NMD_efficiency")))) +
      geom_boxplot() + coord_flip(ylim = c(-5,5)) + #coord_cartesian(ylim = c(-5,5)) +
      guides(fill = guide_legend(title = "")) +
      geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
      ylab("NMD efficiency") + ggtitle(title) + xlab("") +
      theme_bw() + scale_fill_brewer(palette = "Dark2") +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text.x = element_text(color="black", size=x.size),# angle = 90, hjust = 1),
            axis.text.y = element_text(color="black", size=18),
            legend.position='top', legend.text = element_text(size = 18))
  }
  print(p)
  dev.off()

  return(list(tissues_NMDeff_SD = tissues_NMDeff_median_and_SD_df_stack, tissues_NMDeff_with_rand = NMD_efficiencies_df_final_sort))
}

NMD_geneset_boxplots_rand_loop <- function(NMD_efficiencies_df, dataset, method, NMD_genesets, n_rand, rand_iteration) {

  if (dataset == "GTEx") {
    tissue <- "acronyms"
  } else if (dataset == "TCGA") {
    tissue <- "cancer_type_strat"
  }
  # 1)
  randomization_pvals <- matrix(nrow=0, ncol=length(NMD_genesets))
  colnames(randomization_pvals) <- NMD_genesets
  randomization_pvals <- as.data.frame(randomization_pvals)
  # 2) Inter-tissue randomization test initialization
  sd_of_medians_diff <- matrix(nrow=length(NMD_genesets), ncol=4)
  rownames(sd_of_medians_diff) <- NMD_genesets
  colnames(sd_of_medians_diff) <- c("median_SD_of_medians_diff","perc5_SD_of_medians_diff", "perc95_SD_of_medians_diff","p_value")
  sd_of_medians_diff <- as.data.frame(sd_of_medians_diff)
  randomization_sd_of_medians <- matrix(nrow = n_rand, ncol = length(NMD_genesets))
  colnames(randomization_sd_of_medians) <- NMD_genesets
  randomization_sd_of_medians <- as.data.frame(randomization_sd_of_medians)
  obs_sd_of_medians <- matrix(nrow = 1, ncol = length(NMD_genesets))
  colnames(obs_sd_of_medians) <- NMD_genesets
  obs_sd_of_medians <- as.data.frame(obs_sd_of_medians)
  # 3)
  randomization_tissue_sd <- matrix(nrow=length(unique(NMD_efficiencies_df[,tissue])), ncol=length(NMD_genesets)*4)
  rownames(randomization_tissue_sd) <- unique(NMD_efficiencies_df[,tissue])
  colnames(randomization_tissue_sd) <- sort(c(paste0(NMD_genesets,"_log2fold"),paste0(NMD_genesets,"_obs_sd"),paste0(NMD_genesets,"_rand_sd"),paste0(NMD_genesets,"_p_val")))
    
  randomization_list[[dataset]][[method]] <<- list(
        randomization_pvals = randomization_pvals, 
        randomization_sd_of_medians = list(sd_of_medians_diff = sd_of_medians_diff, randomization_sd_of_medians = randomization_sd_of_medians), 
        randomization_tissue_sd = randomization_tissue_sd, obs_sd_of_medians = obs_sd_of_medians
        )

  res_list_genesets <- list()
  for (NMD_geneset in NMD_genesets) {
    res_list <- boxplots_by_NMDgeneset(NMD_geneset = NMD_geneset, dataset = dataset, method = method, 
                            NMD_efficiencies_df = NMD_efficiencies_df, n_rand = n_rand)
    res_list_genesets[[NMD_geneset]] <- res_list
  }
  
  # Plot path
  plot_path <- eval(parse(text=paste0(dataset,"_boxplots_",method,"_NMD_genesets_path")))

  # Adjust by FDR
  randomization_list[[dataset]][[method]]$randomization_pvals <- apply(randomization_list[[dataset]][[method]]$randomization_pvals,2,function(x){p.adjust(x,method = "fdr")})
  
  write.table(randomization_list[[dataset]][[method]]$randomization_pvals, file = gsub(".txt",paste0("_random_iteration_",rand_iteration,".txt"),paste0(plot_path, paths[paths$folder_or_object==paste0(dataset,"_",method,"_randomization_pvals"),"path_or_filename"])), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  saveRDS(randomization_list[[dataset]][[method]]$randomization_sd_of_medians,gsub(".txt",paste0("_random_iteration_",rand_iteration,".RData"),paste0(plot_path, paths[paths$folder_or_object==paste0(dataset,"_",method,"_randomization_median_SD"),"path_or_filename"])))
  # write.table(randomization_list[[dataset]][[method]]$randomization_sd_of_medians, file = gsub(".txt",paste0("_random_iteration_",rand_iteration,".txt"),paste0(plot_path, paths[paths$folder_or_object==paste0(dataset,"_",method,"_randomization_median_SD"),"path_or_filename"])), 
  #             sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(randomization_list[[dataset]][[method]]$randomization_tissue_sd, file = gsub(".txt",paste0("_random_iteration_",rand_iteration,".txt"),paste0(plot_path, paths[paths$folder_or_object==paste0(dataset,"_",method,"_randomization_tissues_SD"),"path_or_filename"])), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  return(res_list_genesets)
}

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# Arguments and paths

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/analysis_results_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#args <- commandArgs(trailingOnly=TRUE)
#paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
TMB_TCGA_path <- paths[paths$folder_or_object=="TMB_TCGA_path","path_or_filename"]
TCGA_names_path <- paths[paths$folder_or_object=="TCGA_names_path","path_or_filename"]
TCGA_samples_metadata_path <- paths[paths$folder_or_object=="TCGA_metadata_path","path_or_filename"]
GTEx_names_path <- paths[paths$folder_or_object=="GTEx_names_path","path_or_filename"]
GTEx_NB_res_tissue_path <- paths[paths$folder_or_object=="GTEx_NB_res_tissue_path","path_or_filename"]

# Pancancer Boxplots
TCGA_boxplot_endogenous_pancancer_path <- paths[paths$folder_or_object=="TCGA_boxplot_endogenous_pancancer_path","path_or_filename"]
TCGA_boxplot_PTCs_pancancer_path <- paths[paths$folder_or_object=="TCGA_boxplot_PTCs_pancancer_path","path_or_filename"]
TCGA_boxplot_ASE_pancancer_path <- paths[paths$folder_or_object=="TCGA_boxplot_ASE_pancancer_path","path_or_filename"]
GTEx_boxplot_endogenous_pantissue_path <- paths[paths$folder_or_object=="GTEx_boxplot_endogenous_pantissue_path","path_or_filename"]
GTEx_boxplot_ASE_pantissue_path <- paths[paths$folder_or_object=="GTEx_boxplot_ASE_pantissue_path","path_or_filename"]
GTEx_boxplot_PTCs_pantissue_path <- paths[paths$folder_or_object=="GTEx_boxplot_PTCs_pantissue_path","path_or_filename"]

# NMD geneset boxplots
TCGA_boxplots_endogenous_NMD_genesets_path <- paths[paths$folder_or_object=="TCGA_boxplots_endogenous_NMD_genesets_path","path_or_filename"]
TCGA_boxplots_PTCs_NMD_genesets_path <- paths[paths$folder_or_object=="TCGA_boxplots_PTCs_NMD_genesets_path","path_or_filename"]
TCGA_boxplots_ASE_NMD_genesets_path <- paths[paths$folder_or_object=="TCGA_boxplots_ASE_NMD_genesets_path","path_or_filename"]
TCGA_boxplots_consensus_NMD_genesets_path <- paths[paths$folder_or_object=="TCGA_boxplots_consensus_NMD_genesets_path","path_or_filename"]
GTEx_boxplots_endogenous_NMD_genesets_path <- paths[paths$folder_or_object=="GTEx_endogenous_boxplots_NMD_genesets_path","path_or_filename"]
GTEx_boxplots_ASE_NMD_genesets_path <- paths[paths$folder_or_object=="GTEx_ASE_boxplots_NMD_genesets_path","path_or_filename"]
GTEx_boxplots_PTCs_NMD_genesets_path <- paths[paths$folder_or_object=="GTEx_PTCs_boxplots_NMD_genesets_path","path_or_filename"]
GTEx_boxplots_consensus_NMD_genesets_path <- paths[paths$folder_or_object=="GTEx_boxplots_consensus_NMD_genesets_path","path_or_filename"]

# Tissues boxplots all NMD genesets
# TCGA_endogenous_boxplots_tissues_path <- paths[paths$folder_or_object=="TCGA_endogenous_boxplots_tissues_path","path_or_filename"]
# TCGA_PTCs_boxplots_tissues_path <- paths[paths$folder_or_object=="TCGA_PTCs_boxplots_tissues_path","path_or_filename"]
# TCGA_ASE_boxplots_tissues_path <- paths[paths$folder_or_object=="TCGA_ASE_boxplots_tissues_path","path_or_filename"]
# GTEx_endogenous_boxplots_tissues_path <- paths[paths$folder_or_object=="GTEx_endogenous_boxplots_tissues_path","path_or_filename"]
# GTEx_ASE_boxplots_tissues_path <- paths[paths$folder_or_object=="GTEx_ASE_boxplots_tissues_path","path_or_filename"] 

# 1) Data

# 1.1) sample NMD efficiencies TCGA

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

endogenous_NMD_genesets <-  c("NMD Colombo","NMD Karousis","NMD Tani","NMD Courtney","NMD Ensembl",
                      "NMD All","NMD Consensus","NMD SMG6","NMD SMG7",
                      "RandomGenes without NMD features","RandomGenes with NMD features")
ASE_NMD_genesets <- c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01","Synonymous 0.01",
                      "PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2")

# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# # 1.2) PTC NMD efficiencies TCGA
# # PTC
# PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_TCGA_dataset/germline_PTCs_all_TCGA_confident.txt"
# PTCs_NMD_efficiencies_TCGA <- read.table(file = PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # ASE
# PTCs_ASE_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_TCGA_dataset/germline_PTCs_ASE_all_TCGA_confident.txt"
# PTCs_ASE_NMD_efficiencies_TCGA <- read.table(file = PTCs_ASE_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # Add sample metadata
# PTCs_ASE_NMD_efficiencies_TCGA <- merge(PTCs_ASE_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)
# PTCs_NMD_efficiencies_TCGA <- merge(PTCs_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)

# 1.3) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = FALSE)

# 1.4) Missing samples
# Sample size TCGA
# TCGA_RNAseq_number_samples <- read.table(file = paste0(TCGA_names_path,paths[paths$folder_or_object=="TCGA_number_samples","path_or_filename"]),
#                                 header = TRUE, stringsAsFactors = FALSE)
# rownames(TCGA_RNAseq_number_samples) <- paste0("TCGA-",TCGA_RNAseq_number_samples$cancer)
# cancer_types <- rownames(TCGA_RNAseq_number_samples)
# TCGA_number_samples <- data.frame(cancer = cancer_types, total_samples = NA, 
#                                   samples_endogenous = NA, error_samples_endogenous = NA, missing_samples_endogenous = NA,
#                                   samples_PTCs = NA, error_samples_PTCs = NA, missing_samples_PTCs = NA,
#                                   samples_ASE = NA, error_samples_ASE = NA, missing_samples_ASE = NA)
# for (i in 1:nrow(TCGA_number_samples)) {
#   cancer <- as.character(TCGA_number_samples[i,"cancer"])
#   cancer_original <- gsub("(TCGA-\\w{2,4})_.*","\\1",cancer)
#   RNAseq_sample_size <- TCGA_RNAseq_number_samples[rownames(TCGA_RNAseq_number_samples)%in%cancer_original,"RNAseq_samples"]
#   TCGA_number_samples[i,"total_samples"] <- RNAseq_sample_size
#   filter <- sample_NMD_efficiencies_TCGA$cancer_type %in% cancer_original
#   # Endogenous
#   NMD_eff <- sample_NMD_efficiencies_TCGA[filter,"endogenous_NMD_global_2_shared"]
#   TCGA_number_samples[i,"samples_endogenous"] <- length(na.omit(NMD_eff))
#   TCGA_number_samples[i,"error_samples_endogenous"] <- sum(is.na(NMD_eff))
#   TCGA_number_samples[i,"missing_samples_endogenous"] <- round(1-length(na.omit(NMD_eff))/RNAseq_sample_size,2)
#   # PTCs
#   NMD_eff <- sample_NMD_efficiencies_TCGA[filter,"PTCs_stopgain_NMD_triggering"]
#   TCGA_number_samples[i,"samples_PTCs"] <- length(na.omit(NMD_eff))
#   TCGA_number_samples[i,"error_samples_PTCs"] <- sum(is.na(NMD_eff))
#   TCGA_number_samples[i,"missing_samples_PTCs"] <- round(1-length(na.omit(NMD_eff))/RNAseq_sample_size,2)
#   # ASE
#   NMD_eff <- sample_NMD_efficiencies_TCGA[filter,"ASE_stopgain_0.2"]
#   TCGA_number_samples[i,"samples_ASE"] <- length(na.omit(NMD_eff))
#   TCGA_number_samples[i,"error_samples_ASE"] <- sum(is.na(NMD_eff))
#   TCGA_number_samples[i,"missing_samples_ASE"] <- round(1-length(na.omit(NMD_eff))/RNAseq_sample_size,2)
# } 

# original_cancers <- TCGA_number_samples$cancer%in%rownames(TCGA_RNAseq_number_samples)
# TCGA_number_samples_filt <- TCGA_number_samples[original_cancers,c("cancer","missing_samples_endogenous","missing_samples_PTCs","missing_samples_ASE")]
# TCGA_number_samples_stack <- stack(TCGA_number_samples_filt)
# TCGA_number_samples_stack$tissues <- rep(TCGA_number_samples_filt$cancer,4)
# TCGA_number_samples_stack <- TCGA_number_samples_stack[!TCGA_number_samples_stack$ind %in% "cancer",]

# png(paste0(TCGA_NB_res_strat_list_path,paths[paths$folder_or_object=="TCGA_missing_samples_boxplot","path_or_filename"]), width = 4000, height = 3000, res = 300)
# p <- ggplot(data = TCGA_number_samples_stack, aes(x = tissues, y = values, fill = ind)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   #geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
#   ylab("% missing samples") + xlab("tissues") + ggtitle(paste0("Barplot")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         axis.text.x = element_text(color="black", size=10, angle = 75),
#         axis.title.y = element_text(color="black", size=13, face="bold"),
#         axis.title.x = element_text(color="black", size=13, face="bold"))
# print(p)
# dev.off()

# 2.1) Pancancer 
NMD_genesets_consensus <- c("NMDeff_mean")
#NMD_genesets_consensus <- c("NMDeff_mean","NMDeff_PCA")
# NMD efficiency of all NMD genesets
### TCGA ###
# df_results <- pantissue_boxplot_NMDeff(dataset = "TCGA", method = "endogenous", NMD_genesets = endogenous_NMD_genesets)
# write.table(df_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1A.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_results, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1A.RData")
# write.table(df_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2A.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_results, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2A.RData")
# # df_results <- pantissue_boxplot_NMDeff(dataset = "TCGA", method = "ASE", NMD_genesets = ASE_NMD_genesets[grep("0.01",ASE_NMD_genesets)], VAF = "0.01")
# df_results <- pantissue_boxplot_NMDeff(dataset = "TCGA", method = "ASE", NMD_genesets = ASE_NMD_genesets[grep("0.2",ASE_NMD_genesets)], VAF = "0.2")
# write.table(df_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1B.txt", 
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# #pantissue_boxplot_NMDeff(dataset = "TCGA", method = "PTCs", NMD_genesets = PTC_sets)
# ### GTEx ###
# df_results <- pantissue_boxplot_NMDeff(dataset = "GTEx", method = "endogenous", NMD_genesets = endogenous_NMD_genesets)
# write.table(df_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2B.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_results, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2B.RData")
# # df_results <- pantissue_boxplot_NMDeff(dataset = "GTEx", method = "ASE", NMD_genesets = ASE_NMD_genesets[grep("0.01",ASE_NMD_genesets)], VAF = "0.01")
# df_results <- pantissue_boxplot_NMDeff(dataset = "GTEx", method = "ASE", NMD_genesets = ASE_NMD_genesets[grep("0.2",ASE_NMD_genesets)], VAF = "0.2")
# write.table(df_results, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2C.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_results, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2C.RData")
# #pantissue_boxplot_NMDeff(dataset = "GTEx", method = "PTCs", NMD_genesets = PTC_sets)

# 2.2) By tissues

# 2.2.1) NMD efficiency separately for each NMD geneset and stratified by tissues + its distribution + randomization test

randomization_list <- list()
# Randomize 3 times
n <- 2000
for (random_iteration in 1:3) {
  print(paste0("-----------------RANDOM ITERATION #",random_iteration,"-----------------"))
  endogenous_NMD_genesets <- c("NMD All","NMD Consensus", "RandomGenes with NMD features","RandomGenes without NMD features")
  TCGA_end_variability <- NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_TCGA, dataset = "TCGA", method = "endogenous", 
                            NMD_genesets = endogenous_NMD_genesets, n_rand = n, rand_iteration = random_iteration)
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_end_variability_random_iteration_",random_iteration,".rds")
  saveRDS(TCGA_end_variability, output_path)
  TCGA_ASE_variability <- NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_TCGA, dataset = "TCGA", method = "ASE", 
                            NMD_genesets = ASE_NMD_genesets[grep("0.2",ASE_NMD_genesets)], n_rand = n, rand_iteration = random_iteration)
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_ASE_variability_random_iteration_",random_iteration,".rds")
  saveRDS(TCGA_ASE_variability, output_path)
  # Subsample 139 rows per tissue
  # set.seed(123)
  # subsampled_data <- sample_NMD_efficiencies_GTEx %>%
  #   group_by(acronyms) %>%
  #   slice_sample(n = 139) %>%
  #   ungroup()
  # summary(subsampled_data[,"NMD Consensus"])
  GTEx_end_variability <- NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = subsampled_data, dataset = "GTEx", method = "endogenous", 
                            NMD_genesets = endogenous_NMD_genesets, n_rand = n, rand_iteration = random_iteration)
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_end_variability_random_iteration_",random_iteration,".rds")
  # output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_end_variability_random_iteration_",random_iteration,"_brain_subsampling.rds")
  saveRDS(GTEx_end_variability, output_path)
  GTEx_ASE_variability <- NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = subsampled_data, dataset = "GTEx", method = "ASE", 
                            NMD_genesets = ASE_NMD_genesets[grep("0.2",ASE_NMD_genesets)], n_rand = n, rand_iteration = random_iteration)
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_ASE_variability_random_iteration_",random_iteration,".rds")
  # output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_ASE_variability_random_iteration_",random_iteration,"_brain_subsampling.rds")
  saveRDS(GTEx_ASE_variability, output_path)
}

stop("ADIOS")

#NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_TCGA, dataset = "TCGA", method = "PTCs", NMD_genesets = PTC_sets, n_rand = n)
#NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_TCGA, dataset = "TCGA", method = "consensus", NMD_genesets = NMD_genesets_consensus, n_rand = n)
#NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_GTEx, dataset = "GTEx", method = "PTCs", NMD_genesets = PTC_sets, n_rand = n)
#NMD_geneset_boxplots_rand_loop(NMD_efficiencies_df = sample_NMD_efficiencies_GTEx, dataset = "GTEx", method = "consensus", NMD_genesets = NMD_genesets_consensus, n_rand = n)

# 3) Plots

# 3.1) Inter-tissue variability NMDeff Obs-Rand difference
NMD_method_End <- "NMD Consensus"
NMD_method_ASE <- "PTC NMD-triggering 0.2"

TCGA_pan_organ_system <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_pan_organ_system.txt", 
          header = TRUE, stringsAsFactors = FALSE)
TCGA_pan_organ_system$cancer_type <- gsub("TCGA-","",TCGA_pan_organ_system$cancer_type)

## TMP ##
TCGA_end_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_end_variability_random_iteration_1.rds")
TCGA_ASE_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_ASE_variability_random_iteration_1.rds")
#########

# 3.1.1) TCGA
TCGA_end_NMDeff_diff <- TCGA_end_variability[[NMD_method_End]][["tissues_NMDeff_SD"]]
TCGA_end_NMDeff_diff$NMD_method <- "ETG"
TCGA_ASE_NMDeff_diff <- TCGA_ASE_variability[[NMD_method_ASE]][["tissues_NMDeff_SD"]]
TCGA_ASE_NMDeff_diff$NMD_method <- "ASE"
TCGA_NMDeff_diff <- rbind(TCGA_end_NMDeff_diff,TCGA_ASE_NMDeff_diff)
TCGA_NMDeff_diff$NMD_method <- factor(TCGA_NMDeff_diff$NMD_method, levels = c("ETG","ASE"))
TCGA_NMDeff_diff <- TCGA_NMDeff_diff[TCGA_NMDeff_diff$ind == "NMDeff median",]
# Organ system
TCGA_NMDeff_diff$pan_organ_system <- "Other"
TCGA_NMDeff_diff[TCGA_NMDeff_diff$tissues %in% c("KIRC","KIRP"),"pan_organ_system"] <- "Pan-Kidney"
TCGA_NMDeff_diff[grep("COAD|STAD|READ|ESCA",TCGA_NMDeff_diff$tissue),"pan_organ_system"] <- "Pan-GI"
TCGA_NMDeff_diff[grep("BRCA|OV|UCS|UCEC|CESC|TGCT|PRAD",TCGA_NMDeff_diff$tissue),"pan_organ_system"] <- "Pan-gyn"
TCGA_NMDeff_diff[TCGA_NMDeff_diff$tissues %in% c("LGG","PCPG","GBM"),"pan_organ_system"] <- "Brain"
TCGA_NMDeff_diff[grep("HNSC|LUSC|BLCA",TCGA_NMDeff_diff$tissue),"pan_organ_system"] <- "Pan-Squamous"
TCGA_NMDeff_diff_filt <- TCGA_NMDeff_diff[TCGA_NMDeff_diff$pan_organ_system != "Other",]
# Order
df_organs <- aggregate(values ~ pan_organ_system + NMD_method, data = TCGA_NMDeff_diff_filt , median)
df_organs <- df_organs[df_organs$NMD_method == "ETG",]
df_organs <- df_organs %>%
      arrange(values)
df_tissues <- aggregate(values ~ pan_organ_system + tissues + NMD_method, data = TCGA_NMDeff_diff_filt , median)
df_tissues <- df_tissues[df_tissues$NMD_method == "ETG",]
df_tissues$pan_organ_system <- factor(df_tissues$pan_organ_system, levels = df_organs$pan_organ_system)
df_tissues <- df_tissues[order(df_tissues$pan_organ_system),]
df_tissues <- df_tissues %>%
      arrange(pan_organ_system,values)
TCGA_NMDeff_diff_filt$tissues <- factor(TCGA_NMDeff_diff_filt$tissues, levels = as.character(df_tissues$tissues))
TCGA_NMDeff_diff_filt$pan_organ_system <- factor(TCGA_NMDeff_diff_filt$pan_organ_system, levels = unique(df_tissues$pan_organ_system))

# p <- ggplot(data = TCGA_NMDeff_diff, aes(x = tissues, y = values, color = color)) +
p <- TCGA_NMDeff_diff_filt %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.4,0.4)) +
    #geom_point(alpha=0.25) + xlab("") +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    ylab("Observed - Randomized Difference") + xlab("") +
    theme_bw(base_size = 15) + #guides(color = guide_legend(title = "P-value < 0.05")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=30),
          axis.text.x = element_text(color="black", size=30),
          axis.text.y = element_text(color="black", size=28),
          strip.text = element_text(size = 30),
          panel.spacing = unit(2, "cm"),
          legend.position = "top", 
          legend.text = element_text(size = 26),
          legend.title=element_text(size=28)) +
    scale_color_brewer(palette = "Accent", direction = -1) +
    guides(fill = guide_legend(title = "Pan-Organ System",nrow = 3,override.aes = list(size = 12)))
plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_NMDeff_median_diff_by_tissue.png")  
png(plot_path, width = 4250, height = 4500, res = 300)
print(p)
dev.off()

# Save
write.table(TCGA_NMDeff_diff_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_NMDeff_diff_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2C.RData")

# Order
df_organs <- aggregate(values ~ pan_organ_system + NMD_method, data = TCGA_NMDeff_diff , median)
df_organs <- df_organs[df_organs$NMD_method == "ETG",]
df_organs <- df_organs %>%
      arrange(values)
df_tissues <- aggregate(values ~ pan_organ_system + tissues + NMD_method, data = TCGA_NMDeff_diff , median)
df_tissues <- df_tissues[df_tissues$NMD_method == "ETG",]
df_tissues$pan_organ_system <- factor(df_tissues$pan_organ_system, levels = df_organs$pan_organ_system)
df_tissues <- df_tissues[order(df_tissues$pan_organ_system),]
df_tissues <- df_tissues %>%
      arrange(pan_organ_system,values)
TCGA_NMDeff_diff$tissues <- factor(TCGA_NMDeff_diff$tissues, levels = as.character(df_tissues$tissues))
TCGA_NMDeff_diff$pan_organ_system <- factor(TCGA_NMDeff_diff$pan_organ_system, levels = unique(df_tissues$pan_organ_system))

# Save
write.table(TCGA_NMDeff_diff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_NMDeff_diff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5A.RData")

# 3.2.2) GTEx

## TMP ##
GTEx_end_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_end_variability_random_iteration_1.rds")
GTEx_ASE_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_ASE_variability_random_iteration_1.rds")
#Subset
GTEx_end_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_end_variability_random_iteration_1_brain_subsampling.rds")
GTEx_ASE_variability <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_ASE_variability_random_iteration_1_brain_subsampling.rds")
#########

GTEx_end_NMDeff_diff <- GTEx_end_variability[[NMD_method_End]][["tissues_NMDeff_SD"]]
GTEx_end_NMDeff_diff$NMD_method <- "ETG"
GTEx_ASE_NMDeff_diff <- GTEx_ASE_variability[[NMD_method_ASE]][["tissues_NMDeff_SD"]]
GTEx_ASE_NMDeff_diff$NMD_method <- "ASE"
GTEx_NMDeff_diff <- rbind(GTEx_end_NMDeff_diff,GTEx_ASE_NMDeff_diff)
GTEx_NMDeff_diff$NMD_method <- factor(GTEx_NMDeff_diff$NMD_method, levels = c("ETG","ASE"))
GTEx_NMDeff_diff <- GTEx_NMDeff_diff[GTEx_NMDeff_diff$ind == "NMDeff median",]

# Organ system
GTEx_NMDeff_diff$pan_organ_system <- "Other"
GTEx_NMDeff_diff[GTEx_NMDeff_diff$tissues %in% c("KDNCTX","KDNMDL"),"pan_organ_system"] <- "Pan-Kidney"
GTEx_NMDeff_diff[GTEx_NMDeff_diff$tissues %in% c("CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","STMACH","SNTTRM"),"pan_organ_system"] <- "Pan-GI"
GTEx_NMDeff_diff[GTEx_NMDeff_diff$tissues %in% c("BREAST","UTERUS","TESTIS","PRSTTE","VAGINA","OVARY","UCS","CVXECT","CVSEND","FLLPNT"),"pan_organ_system"] <- "Pan-gyn"
GTEx_NMDeff_diff[GTEx_NMDeff_diff$tissues %in% c("NERVET","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG"),"pan_organ_system"] <- "Brain"
GTEx_NMDeff_diff[GTEx_NMDeff_diff$tissues %in% c("LUNG","BLDDER"),"pan_organ_system"] <- "Pan-Squamous"
GTEx_NMDeff_diff_filt <- GTEx_NMDeff_diff[GTEx_NMDeff_diff$pan_organ_system != "Other",]

# Order
df_organs <- aggregate(values ~ pan_organ_system + NMD_method, data = GTEx_NMDeff_diff_filt , median)
df_organs <- df_organs[df_organs$NMD_method == "ETG",]
df_organs <- df_organs %>%
      arrange(values)
df_tissues <- aggregate(values ~ pan_organ_system + tissues + NMD_method, data = GTEx_NMDeff_diff_filt , median)
df_tissues <- df_tissues[df_tissues$NMD_method == "ETG",]
df_tissues$pan_organ_system <- factor(df_tissues$pan_organ_system, levels = df_organs$pan_organ_system)
df_tissues <- df_tissues[order(df_tissues$pan_organ_system),]
df_tissues <- df_tissues %>%
      arrange(pan_organ_system,values)
GTEx_NMDeff_diff_filt$tissues <- factor(GTEx_NMDeff_diff_filt$tissues, levels = as.character(df_tissues$tissues))
GTEx_NMDeff_diff_filt$pan_organ_system <- factor(GTEx_NMDeff_diff_filt$pan_organ_system, levels = unique(df_tissues$pan_organ_system))

p <- GTEx_NMDeff_diff_filt %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.6,0.6)) +
    #geom_point(alpha=0.25) + xlab("") +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    ylab("Observed - Randomized Difference") + xlab("") +
    theme_bw(base_size = 15) + #guides(color = guide_legend(title = "P-value < 0.05")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=30),
          axis.text.x = element_text(color="black", size=30),
          axis.text.y = element_text(color="black", size=28),
          strip.text = element_text(size = 30),
          panel.spacing = unit(2, "cm"),
          legend.position = "top", 
          legend.text = element_text(size = 26),
          legend.title=element_text(size=28)) +
    scale_color_brewer(palette = "Accent", direction = -1) +
    guides(fill = guide_legend(title = "Pan-Organ System",nrow = 3,override.aes = list(size = 12)))
plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_NMDeff_median_diff_by_tissue.png")  
png(plot_path, width = 4250, height = 4500, res = 300)
print(p)
dev.off()

# Save
# write.table(GTEx_NMDeff_diff_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2D.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# saveRDS(GTEx_NMDeff_diff_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2D.RData")
#Subset
write.table(GTEx_NMDeff_diff_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig10/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_NMDeff_diff_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig10/panel_A.RData")

# Order
df_organs <- aggregate(values ~ pan_organ_system + NMD_method, data = GTEx_NMDeff_diff, median)
df_organs <- df_organs[df_organs$NMD_method == "ETG",]
df_organs <- df_organs %>%
      arrange(values)
df_tissues <- aggregate(values ~ pan_organ_system + tissues + NMD_method, data = GTEx_NMDeff_diff , median)
df_tissues <- df_tissues[df_tissues$NMD_method == "ETG",]
df_tissues$pan_organ_system <- factor(df_tissues$pan_organ_system, levels = df_organs$pan_organ_system)
df_tissues <- df_tissues[order(df_tissues$pan_organ_system),]
df_tissues <- df_tissues %>%
      arrange(pan_organ_system,values)
GTEx_NMDeff_diff$tissues <- factor(GTEx_NMDeff_diff$tissues, levels = as.character(df_tissues$tissues))
GTEx_NMDeff_diff$pan_organ_system <- factor(GTEx_NMDeff_diff$pan_organ_system, levels = unique(df_tissues$pan_organ_system))

# Save
write.table(GTEx_NMDeff_diff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_NMDeff_diff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5B.RData")


# 3.2) Inter-tissue variability NMDeff Obs and Rand in parallel
# 3.2.1) TCGA
TCGA_end_NMDeff_rand <- TCGA_end_variability[[NMD_method_End]][["tissues_NMDeff_with_rand"]]
TCGA_end_NMDeff_rand$NMD_method <- "ETG"
TCGA_ASE_NMDeff_rand <- TCGA_ASE_variability[[NMD_method_ASE]][["tissues_NMDeff_with_rand"]]
TCGA_ASE_NMDeff_rand$NMD_method <- "ASE"
TCGA_NMDeff_rand <- rbind(TCGA_end_NMDeff_rand,TCGA_ASE_NMDeff_rand)
TCGA_NMDeff_rand$NMD_method <- factor(TCGA_NMDeff_rand$NMD_method, levels = c("ETG","ASE"))
# Remove sample size
TCGA_NMDeff_rand$tissue <- gsub(" \\([0-9].*\\)","",TCGA_NMDeff_rand$tissue_ss)
# Order
NMDeff_medians <- aggregate(NMD_efficiency ~ tissue + NMD_method + type, data = TCGA_NMDeff_rand, median)
NMDeff_medians <- NMDeff_medians[grep("ETG",NMDeff_medians$NMD_method),]
NMDeff_medians <- NMDeff_medians[grep("NMD_efficiency",NMDeff_medians$type),]
NMDeff_medians <- NMDeff_medians[order(NMDeff_medians$V1),]
TCGA_NMDeff_rand$tissue <- factor(TCGA_NMDeff_rand$tissue, levels = unique(NMDeff_medians$tissue))
# Subset for Brain tissues
# TCGA_NMDeff_rand <- TCGA_NMDeff_rand[grep("LGG|GBM",TCGA_NMDeff_rand$tissue),]
#TCGA_NMDeff_rand <- TCGA_NMDeff_rand[grep("LGG|GBM|PCPG|TGCT|UCEC|OV|CESC",TCGA_NMDeff_rand$tissue),]

p <- ggplot(data = TCGA_NMDeff_rand, aes(x=tissue, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","NMD_efficiency")))) +
  geom_boxplot() + coord_flip(ylim = c(-2.5,2.5)) + #coord_cartesian(ylim = c(-5,5)) +
  facet_wrap(. ~ NMD_method) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("NMD efficiency") + ggtitle("") + xlab("") +
  theme_bw() + scale_fill_brewer(palette = "Dark2") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(color="black", size=25),
        strip.text = element_text(size = 30),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(color="black", size=25),# angle = 90, hjust = 1),
        axis.text.y = element_text(color="black", size=22),
        legend.position='top', legend.text = element_text(size = 25)) +
        guides(fill = guide_legend(title = "", override.aes = list(size = 12)))
# plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_NMDeff_and_rand_by_tissue_brain.png")  
plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/TCGA_NMDeff_and_rand_by_tissue.png")  
png(plot_path, width = 3500, height = 4500, res = 300)
print(p)
dev.off()

# Save
write.table(TCGA_NMDeff_rand, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4A_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_NMDeff_rand, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4A_B.RData")

# 3.2.2) GTEx
GTEx_end_NMDeff_rand <- GTEx_end_variability[[NMD_method_End]][["tissues_NMDeff_with_rand"]]
GTEx_end_NMDeff_rand$NMD_method <- "ETG"
GTEx_ASE_NMDeff_rand <- GTEx_ASE_variability[[NMD_method_ASE]][["tissues_NMDeff_with_rand"]]
GTEx_ASE_NMDeff_rand$NMD_method <- "ASE"
GTEx_NMDeff_rand <- rbind(GTEx_end_NMDeff_rand,GTEx_ASE_NMDeff_rand)
GTEx_NMDeff_rand$NMD_method <- factor(GTEx_NMDeff_rand$NMD_method, levels = c("ETG","ASE"))
# Remove sample size
GTEx_NMDeff_rand$tissue <- gsub(" \\([0-9].*\\)","",GTEx_NMDeff_rand$tissue_ss)
# Order
NMDeff_medians <- aggregate(NMD_efficiency ~ tissue + NMD_method + type, data = GTEx_NMDeff_rand, median)
NMDeff_medians <- NMDeff_medians[grep("ETG",NMDeff_medians$NMD_method),]
NMDeff_medians <- NMDeff_medians[grep("NMD_efficiency",NMDeff_medians$type),]
NMDeff_medians <- NMDeff_medians[order(NMDeff_medians$V1),]
GTEx_NMDeff_rand$tissue <- factor(GTEx_NMDeff_rand$tissue, levels = unique(NMDeff_medians$tissue))
# Subset for Brain tissues
#GTEx_NMDeff_rand <- GTEx_NMDeff_rand[grep("CV|BRN|NERVET|OVARY|UTERUS|PRSTTE|TESTIS",GTEx_NMDeff_rand$tissue),]
# GTEx_NMDeff_rand <- GTEx_NMDeff_rand[grep("BRN|NERVET",GTEx_NMDeff_rand$tissue),]

p <- ggplot(data = GTEx_NMDeff_rand, aes(x=tissue, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","NMD_efficiency")))) +
  geom_boxplot() + coord_flip(ylim = c(-2.5,2.5)) + #coord_cartesian(ylim = c(-5,5)) +
  facet_wrap(. ~ NMD_method) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  ylab("NMD efficiency") + ggtitle("") + xlab("") +
  theme_bw() + scale_fill_brewer(palette = "Dark2") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(color="black", size=35),
        strip.text = element_text(size = 30),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(color="black", size=30),# angle = 90, hjust = 1),
        axis.text.y = element_text(color="black", size=22),
        legend.position='top', legend.text = element_text(size = 35)) +
        guides(fill = guide_legend(title = "", override.aes = list(size = 14)))
plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/GTEx_NMDeff_and_rand_by_tissue_brain.png")  
png(plot_path, width = 4250, height = 4500, res = 300)
print(p)
dev.off()

# Save
write.table(GTEx_NMDeff_rand, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4C_D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_NMDeff_rand, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4C_D.RData")

# 3.2.3)  GTEx + TCGA
# TCGA_NMDeff_rand$dataset <- "TCGA"
# GTEx_NMDeff_rand$dataset <- "GTEx"
# All_NMDeff_rand <- rbind(TCGA_NMDeff_rand,GTEx_NMDeff_rand)
# # Order
# NMDeff_medians <- aggregate(NMD_efficiency ~ tissue + NMD_method + type, data = All_NMDeff_rand, median)
# NMDeff_medians <- NMDeff_medians[grep("ETG",NMDeff_medians$NMD_method),]
# NMDeff_medians <- NMDeff_medians[grep("NMD_efficiency",NMDeff_medians$type),]
# NMDeff_medians <- NMDeff_medians[order(NMDeff_medians$V1),]
# All_NMDeff_rand$tissue <- factor(All_NMDeff_rand$tissue, levels = unique(NMDeff_medians$tissue))

# p <- ggplot(data = All_NMDeff_rand, aes(x=tissue, y = NMD_efficiency, 
#       fill = factor(type, levels = c("Randomization","NMD_efficiency")), 
#       color = factor(dataset))
#       ) +
#   geom_boxplot(linewidth = 0.7) + coord_flip(ylim = c(-2.5,2.5)) + #coord_cartesian(ylim = c(-5,5)) +
#   facet_wrap(. ~ NMD_method) +
#   geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
#   ylab("NMD efficiency") + ggtitle("") + xlab("") +
#   theme_bw() + scale_fill_brewer(palette = "Dark2") +
#   scale_color_brewer(palette = "Accent") +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         axis.title.x = element_text(color="black", size=35),
#         strip.text = element_text(size = 30),
#         axis.title.y = element_text(color="black", size=20),
#         axis.text.x = element_text(color="black", size=30),# angle = 90, hjust = 1),
#         axis.text.y = element_text(color="black", size=22),
#         legend.position='top', legend.text = element_text(size = 35)) +
#         guides(fill = guide_legend(title = "", override.aes = list(size = 14)),
#               color = guide_legend(title = "", override.aes = list(size = 14)))
# plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/inter_tissue/All_NMDeff_and_rand_by_tissue_brain.png")  
# png(plot_path, width = 4250, height = 4500, res = 300)
# print(p)
# dev.off()











