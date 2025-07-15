#################################### FUNCTIONS ######################################
#################################### FUNCTIONS ######################################
#################################### FUNCTIONS ######################################

somatic_mut_NMDeff_associations <- function(NMD_method_VAF, ensembl_v88_gtf, PCs) {
  CGC_results_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/pancancer_",NMD_method_VAF,"_CGC_somatic_mut_CNV_PCs_",PCs,".txt")
  CGC_results <- read.csv(file = CGC_results_path, header = TRUE, stringsAsFactors = FALSE, sep ="\t")
  CGC_results <- merge(CGC_results, ensembl_v88_gtf, by.x = "Gene_symbol", by.y ="gene_name", all.x = TRUE)
  CGC_results$chr <- gsub("(.*)\\:.*","\\1",CGC_results$genome_location)
  return(CGC_results)
}

manhattan_fine_mapping <- function(somatic_mut_NMDeff_associations, chr, NMD_method_VAF, adjust = "no", type = "no") {

  # Filter by Chr X
  somatic_mut_NMDeff_associations <- somatic_mut_NMDeff_associations[which(somatic_mut_NMDeff_associations$chr %in% as.character(chr)),]
  # Cytogenetic locations of genes
  somatic_mut_NMDeff_associations <- merge(somatic_mut_NMDeff_associations, anno_gene[,c("ensembl_gene_id","band")], by.x = "ENSEMBL_gene", by.y = "ensembl_gene_id", all.x = TRUE)
  somatic_mut_NMDeff_associations$chr_arm <- factor(substr(somatic_mut_NMDeff_associations$band,1,1))
  # Order by Genome Location
  somatic_mut_NMDeff_associations$genome_location <- as.numeric(gsub(".*:(.*)-.*","\\1",somatic_mut_NMDeff_associations$genome_location))
  somatic_mut_NMDeff_associations <- somatic_mut_NMDeff_associations[order(somatic_mut_NMDeff_associations$genome_location),]
  # Remove NA
  somatic_mut_NMDeff_associations <- somatic_mut_NMDeff_associations[!is.na(somatic_mut_NMDeff_associations$som_CNV_amp_pval),]
  # Adjust pvalues
  error <- FALSE
  tryCatch({
    somatic_mut_NMDeff_associations[,"som_CNV_amp_qval"]  <- qvalue(somatic_mut_NMDeff_associations$som_CNV_amp_pval)$qvalues
  }, error = function(e){
      print(e)
      error <<- TRUE
    }
  )
  somatic_mut_NMDeff_associations[,"som_CNV_amp_pval_adjusted"] <- p.adjust(somatic_mut_NMDeff_associations$som_CNV_amp_pval,method="fdr")
  if (adjust == "yes") {
    if (type == "qval") {
      char <- "som_CNV_amp_pval_adjusted"
    } else if ( type == "FDR") {
      char <- "som_CNV_amp_pval_adjusted"
    } 
    if (isTRUE(error)) {
      return(NA)
    }
  } else {
    char <- "som_CNV_amp_pval"
    type <- "no"
  }

  # Plot 
  png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/manhattan_by_chr/",NMD_method_VAF,"_chr",chr,"_genes_pvals_adjust_",type,".png"), width = 4500, height = 3500, res = 300)
  p <- ggplot(data = somatic_mut_NMDeff_associations, aes(x = genome_location, y = -log10(eval(parse(text=char))), color = chr_arm)) +
    geom_point(size = 5) + # ggtitle(paste0("Spearman rank-coefficient --> ", paste0(round(cor_res$estimate,3))," / p-value = ",round(cor_res$p_value,3))) +
    xlab(paste0("Chr ",chr," Genome Location")) +
    ylab(paste0("-log10(P-value Adjusted ",type)) +
    #theme_classic() + 
    #geom_text_repel(aes(label=Gene_symbol),hjust=0.5, vjust=0.5, size = 6, max.overlaps = nrow(somatic_mut_NMDeff_associations), alpha = 1) +
    geom_text_repel(data = subset(somatic_mut_NMDeff_associations, eval(parse(text=char)) < 0.25 & !Gene_symbol %in% c("SMG5","SMG7")), 
                    aes(label = Gene_symbol), color = "black", size = 6) +
    geom_label_repel(data = subset(somatic_mut_NMDeff_associations, Gene_symbol %in% c("SMG5","SMG7")), 
                    aes(label = Gene_symbol), color = "black", size = 6) +
    geom_hline(yintercept=-log10(0.25), size = 2, color = "black", linetype="dashed") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          #axis.text.x = element_text(color="black", size=7, angle = 45),
          axis.text.x = element_blank(),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.y = element_text(color="black", size=20),
          legend.title=element_text(size=30),
          legend.position='top') #+ scale_x_discrete(guide=guide_axis(n.dodge=2))
  #geom_abline(intercept = 0 , slope = 1)
  print(p)
  dev.off()

}

genes_update <- function(TCGA_CNV_genes, ensembl_v88_gene_transcript_genesymbol) {

  # Add gene symbols / ensembl gene id to TCGA_CNV_genes
  genes_symbols <- unlist(lapply(strsplit(rownames(TCGA_CNV_genes),"\\|"), function(x){
    x[1]
  }))
  ensembl_gene_id <- unlist(lapply(strsplit(rownames(TCGA_CNV_genes),"\\|"), function(x){
    x[2]
  }))
  ensembl_gene_id <- gsub("(.*)\\..*","\\1",ensembl_gene_id)
  TCGA_CNV_genes$gene_symbols <- genes_symbols
  TCGA_CNV_genes$ensembl_gene_id <- ensembl_gene_id

  # 1.8.1) Update TCGA_CNV_genes ENSEMBLE GENE IDs (90%)
  #### Genes with no ENSEMBL gene ID are UNIQUE, so let's add it from our table
  missing_rows <- which(is.na(TCGA_CNV_genes$ensembl_gene_id))
  matching_genes <- ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name %in% TCGA_CNV_genes[,"gene_symbols"],c("gene_id","gene_name")]
  matching_genes <- matching_genes[!duplicated(matching_genes),]
  TCGA_CNV_genes_1 <- TCGA_CNV_genes[missing_rows,]
  TCGA_CNV_genes_1 <- merge(TCGA_CNV_genes_1,matching_genes, by.x = "gene_symbols", by.y = "gene_name", all.x = TRUE)
  TCGA_CNV_genes_1$ensembl_gene_id <- TCGA_CNV_genes_1$gene_id
  TCGA_CNV_genes_1$gene_id <- NULL
  TCGA_CNV_genes_1 <- TCGA_CNV_genes_1[-which(duplicated(TCGA_CNV_genes_1$gene_symbols)),]
  TCGA_CNV_genes_2 <- TCGA_CNV_genes[-missing_rows,]
  TCGA_CNV_genes_2 <- merge(TCGA_CNV_genes_2,matching_genes, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_2$gene_name <- NULL
  TCGA_CNV_genes_updated <- rbind(TCGA_CNV_genes_1,TCGA_CNV_genes_2)
  # Add Chromosome and genome location
  TCGA_CNV_genes_updated <- merge(TCGA_CNV_genes_updated,ensembl_v88_gtf_filt, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_updated <- TCGA_CNV_genes_updated[!is.na(TCGA_CNV_genes_updated$genome_location),]
  TCGA_CNV_genes_updated$genome_location_str <- substr(TCGA_CNV_genes_updated$genome_location,1,8)
  TCGA_CNV_genes_updated$chr <- gsub("(.*):.*","\\1",TCGA_CNV_genes_updated$genome_location)
  # Cytogenetic locations of genes
  TCGA_CNV_genes_updated <- merge(TCGA_CNV_genes_updated, anno_gene[,c("ensembl_gene_id","band")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  TCGA_CNV_genes_updated$chr_arm <- factor(substr(TCGA_CNV_genes_updated$band,1,1))
  # Colnames
  colnames(TCGA_CNV_genes_updated) <- gsub("\\.","-",colnames(TCGA_CNV_genes_updated))
  return(TCGA_CNV_genes_updated)

}

CNV_state_in_chr_for_PC <- function(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC, chr, sample_NMD_efficiencies_TCGA, quantiles = NULL) {

  list_results <- list()

  # 1) % of samples in each CNA-PC bin by cancer type
  # We want to know which cancer types have the most weighted samples
  TCGA_CNV_PC <- TCGA_CNV_PCA_ind[,paste0("Dim.",PC), drop = FALSE]

  # # Remove samples with 0 weights in the PC
  # TCGA_CNV_PC <- TCGA_CNV_PC[-which(TCGA_CNV_PC[,paste0("Dim.",PC)]  == 0),]
  # Obtain average CNV state per position for each bin of PC[X] samples
  TCGA_CNV_PC <- TCGA_CNV_PC[order(TCGA_CNV_PC[,1]),,drop=FALSE]
  thres <- quantile(TCGA_CNV_PC[,paste0("Dim.",PC)],seq(0,1,0.05))
  print(thres)
  thres["0%"] <- thres["0%"]-0.01
  thres["100%"] <- thres["100%"]+0.01
  if (is.null(quantiles)) {
    quantiles <- c("0%","33%","66%","100%")
  }
  bins <- cut(TCGA_CNV_PC[,paste0("Dim.",PC)], breaks = thres[quantiles])
  TCGA_CNV_PC$bins <- bins
  print(table(bins))

  if (PC %in% c(3,52)) {
    bin_order <- 1
    # Reverse and negate
    TCGA_CNV_PC <- TCGA_CNV_PC %>%
      mutate(bins = bins %>% 
              str_extract_all("[-]?[0-9.]+") %>%  # Extract numeric values
              map(~as.numeric(.x)) %>%            # Convert to numeric
              map(~-rev(.x)) %>%                  # Negate and reverse
              map_chr(~sprintf("[%.4f,%.4f)", .x[1], .x[2]))) # Construct new range
    bins <- factor(TCGA_CNV_PC$bins, levels = rev(names(table(TCGA_CNV_PC$bins))))
    TCGA_CNV_PC$bins <- bins 
    print(table(bins))
  } else if (PC == 86) {
    bin_order <- 3
  } else {
    bin_order <- 1
  }
  
  # Barplot of % of samples in each Bin for each cancer type
  NMDeff_PC <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_PC[,-1, drop = FALSE], by.x ="sample", by.y = "row.names", all.x = TRUE)
  NMDeff_PC <- NMDeff_PC[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2","bins","cancer_type","sample")]
  colnames(NMDeff_PC) <- c("Endogenous","ASE","bins","cancer_type","sample")
  PC_bins_cancer_type <- data.frame(table(NMDeff_PC$bins,NMDeff_PC$cancer_type))
  colnames(PC_bins_cancer_type) <- c(paste0("CNA_PC_bins"),"cancer_type","samples")
  # % of samples
  PC_bins_cancer_type <- PC_bins_cancer_type %>%
    dplyr::group_by(cancer_type) %>%
    dplyr::mutate(samples_percentage = samples / sum(samples) * 100)
  # Calculate the order of the cancer_type based on the percentage in (-65.5,-31.1] category
  cancer_type_order <- PC_bins_cancer_type %>%
    dplyr::filter(CNA_PC_bins == names(table(NMDeff_PC$bins))[bin_order]) %>%
    dplyr::arrange(desc(samples_percentage)) %>%
    dplyr::pull(cancer_type)
  PC_bins_cancer_type$cancer_type <- factor(PC_bins_cancer_type$cancer_type, levels = as.character(cancer_type_order))
  # Group variable reversed if necessary
  if (PC %in% c(3,52)) {
    PC_bins_cancer_type$CNA_PC_bins <- fct_rev(PC_bins_cancer_type$CNA_PC_bins)
  }
  
  # Barplot
  list_results[["PC_incidence_barplot"]] <- PC_bins_cancer_type
  png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/CNA_PC",PC,"_by_bins_and_cancer_type.png"), 
      width = 4500, height = 5000, res = 300)
  p <- ggplot(PC_bins_cancer_type, aes(y = round(samples_percentage,2), x = cancer_type, fill = CNA_PC_bins)) + 
    geom_bar(position="stack", stat="identity")  + theme_classic() + #ylim(c(0,100)) +
    labs(title = paste0("CNA-PC",PC," cancer type incidence"), x = "", y = "% of individuals", fill = paste0("Group")) + 
    scale_fill_brewer(labels = c("Low", "Mid", "High"), palette = "Set2", direction = -1) +
    theme(plot.title = element_text(hjust = 0.5, size = 55, margin = margin(t = 20, b = 20)),
          axis.title.x = element_text(color="black", size=50),
          axis.title.y = element_text(color="black", size=55),
          axis.text.y = element_text(color="black", size=50),
          axis.text.x = element_text(color="black", size=30, angle = 90, hjust=0.5, vjust=0.5),
          legend.text = element_text(size = 50),         # Increase legend label font size
          legend.title = element_text(size = 55),
          legend.title.align = 0.5,
          legend.position = "top") +
          guides(fill = guide_legend(override.aes = list(size = 14), nrow = 1))
  print(p)
  dev.off()

  # 2) Boxplot of NMDeff for each CNA-PC bin
  NMDeff_PC_stacked <- stack(NMDeff_PC)
  NMDeff_PC_stacked$bins <- rep(NMDeff_PC$bins,2)
  NMDeff_PC_stacked$cancer_type <- rep(NMDeff_PC$cancer_type,2)
  NMDeff_PC_stacked$sample <- rep(NMDeff_PC$sample,2)
  NMDeff_PC_stacked <- na.omit(NMDeff_PC_stacked)
  colnames(NMDeff_PC_stacked) <- c("NMDeff","NMD_method","bins","cancer_type","sample")
  NMDeff_PC_stacked <- na.omit(NMDeff_PC_stacked)
  combinations <- combn(names(table(NMDeff_PC_stacked$bins)), 2, simplify = FALSE)
  if (PC == 86) {
    NMDeff_PC_stacked$bins <- fct_rev(NMDeff_PC_stacked$bins)
  }

  list_results[["NMDeff_boxplot"]] <- NMDeff_PC_stacked
  png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/NMDeff_by_CNA_PC",PC,"_bins_pancancer.png"), 
    width = 5000, height = 3500, res = 300)
  p <- NMDeff_PC_stacked %>% 
        #filter(cancer_type == "COAD") %>%
        group_by(NMD_method,bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
            geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("NMD efficiency") +
            ggtitle(paste0("iNMDeff for the CNA-PC",PC," groups")) +
            geom_text(fontface = "bold", size = 10) +
            scale_x_discrete(labels = c("High", "Mid", "Low")) +
            facet_wrap( ~ NMD_method, scales = "free") + coord_cartesian(ylim = c(-2,2)) + 
            geom_boxplot(width=0.3, color="black", alpha=0.2) +
            theme_bw(base_size = 30) +
            theme(plot.title = element_text(hjust = 0.5, size = 55, margin = margin(t = 20, b = 20)),,
                  axis.text.x = element_text(size = 50, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 55),
                  axis.text.y = element_text(size = 50),
                  axis.title.y = element_text(size = 55),
                  strip.text = element_text(size = 50),
                  legend.position = "none") +
            stat_compare_means(comparisons = combinations, size = 10,
                              label.y = c(0.25,0.5,0.75),
                              label = "p.format", method = "wilcox.test", hide.ns = TRUE)
  print(p)
  dev.off()
  
  # Which is the CNV state for the SMG5/7 genes for the most weighted samples
  TCGA_CNV_filt <- TCGA_CNV_genes_updated[TCGA_CNV_genes_updated$gene_name %in% c("SMG5","SMG7"),grep("TCGA",colnames(TCGA_CNV_genes_updated))]
  TCGA_CNV_filt <- data.frame(t(TCGA_CNV_filt))
  rownames(TCGA_CNV_filt) <- gsub("\\.","-",rownames(TCGA_CNV_filt))
  TCGA_CNV_filt$CNV <- TCGA_CNV_filt[,1] + TCGA_CNV_filt[,2]
  table(TCGA_CNV_filt$CNV)
  sample_NMD_efficiencies_TCGA_filt <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_filt, by.x="sample", by.y="row.names")

  # 3.2.1) Visualization of arm-level CNV in chr[X] for each Percentile of samples from PC[X]

  # Obtain genes from Chr[X]
  TCGA_CNV_genes_filt <- TCGA_CNV_genes_updated[TCGA_CNV_genes_updated$chr %in% chr,]
  # Obtain average CNV state per position pancancer samples
  CNV_genes_mean_samples <- rowMeans(TCGA_CNV_genes_filt[,grep("TCGA",colnames(TCGA_CNV_genes_filt))])
  TCGA_CNV_genes_filt$CNV_samples_average <- data.frame(CNV_genes_mean_samples)$CNV_genes_mean_samples
  # Distribution of CNV-PC signature?
  png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/chr",chr,"_PC",PC,"_TCGA_CNV_density_plot.png"), 
      width = 4000, height = 5000, res = 300)
  p <- ggplot(  TCGA_CNV_PC, aes(x = eval(parse(text=paste0("Dim.",PC))))) +
        geom_density() + theme_classic() +
        labs(title = "Density Plot", x = "Values", y = "Density") +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.title.x = element_text(color="black", size=45),
              axis.title.y = element_text(color="black", size=45),
              axis.text.y = element_text(color="black", size=40),
              axis.text.x = element_text(color="black", size=40))
  print(p)
  dev.off()

  # Average CNV state
  colnames(TCGA_CNV_genes_filt) <- gsub("\\.","-",colnames(TCGA_CNV_genes_filt))
  count <- 0
  samples_PC_bins <- list()
  if (PC %in% c(3,52)) {
    bins_levels <- levels(bins)
  } else if (PC == 86) {
    bins_levels <- rev(levels(bins))
  }
  for (bin in bins_levels) {
    count <- count + 1
    TCGA_CNV_PC_filt <- TCGA_CNV_PC[TCGA_CNV_PC$bins  %in% bin,]
    samples_PC <- rownames(TCGA_CNV_PC_filt)
    CNV_genes_mean_samples <- rowMeans(TCGA_CNV_genes_filt[,colnames(TCGA_CNV_genes_filt) %in% samples_PC])
    TCGA_CNV_genes_filt[,paste0("PC_bin",count)] <- NA
    TCGA_CNV_genes_filt[,paste0("PC_bin",count)] <- data.frame(CNV_genes_mean_samples)$CNV_genes_mean_samples
    # Take top 5 samples for each bin
    if (count == 1) {
      samples_PC <- samples_PC[1:5]
    } else {
      samples_PC <- samples_PC[(length(samples_PC)-4):(length(samples_PC))]
    }
    samples_PC_bins[[bin]] <- samples_PC
  }
  colnames(TCGA_CNV_genes_filt) <- gsub("\\.","-",colnames(TCGA_CNV_genes_filt))

  #for (type in c("average","top_5")) {
  for (type in c("average")) {

    if (type == "top_5") {
      char_grep_1 <- paste0("TCGA|PC")
      # Filter by Top 5 samples from bins
      top5_samples_1 <- which(colnames(TCGA_CNV_genes_filt) %in% samples_PC_bins[[1]])
      colnames(TCGA_CNV_genes_filt)[top5_samples_1] <- paste0(colnames(TCGA_CNV_genes_filt)[top5_samples_1],"_bin1")
      top5_samples_2 <- which(colnames(TCGA_CNV_genes_filt) %in% samples_PC_bins[[2]])
      colnames(TCGA_CNV_genes_filt)[top5_samples_2] <- paste0(colnames(TCGA_CNV_genes_filt)[top5_samples_2],"_bin2")
      top5_samples_3 <- which(colnames(TCGA_CNV_genes_filt) %in% samples_PC_bins[[3]])
      colnames(TCGA_CNV_genes_filt)[top5_samples_3] <- paste0(colnames(TCGA_CNV_genes_filt)[top5_samples_3],"_bin3")
      # top5_samples_4 <- which(colnames(TCGA_CNV_genes_filt) %in% samples_PC_bins[[4]])
      # colnames(TCGA_CNV_genes_filt)[top5_samples_4] <- paste0(colnames(TCGA_CNV_genes_filt)[top5_samples_4],"_bin4")
      top5_samples <- c(top5_samples_1,top5_samples_2,top5_samples_3)
    } else if (type == "average") {
      char_grep_1 <- "TCGA"
      top5_samples <- c()
    }

    cols <- grep(char_grep_1,colnames(TCGA_CNV_genes_filt))
    cols <- cols[!cols %in% top5_samples]
    TCGA_CNV_genes_filt2 <- TCGA_CNV_genes_filt[,-c(cols)]
    # Order by Genome Location
    TCGA_CNV_genes_filt2$genome_location <- as.numeric(gsub(".*:(.*)-.*","\\1",TCGA_CNV_genes_filt2$genome_location))
    TCGA_CNV_genes_filt2 <- TCGA_CNV_genes_filt2[order(TCGA_CNV_genes_filt2$genome_location),]
    # Prepare DF for plotting
    #cols <- grep("TCGA-|CNV_samples_average|bin",colnames(TCGA_CNV_genes_filt2))
    cols <- grep("TCGA-|bin",colnames(TCGA_CNV_genes_filt2))
    n_samples <- length(cols)
    final_df <- stack(TCGA_CNV_genes_filt2[,cols])
    final_df$ensembl_gene_id <- as.character(rep(TCGA_CNV_genes_filt2$ensembl_gene_id,n_samples))
    final_df$gene_name <- rep(TCGA_CNV_genes_filt2$gene_name,n_samples)
    final_df$genome_location <- rep(TCGA_CNV_genes_filt2$genome_location,n_samples)
    final_df$genome_location_str <- as.character(rep(TCGA_CNV_genes_filt2$genome_location_str,n_samples))
    final_df$chr_arm <- as.character(rep(TCGA_CNV_genes_filt2$chr_arm,n_samples))
    # Renaming
    #final_df$ind <- gsub("PC_",paste0("CNA-PC",PC," "),final_df$ind)
    final_df$ind <- as.character(final_df$ind)
    final_df$ind <- ifelse(final_df$ind == "PC_bin1","High",final_df$ind)
    final_df$ind <- ifelse(final_df$ind == "PC_bin2","Mid",final_df$ind)
    final_df$ind <- ifelse(final_df$ind == "PC_bin3","Low",final_df$ind)
    # Remove NA
    final_df <- final_df[!is.na(final_df$chr_arm),]
    # Factor variable
    final_df$ind <- factor(final_df$ind, levels = c("High", "Mid", "Low"))
    SMG5_genome_location <- unique(final_df[final_df$gene_name %in% c("NRAS","RBM8A","SMG5","SMG7"),"genome_location"])
    list_results[["CNA_signature"]] <- final_df
    p <- ggplot(data = final_df, aes(x = genome_location, y = values, 
                group = ind, color = ind)) +
      geom_point(size = 3) + ggtitle(paste0("CNA-PC",PC," stratified by groups of individuals")) +
      geom_line() +  
      scale_color_brewer(palette="Set2") +
      xlab(paste0("Chr ",chr," Genome Location")) +
      ylab(paste0("CNA state")) +
      theme_classic(base_size = 30) +
      geom_hline(yintercept=0) +
      #geom_vline(xintercept=SMG5_genome_location, size = 2, color = "red") +
      #facet_grid(ind ~ ., scales = "fixed") +
      theme(plot.title = element_text(hjust = 0.5, size = 55),
            axis.title.x = element_text(color="black", size=65),
            axis.text.x = element_blank(),
            strip.text = element_text(size = 55),
            axis.title.y = element_text(color="black", size=65),
            axis.text.y = element_text(color="black", size=55),
            legend.position='top',
            legend.title = element_text(size = 55),
            legend.text = element_text(size = 50)) + 
      scale_x_discrete(guide=guide_axis(n.dodge = 3)) +
      guides(color = guide_legend(override.aes = list(size = 16), nrow = 1)) +
      labs(color = "Group")
    png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/chr",chr,"_PC",PC,"_TCGA_CNV_",type,".png"), 
                  width = 5100, height = 5000, res = 300)
    print(p)
    dev.off()
  }
  return(list_results)
}

PCs_CNV_genes_heatmap <- function(TCGA_CNV_PCA_var_updated, CNV_type, num_PCs, scale_PCs) {
  # Scale
  if (scale_PCs == "yes") {
    TCGA_CNV_PCA_var_filt <- data.frame(scale(TCGA_CNV_PCA_var_updated, center = TRUE, scale = TRUE))
  } else {
    TCGA_CNV_PCA_var_filt <- TCGA_CNV_PCA_var_updated
  }
  # Filter by num of PCs to show
  if (num_PCs != "all") {
    TCGA_CNV_PCA_var_filt <- TCGA_CNV_PCA_var_filt[,1:num_PCs]
  }
  # Filter by CNV type
  TCGA_CNV_PCA_var_df <- TCGA_CNV_PCA_var_filt[TCGA_CNV_PCA_var_filt$CNV_type %in% CNV_type,]
  # Remove NA
  TCGA_CNV_PCA_var_df <- TCGA_CNV_PCA_var_df[!is.na(TCGA_CNV_PCA_var_df$chr),]
  TCGA_CNV_PCA_var_df$chr <- ifelse(TCGA_CNV_PCA_var_df$chr == "X","23",TCGA_CNV_PCA_var_df$chr)
  TCGA_CNV_PCA_var_df$chr <- ifelse(TCGA_CNV_PCA_var_df$chr == "Y","24",TCGA_CNV_PCA_var_df$chr)
  # Order by Genome Location (inside each chromosome)
  TCGA_CNV_PCA_var_df <- order_genome_location_by_chr(TCGA_CNV_PCA_var_df)
  # Genes of interest
  TCGA_CNV_PCA_var_df$genes_of_interest <- NA
  TCGA_CNV_PCA_var_df[TCGA_CNV_PCA_var_df$gene_symbols %in% c("MYC","ERBB2","CCNE1","EGFR"),"genes_of_interest"] <- TCGA_CNV_PCA_var_df[TCGA_CNV_PCA_var_df$gene_symbols %in% c("MYC","ERBB2","CCNE1","EGFR"),"gene_symbols"]
  cols <- grep("Dim",colnames(TCGA_CNV_PCA_var_df))
  # Change names of PCs
  colnames(TCGA_CNV_PCA_var_df)[grep("Dim.",colnames(TCGA_CNV_PCA_var_df))] <- paste0("PC ",1:86) 
  # Color chr
  n_colors <- length(unique(TCGA_CNV_PCA_var_df$chr))
  set.seed(123)
  color_vector <- sample(rgb(runif(100), runif(100), runif(100)), n_colors)
  #color_vector <- viridis(n_colors)
  colors <- setNames(color_vector, unique(TCGA_CNV_PCA_var_df$chr))
  TCGA_CNV_PCA_var_df$color <- mapvalues(
    TCGA_CNV_PCA_var_df$chr,
    from = names(colors),
    to = colors
  )
  color_df <- data.frame(Chr = as.character(TCGA_CNV_PCA_var_df$chr))
  rownames(color_df) <- rownames(TCGA_CNV_PCA_var_df)
  color_list <- list(Chr = colors)

  # Plot
  # heatmap_color <- colorRampPalette(c("grey", "black"))(100)
  # heatmap_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/CNV_",CNV_type,"_genes_",num_PCs,"PCs_heatmap_scaled_",scale_PCs,".png")
  # png(heatmap_plot, width = 5500, height = 3500, res = 300)
  # p <- heatmap(as.matrix((TCGA_CNV_PCA_var_df[,cols])),  
  #             Rowv = NA, 
  #             title = "Somatic pan-cancer CNA signatures",
  #             col = heatmap_color,
  #             Colv = NA,
  #             labRow = TCGA_CNV_PCA_var_df$genes_of_interest,
  #             RowSideColors = TCGA_CNV_PCA_var_df$color,
  #             cexRow=2,
  #             cexCol=2)
  # print(p)
  # dev.off()

  # Plot
  heatmap_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/CNV_",CNV_type,"_genes_",num_PCs,"PCs_heatmap_scaled_",scale_PCs,".png")
  png(heatmap_plot, width = 6000, height = 4000, res = 300)
  pheatmap(as.matrix((TCGA_CNV_PCA_var_df[,cols])),
          #main = "Somatic pan-cancer CNA signatures",
          fontsize_col = 16,
          # fontsize_row = 2,
          annotation_row = color_df,
          annotation_colors = color_list,
          show_colnames = TRUE,
          show_rownames = FALSE,
          cellwidth=14, 
          cellheight=0.035,
          cluster_rows = FALSE,
          na_col = "#DDDDDD",
          cluster_cols = FALSE,
          fontsize = 23)
  dev.off()

  return(list(TCGA_CNV_PCA_var_df = TCGA_CNV_PCA_var_df, color_df = color_df, color_list = color_list))

}

order_genome_location_by_chr <- function(genes_chr) {
  genes_chr$genome_location_chr <- as.numeric(gsub(".*:(.*)-.*","\\1",genes_chr$genome_location))
  genes_chr$chr <- as.numeric(genes_chr$chr)
  genes_chr <- genes_chr[order(genes_chr$chr, genes_chr$genome_location_chr, decreasing = TRUE),]
  return(genes_chr)
}

gene_coCNV_frequency <- function(gene, chr, TCGA_CNV_genes_updated_tmp, CNV_level_type, CNV_type) {
  # Subset of genes from Chr [X]
  TCGA_CNV_genes_filt <- TCGA_CNV_genes_updated_tmp[TCGA_CNV_genes_updated_tmp$chr %in% chr,]
  # Obtain ensembl gene id
  ensembl_gene_id <- as.character(TCGA_CNV_genes_filt[TCGA_CNV_genes_filt$gene_symbols %in% gene,"ensembl_gene_id"])
  # Prepare dataframe
  TCGA_CNV_genes_filt2 <- TCGA_CNV_genes_filt[,grep("TCGA",colnames(TCGA_CNV_genes_filt))]
  TCGA_CNV_genes_filt2 <- data.frame(t(TCGA_CNV_genes_filt2[,-c(1:2)]))
  colnames(TCGA_CNV_genes_filt2) <- as.character(TCGA_CNV_genes_filt$ensembl_gene_id)
  if (CNV_type == "amp") {
    coCNV_char_exp_1 <- "sum(.x > 0 & TCGA_CNV_genes_filt2[[ensembl_gene_id]] > 0) / nrow(TCGA_CNV_genes_filt2)"
    coCNV_char_exp_2 <- "sum(.x > 1.5 & TCGA_CNV_genes_filt2[[ensembl_gene_id]] > 1.5) / nrow(TCGA_CNV_genes_filt2)"
    CNV_char_exp_1 <- "sum(.x > 0) / nrow(TCGA_CNV_genes_filt2)"
    CNV_char_exp_2 <- "sum(.x > 1.5) / nrow(TCGA_CNV_genes_filt2)"
  } else if (CNV_type == "del") {
    coCNV_char_exp_1 <- "sum(.x < 0 & TCGA_CNV_genes_filt2[[ensembl_gene_id]] < 0) / nrow(TCGA_CNV_genes_filt2)"
    coCNV_char_exp_2 <- "sum(.x < -1.5 & TCGA_CNV_genes_filt2[[ensembl_gene_id]] < -1.5) / nrow(TCGA_CNV_genes_filt2)"
    CNV_char_exp_1 <- "sum(.x < 0) / nrow(TCGA_CNV_genes_filt2)"
    CNV_char_exp_2 <- "sum(.x < -1.5) / nrow(TCGA_CNV_genes_filt2)"
  }
  
  # CoCNV frequency
  coCNV_freq_1 <- TCGA_CNV_genes_filt2 %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::summarise(across(everything(), ~ eval(rlang::parse_expr(coCNV_char_exp_1)), na.rm = TRUE)) %>%
    unlist()
  coCNV_freq_2 <- TCGA_CNV_genes_filt2 %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::summarise(across(everything(), ~ eval(rlang::parse_expr(coCNV_char_exp_2)), na.rm = TRUE)) %>%
    unlist()
  # single CNV frequency
  CNV_freq_1 <- TCGA_CNV_genes_filt2 %>%
    dplyr::summarise(across(everything(), ~ eval(rlang::parse_expr(CNV_char_exp_1)), na.rm = TRUE)) %>%
    unlist()
  CNV_freq_2 <- TCGA_CNV_genes_filt2 %>%
    dplyr::summarise(across(everything(), ~ eval(rlang::parse_expr(CNV_char_exp_2)), na.rm = TRUE)) %>%
    unlist()
  # Return
  coCNV_freq_df <- data.frame(coCNV_freq_1 = coCNV_freq_1, coCNV_freq_2 = coCNV_freq_2)
  CNV_freq_df <- data.frame(CNV_freq_1 = CNV_freq_1, CNV_freq_2 = CNV_freq_2)
  final_CNV_freq_df <- merge(CNV_freq_df,coCNV_freq_df, by = "row.names", all.x = TRUE)
  rownames(final_CNV_freq_df) <- final_CNV_freq_df$Row.names
  final_CNV_freq_df$Row.names <- NULL
  return(final_CNV_freq_df)
}

NMDeff_average_by_gene <- function(sample_NMD_efficiencies_TCGA, chr, TCGA_CNV_genes_updated_tmp, CNV_type, gene) {

  NMDeff_estimate <- function(CNV_values, NMD_method, CNV_type) {
    # CNV type threshold
    if (CNV_type == "amp") {
      filter1 <- paste0("(CNV_values > 0)")
      filter2 <- paste0("(CNV_values >= 1.5)")
    } else if (CNV_type == "del") {
      filter1 <- paste0("(CNV_values < 0)")
      filter2 <- paste0("(CNV_values <= -1.5)")
    }
    # Average NMD efficiency of samples with CNV in the gene
    NMDeff_1 <- median(final_df[eval(parse(text=filter1)), NMD_method], na.rm = TRUE)
    NMDeff_2 <- median(final_df[eval(parse(text=filter2)), NMD_method], na.rm = TRUE)
    #NMDeff <- data.frame(NMDeff_1 = NMDeff_1, NMDeff_2 = NMDeff_2)
    return(c(NMDeff_1,NMDeff_2))
  }

  # Subset of genes from Chr [X]
  TCGA_CNV_genes_filt <- TCGA_CNV_genes_updated_tmp[TCGA_CNV_genes_updated_tmp$chr %in% chr,]
  # Obtain ensembl gene id
  ensembl_gene_id <- as.character(TCGA_CNV_genes_filt[TCGA_CNV_genes_filt$gene_symbols %in% gene,"ensembl_gene_id"])
  # Prepare dataframe
  TCGA_CNV_genes_filt2 <- TCGA_CNV_genes_filt[,grep("TCGA",colnames(TCGA_CNV_genes_filt))]
  TCGA_CNV_genes_filt2 <- data.frame(t(TCGA_CNV_genes_filt2[,-c(1:2)]))
  colnames(TCGA_CNV_genes_filt2) <- as.character(TCGA_CNV_genes_filt$ensembl_gene_id)
  # Add NMDeff per sample
  final_df <- merge(TCGA_CNV_genes_filt2, sample_NMD_efficiencies_TCGA[,c("sample","endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")], 
                by.x = "row.names", by.y = "sample", all.x = TRUE)
  # Calculate average NMDeff of samples with the CNV[X] for each gene
  NMDeff_endogenous <- final_df %>%
    dplyr::summarise(across(starts_with("ENSG"), ~ NMDeff_estimate(.x,"endogenous_NMD_Consensus",CNV_type)))
  NMDeff_ASE <- final_df %>%
    dplyr::summarise(across(starts_with("ENSG"), ~ NMDeff_estimate(.x,"ASE_PTC_NMD_triggering_0.2",CNV_type)))
  NMDeff_endogenous <- data.frame(NMDeff_endogenous = t(NMDeff_endogenous))
  NMDeff_ASE <- data.frame(NMDeff_ASE = t(NMDeff_ASE))
  NMDeff_genes <- merge(NMDeff_endogenous, NMDeff_ASE, by = "row.names", all.x = TRUE)
  colnames(NMDeff_genes)[1] <- "ensembl_gene_id"
  colnames(NMDeff_genes) <- gsub("\\.","_",colnames(NMDeff_genes))
  return(NMDeff_genes)
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

#################################### LIBRARIES ######################################
#################################### LIBRARIES ######################################
#################################### LIBRARIES ######################################

library("GWASTools")
library("ggplot2")
library("biomaRt")
#library("Karyoplotter")
library("ggrepel")
library("qvalue")
library("ggbreak")
library("cowplot")
library("viridis")
library("ggpubr")
library("RColorBrewer")
library("stringr")
library("tidyverse")
# reverse factor
library("forcats")
library("gplots")
library("plyr")
library("dplyr")
library("pheatmap")

#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################

# 1) Data
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
# 1.1) NMD efficiencies
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

# 1.2) NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1.3) ENSEMBL gene id and Gene Symbol conversion
ensembl_v88_gene_transcript_genesymbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                                      sep = "\t", header = TRUE)
                      
# 1.4) ENSEMBL gene annotation
# ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$type == "gene",]
ensembl_v88_gtf_filt$genome_location <- gsub("chr","",paste0(ensembl_v88_gtf_filt$seqnames,":",ensembl_v88_gtf_filt$start,"-",ensembl_v88_gtf_filt$end))
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[,c("gene_name","gene_id","genome_location")]
ensembl_v88_gtf_filt$gene_id <- gsub("(.*)\\..*","\\1",ensembl_v88_gtf_filt$gene_id)

# 1.5) TCGA samples metadata
TCGA_samples_metadata <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_patients_info.csv")

# 1.6) Cytogenetic locations of genes
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# anno_gene <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","band", "gene_biotype"),mart=ensembl, useCache = FALSE)
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/biomaRt_ensembl_table.txt")
# write.table(anno_gene, file = output_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
anno_gene <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1.7) somatic CNV PCA data 
TCGA_cancer_names_path <- "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
# 1.7.1) Individuals weights
TCGA_cancer <- "pancancer"
scale <- TRUE
center <- TRUE
alpha <- "3e-04"
num_PCs <- "100"
tryCatch({
    TCGA_CNV_PCA_ind <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer,"_sparse_PCA_ind_",alpha,"_robust_no_num_PCs_",num_PCs,".txt"),
                                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    rownames(TCGA_CNV_PCA_ind) <- gsub("\\.","-",rownames(TCGA_CNV_PCA_ind))
    # Scale PCs
    #TCGA_CNV_PCA_ind <- data.frame(scale(TCGA_CNV_PCA_ind, scale = scale, center = center))
    # Remove PCs with 0
    cols <- colnames(TCGA_CNV_PCA_ind)[which( colSums(TCGA_CNV_PCA_ind) != 0 )]
    TCGA_CNV_PCA_ind <- TCGA_CNV_PCA_ind[,cols]
    print("PCA dimensions --> ")
    print(dim(TCGA_CNV_PCA_ind))
    colnames(TCGA_CNV_PCA_ind) <- paste0("Dim.",1:ncol(TCGA_CNV_PCA_ind))
}, error = function(e){
    print(e)
    }
)
# Merge
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_PCA_ind, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 1.7.2) Genes weights
tryCatch({
    TCGA_CNV_PCA_var <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer,"_sparse_PCA_var_",alpha,"_robust_no_num_PCs_",num_PCs,".txt"),
                                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    # Scale PCs
    #TCGA_CNV_PCA_var <- data.frame(scale(TCGA_CNV_PCA_var, scale = scale, center = center))
    # Remove PCs with 0
    cols <- colnames(TCGA_CNV_PCA_var)[which( colSums(TCGA_CNV_PCA_var) != 0 )]
    TCGA_CNV_PCA_var <- TCGA_CNV_PCA_var[,cols]
    print("PCA dimensions --> ")
    print(dim(TCGA_CNV_PCA_var))
}, error = function(e){
    print(e)
    }
)
# Update genes
CNV_type <- "amp"
TCGA_CNV_PCA_var_amp <- TCGA_CNV_PCA_var[grep(CNV_type,rownames(TCGA_CNV_PCA_var)),]
rownames(TCGA_CNV_PCA_var_amp) <- gsub(paste0("_",CNV_type),"",rownames(TCGA_CNV_PCA_var_amp))
TCGA_CNV_PCA_var_amp_updated <- genes_update(TCGA_CNV_PCA_var_amp,ensembl_v88_gene_transcript_genesymbol)
TCGA_CNV_PCA_var_amp_updated$CNV_type <- CNV_type
CNV_type <- "del"
TCGA_CNV_PCA_var_del <- TCGA_CNV_PCA_var[grep(CNV_type,rownames(TCGA_CNV_PCA_var)),]
rownames(TCGA_CNV_PCA_var_del) <- gsub(paste0("_",CNV_type),"",rownames(TCGA_CNV_PCA_var_del))
TCGA_CNV_PCA_var_del_updated <- genes_update(TCGA_CNV_PCA_var_del,ensembl_v88_gene_transcript_genesymbol)
TCGA_CNV_PCA_var_del_updated$CNV_type <- CNV_type
TCGA_CNV_PCA_var_updated <- rbind(TCGA_CNV_PCA_var_amp_updated,TCGA_CNV_PCA_var_del_updated)
colnames(TCGA_CNV_PCA_var_updated) <- gsub("Dim-","Dim\\.",colnames(TCGA_CNV_PCA_var_updated))

# 1.8) TCGA CNV genes data 
# 1.8.1) Arm-level
TCGA_CNV_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/TCGA_pancancer_CNV.txt")
TCGA_CNV_genes <- read.table(file = TCGA_CNV_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(dim(TCGA_CNV_genes))
# Update Genes
TCGA_CNV_genes_updated <- genes_update(TCGA_CNV_genes,ensembl_v88_gene_transcript_genesymbol)

# 1.8.2) Focal-level
# 1.1) CNV files
# TCGA_CNV_genes_focal <- c()
# for (TCGA_cancer in TCGA_cancers) {
#   print(TCGA_cancer)
#   tryCatch({
#     CNV_file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_",gsub("TCGA-","",TCGA_cancer),"-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/focal_data_by_genes.txt")
#     if (TCGA_cancer == "TCGA-SKCM") {
#       CNV_file_path <- gsub("-TP","-TM",paste0(CNV_file_path))
#     }
#     # CNV
#     CNV_file <- read.table(file = CNV_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#     colnames(CNV_file) <- substr(colnames(CNV_file),1,12)
#     rownames(CNV_file) <- CNV_file$Gene.Symbol
#     CNV_file <- CNV_file[,grep("TCGA*",colnames(CNV_file))]
#     if (length(TCGA_CNV_genes_focal) == 0) {
#       TCGA_CNV_genes_focal <- CNV_file
#     } else {
#       TCGA_CNV_genes_focal <- merge(TCGA_CNV_genes_focal, CNV_file, by.x= "row.names", by.y = "row.names", all.x = TRUE)
#       rownames(TCGA_CNV_genes_focal) <- TCGA_CNV_genes_focal$Row.names
#       TCGA_CNV_genes_focal$Row.names <- NULL
#     }
#   }, error = function(e){
#       print(e)
#       }
#   )
# }

TCGA_CNV_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/TCGA_pancancer_CNV_focal.txt")
# write.table(TCGA_CNV_genes_focal, file = TCGA_CNV_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
TCGA_CNV_genes_focal <- read.table(file = TCGA_CNV_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(dim(TCGA_CNV_genes_focal))
# Update Genes
TCGA_CNV_genes_focal_updated <- genes_update(TCGA_CNV_genes_focal,ensembl_v88_gene_transcript_genesymbol)

# 1.9) Gene-level somatic mut/CNV associations
# 1.9.1) ASE
somatic_mut_NMDeff_associations_ASE <- somatic_mut_NMDeff_associations(NMD_method_VAF = "ASE_0.2", ensembl_v88_gtf = ensembl_v88_gtf_filt, PCs = "0")
# 1.9.2) NMDeff
somatic_mut_NMDeff_associations_endogenous <- somatic_mut_NMDeff_associations(NMD_method_VAF = "endogenous", ensembl_v88_gtf = ensembl_v88_gtf_filt, PCs = "1")

# 1.10) Classification of CL based on arm-level CNV by --> Cohen Sharir Nat 2021
# output_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/Cohen_Sharir_Nat_2021_DepMap_CNV_CL.csv"
# DepMap_CNV_arm_level_type <- read.csv(file = output_path, header = TRUE, sep = ";")

# 1.11) Co-dependency scores
# Calculate pairs of correlations between all genes and UPF1, SMG1
# Do this for every Chr 1q CL classification separately

# # 1.11.1) CRISPR KO Data (Achilles - DepMap)
# CRISPR_achilles_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/Achilles_gene_effect.csv"
# CRISPR_achilles_file <- read.csv(file = CRISPR_achilles_path, header = TRUE)
# CRISPR_achilles_KO <- merge(CRISPR_achilles_file,DepMap_CNV_arm_level_type[,c("DepMap_ID","X1q")], by = "DepMap_ID", all.x = TRUE)
# colnames(CRISPR_achilles_KO) <- gsub("(.*)\\.\\..*","\\1",colnames(CRISPR_achilles_KO))
# NMD_KO_codependency_scores_all <- list()
# for (CL_1q in c(-1,0,1)) {
#     print(CL_1q)
#     CRISPR_achilles_KO_tmp <- CRISPR_achilles_KO[CRISPR_achilles_KO$X1q %in% CL_1q,]
#     CRISPR_achilles_KO_tmp$X1q <- NULL
#     NMD_KO_codependency_scores <- c()
#     NMD_factors_index <- grep("UPF1|^SMG1", colnames(CRISPR_achilles_KO_tmp))
#     for (i in 2:ncol(CRISPR_achilles_KO_tmp)) {
#         if (i %in% NMD_factors_index) {next}
#         CRISPR_achilles_filt <- CRISPR_achilles_KO_tmp[,c(i,NMD_factors_index)]
#         # Co-dependency score
#         cor_result <- cor(CRISPR_achilles_filt, use = "pairwise.complete.obs", method = "pearson")
#         cor_result <- data.frame(cor_result[1,-1, drop = FALSE])
#         colnames(cor_result) <- paste0(colnames(cor_result),"_KO_codependency_score")
#         cor_result$gene_symbol <- rownames(cor_result)
#         # Save
#         if (length(NMD_KO_codependency_scores) == 0) {
#             NMD_KO_codependency_scores <- cor_result
#         } else {
#             NMD_KO_codependency_scores <- rbind(NMD_KO_codependency_scores,cor_result)
#         }
#     }
#     NMD_KO_codependency_scores$CL_chr1q <- CL_1q
#     NMD_KO_codependency_scores$gene_symbol <- rownames(NMD_KO_codependency_scores)
#     NMD_KO_codependency_scores_all[[as.character(CL_1q)]] <- NMD_KO_codependency_scores
# }
# NMD_KO_codependency_scores_by_chr1_CNV <- do.call(rbind,NMD_KO_codependency_scores_all)
# rownames(NMD_KO_codependency_scores_by_chr1_CNV) <- NULL

# # # 1.11.2) siRNA Data (DepMap)
# siRNA_DepMap_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/D2_combined_gene_dep_scores.csv"
# siRNA_DepMap_file <- read.table(file = siRNA_DepMap_path, header = TRUE, sep = ",", check.names=FALSE)
# rownames(siRNA_DepMap_file) <- gsub("(.*) \\(.*","\\1",siRNA_DepMap_file[,1])
# siRNA_DepMap_file[,1] <- NULL
# CCLE_ID <- colnames(siRNA_DepMap_file)
# siRNA_DepMap_file <- data.frame(t(siRNA_DepMap_file))
# siRNA_DepMap_file$CCLE_ID <- CCLE_ID
# siRNA_DepMap_KD <- merge(siRNA_DepMap_file,DepMap_CNV_arm_level_type[,c("CCLE_ID","X1q")], by.x = "CCLE_ID", by.y = "CCLE_ID", all.x = TRUE)
# NMD_KD_codependency_scores_all <- list()
# for (CL_1q in c(-1,0,1)) {
#   print(CL_1q)
#   siRNA_DepMap_KD_tmp <- siRNA_DepMap_KD[siRNA_DepMap_KD$X1q %in% CL_1q,]
#   siRNA_DepMap_KD_tmp[,c("X1q")] <- NULL
#   NMD_KD_codependency_scores <- c()
#   NMD_factors_index <- grep("UPF1$|^SMG1$", colnames(siRNA_DepMap_KD_tmp))
#   for (i in 2:ncol(siRNA_DepMap_KD_tmp)) {
#     if (i %in% NMD_factors_index) {next}
#     siRNA_DepMap_file_filt <- siRNA_DepMap_KD_tmp[,c(i,NMD_factors_index)]
#     # Co-dependency score
#     cor_result <- cor(siRNA_DepMap_file_filt, use = "pairwise.complete.obs", method = "pearson")
#     cor_result <- data.frame(cor_result[1,-1, drop = FALSE])
#     colnames(cor_result) <- paste0(colnames(cor_result),"_KD_codependency_score")
#     cor_result$gene_symbol <- rownames(cor_result)
#     # Save
#     if (length(NMD_KD_codependency_scores) == 0) {
#       NMD_KD_codependency_scores <- cor_result
#     } else {
#       NMD_KD_codependency_scores <- rbind(NMD_KD_codependency_scores,cor_result)
#     }
#   }
#     NMD_KD_codependency_scores$CL_chr1q <- CL_1q
#     NMD_KD_codependency_scores$gene_symbol <- rownames(NMD_KD_codependency_scores)
#     NMD_KD_codependency_scores_all[[as.character(CL_1q)]] <- NMD_KD_codependency_scores
# }
# NMD_KD_codependency_scores_by_chr1_CNV <- do.call(rbind,NMD_KD_codependency_scores_all)
# rownames(NMD_KD_codependency_scores_by_chr1_CNV) <- NULL

# # # 1.10.3) Merge
# NMD_KO_codependency_scores_by_chr1_CNV <- NMD_KO_codependency_scores_by_chr1_CNV %>%
#                       rowwise() %>%
#                       mutate(mean_KO_codependency_score = mean(c_across(c('SMG1_KO_codependency_score', 'UPF1_KO_codependency_score')), na.rm=TRUE))
# NMD_KD_codependency_scores_by_chr1_CNV <- NMD_KD_codependency_scores_by_chr1_CNV %>%
#                       rowwise() %>%
#                       mutate(mean_KD_codependency_score = mean(c_across(c('SMG1_KD_codependency_score', 'UPF1_KD_codependency_score')), na.rm=TRUE))
# NMD_codependency_scores <- merge(NMD_KO_codependency_scores_by_chr1_CNV, NMD_KD_codependency_scores_by_chr1_CNV, 
#                             by = c("gene_symbol","CL_chr1q") , all = TRUE)
# NMD_codependency_scores <- NMD_codependency_scores[!duplicated(NMD_codependency_scores),]
# # Add Ensembl Gene ID
# ensembl_gene_symbol <- na.omit(unique(ensembl_v88_gene_transcript_genesymbol[,c("gene_id","gene_name")]))
# NMD_codependency_scores <- merge(NMD_codependency_scores,ensembl_gene_symbol, by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE)

# # Write
file_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/DepMap/NMD_codependeny_scores.txt"
# write.table(NMD_codependency_scores, file = file_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
NMD_codependency_scores <- read.table(file = file_path, header = TRUE, sep = "\t")

#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################
#################################### ANALYSIS ######################################

# 2) Manhattan plot - fine mapping at each Chromosome
# Gene-level associations of focal somatic CNV and NMDeff, adjusted within each Chr and stratified by arm
for (chr in c(1:22,"X","Y")) {
  print(chr)
  manhattan_fine_mapping(somatic_mut_NMDeff_associations = somatic_mut_NMDeff_associations_ASE, chr = chr, 
                          NMD_method_VAF = "ASE_0.2", adjust = "yes", type = "qval")
  manhattan_fine_mapping(somatic_mut_NMDeff_associations = somatic_mut_NMDeff_associations_endogenous, chr = chr, 
                          NMD_method_VAF = "endogenous", adjust = "yes", type = "qval")
  manhattan_fine_mapping(somatic_mut_NMDeff_associations = somatic_mut_NMDeff_associations_ASE, chr = chr, NMD_method_VAF = "ASE_0.2")
  manhattan_fine_mapping(somatic_mut_NMDeff_associations = somatic_mut_NMDeff_associations_endogenous, chr = chr, NMD_method_VAF = "endogenous")
}

# 3) Identification of affected Chromosome associated to each significant PC 

# 3.1) Heatmap of gene-weights vs PCs
TCGA_CNV_PCA_var_list <- PCs_CNV_genes_heatmap(TCGA_CNV_PCA_var_updated, CNV_type = "amp", num_PCs = "all", scale_PCs = "no")
# Save results
# write.table(TCGA_CNV_PCA_var_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_CNV_PCA_var_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9B.RData")
TCGA_CNV_PCA_var_list <- PCs_CNV_genes_heatmap(TCGA_CNV_PCA_var_updated, CNV_type = "del", num_PCs = "all", scale_PCs = "no")
# Save results
# write.table(TCGA_CNV_PCA_var_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9C.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_CNV_PCA_var_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9C.RData")

# 3.2) CNV state (Y) for each gene vs chromosome location (X)

# 3.2.1) Pancancer significant PCs

########### CNA-PC3 ###########

# Dim 3 and 86 --> Amp in Chr 1q (more focal in 86) ## SMG5, SMG7, RBM8A as candidates
CNV_state_in_chr_for_PC_list <- CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 3, chr = "1", 
                          sample_NMD_efficiencies_TCGA, quantiles = c("0%","30%","40%","100%"))
# Save results
write.table(CNV_state_in_chr_for_PC_list[["CNA_signature"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["CNA_signature"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3B.RData")
write.table(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3E.RData")
# Supp
write.table(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11C.RData")
write.table(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13A.RData")

#### Manual estimates for the manuscript ####
# Do not scale iNMDeff
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_PCA_ind, by.x = "sample", by.y = "row.names", all.x = TRUE)
CNV_state_in_chr_for_PC_list <- CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, chr = "1", 
                          # sample_NMD_efficiencies_TCGA, quantiles = c("0%","30%","40%","100%"), PC = 3)
                          sample_NMD_efficiencies_TCGA, quantiles = c("0%","65%","95%","100%"), PC = 86)
df <- CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]]
aggregate(NMDeff ~ NMD_method + bins, data = df, median)
# df <- df[df$NMD_method == "ASE" & df$bins %in% c("[31.1000,65.5000)","[-0.0100,13.2000)"),]
# wilcox.test(NMDeff ~ bins, data = df, paired = FALSE)
# Linear model for associations between groups
df <- merge(df,sample_NMD_efficiencies_TCGA, by.x = "sample", by.y = "sample", all.x = TRUE)
df[,c("endogenous_purity","age","sample_lib_size")] <- scale(df[,c("endogenous_purity","age","sample_lib_size")])
df <- df[df$NMD_method == "Endogenous",]
glm_char <- "endogenous_NMD_Consensus ~"
CNV_PC <- "bins"
glm_model <- paste0("glm(",glm_char," endogenous_purity + as.factor(cancer_subtype) + 
                          as.factor(sex) + age +  sample_lib_size + relevel(as.factor(bins), ref = \"[-0.0100,13.2000)\"), data = df, 
                          family = \"gaussian\", na.action = na.exclude)")
glm_res <- eval(parse(text=glm_model))
glm_res <- summary(glm_res)
number <- glm_res$coefficients[grep("bins",rownames(glm_res$coefficients)),"Estimate"]
format(number, scientific = FALSE, digits = 3)
# % of samples with less NMDeff than median by cancer type ## Dim3
df2 <- df[df$NMD_method == "Endogenous" & df$cancer_type == "LUAD",]
median_NMDeff <- median(df2[df2$bins == "[-0.0100,13.2000)","NMDeff"]) # group Low
samples_with_amp <- sum(df2$bins == "[31.1000,65.5000)")
samples_with_amp/nrow(df2)
length(which(df2$bins == "[31.1000,65.5000)" & df2$NMDeff <= median_NMDeff))/samples_with_amp # Group High
# % of samples with less NMDeff than median by cancer type ## Dim86
df2 <- df[df$NMD_method == "Endogenous" & df$cancer_type == "OV",]
median_NMDeff <- median(df2[df2$bins == "(-0.01,0.0129]","NMDeff"]) # group Low
samples_with_amp <- sum(df2$bins == "(1.55,3.07]")
samples_with_amp/nrow(df2)
length(which(df2$bins == "(1.55,3.07]" & df2$NMDeff <= median_NMDeff))/samples_with_amp # Group High

########### CNA-PC86 ###########
#CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 86, chr = "1", sample_NMD_efficiencies_TCGA, quantiles = c("0%","5%","35%","100%"))
CNV_state_in_chr_for_PC_list <- CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 86, chr = "1", sample_NMD_efficiencies_TCGA, quantiles = c("0%","65%","95%","100%"))
# Save results
write.table(CNV_state_in_chr_for_PC_list[["CNA_signature"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["CNA_signature"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11A.RData")
write.table(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11D.RData")
write.table(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13B.RData")

########### CNA-PC52 ###########
# Dim 52 --> Amp in Chr 2p
CNV_state_in_chr_for_PC_list <- CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 52, chr = "2", sample_NMD_efficiencies_TCGA, quantiles = c("0%","5%","30%","100%"))
# Save results
write.table(CNV_state_in_chr_for_PC_list[["CNA_signature"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["CNA_signature"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11B.RData")
write.table(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["NMDeff_boxplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/SuppFig11E.RData")
write.table(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(CNV_state_in_chr_for_PC_list[["PC_incidence_barplot"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/SuppFig13C.RData")


# 3.2.2) By cancer type PCs

# # Dim 16 --> Chr 2 ?? Check this on READ cancer
# CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 16, chr = "2", sample_NMD_efficiencies_TCGA)


# # 3.2.3) Other PCs/Chrs
# # Chr6 has a peak in FAM46A (sig <0.1 in both ASE and Endogenous even after qval correction) in manhattan plots. But different directions of NMDeff
# # I think this is in Dim 14 (arm del of q), which is significant in FDR <50% in TCGA-GBM (this gene is overexpressed in GBM https://pubmed.ncbi.nlm.nih.gov/33370626/)
# CNV_state_in_chr_for_PC(TCGA_CNV_PCA_ind, TCGA_CNV_genes_updated, PC = 14, chr = "6", sample_NMD_efficiencies_TCGA)
# # Chr7 has a perfect peak in q

# 4) gene-NMDeff vs CNV for each chromosome

# Using Arm-level CNV and foca-level CNV data

sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

for (CNV_level_type in c("arm","focal")) {

  print(paste0("CNV level type --> ",CNV_level_type))

  if ( CNV_level_type == "arm") {
    TCGA_CNV_genes_updated_tmp <- TCGA_CNV_genes_updated
  } else if (CNV_level_type == "focal") {
    TCGA_CNV_genes_updated_tmp <- TCGA_CNV_genes_focal_updated
  }
  chr <- "2"
  CNV_type <- "amp"
  #gene <- "SMG5" #SF3B1
  gene <- "SF3B1"
  if (chr == "1") {
    genes_label <- c("SMG5","SMG7","RBM8A","SF3B4","INTS3","SMG7","MDM4","MDM2")
  } else if (chr == "2") {
    genes_label <- c("IL1RN","SF3B1","FARSB","SNRNP200","NOP58","EIF2B4","NBAS","CWC22")
  }
  
  # Frequency of co-amplification of selected genes (e.g. SMG5)
  gene_coCNV_freq <- gene_coCNV_frequency(gene = gene, chr = chr, TCGA_CNV_genes_updated_tmp = TCGA_CNV_genes_updated_tmp, 
                                            CNV_level_type = CNV_level_type, CNV_type = CNV_type)
  # Average NMDeff of samples with CNV in [X] gene, for all genes along a [X] chromosome
  NMDeff_genes_CNV <- NMDeff_average_by_gene(sample_NMD_efficiencies_TCGA, chr = chr, gene = gene,
                                          TCGA_CNV_genes_updated_tmp = TCGA_CNV_genes_updated_tmp, CNV_type = CNV_type)
  # Merge
  NMDeff_and_coCNV_genes <- merge(NMDeff_genes_CNV, gene_coCNV_freq, by.x = "ensembl_gene_id", by.y = "row.names", all.x = TRUE)
  # Add Genome location and order
  NMDeff_and_coCNV_genes <- merge(NMDeff_and_coCNV_genes , TCGA_CNV_genes_updated[,-grep("TCGA",colnames(TCGA_CNV_genes_updated))],
                                  by = "ensembl_gene_id", all.x = TRUE)
  # Order
  # Order by Genome Location
  NMDeff_and_coCNV_genes$genome_location <- as.numeric(gsub(".*:(.*)-.*","\\1",NMDeff_and_coCNV_genes$genome_location))
  NMDeff_and_coCNV_genes <- NMDeff_and_coCNV_genes[order(NMDeff_and_coCNV_genes$genome_location),]
  # Colors
  CNV_color <- brewer.pal(n = 2, name = "Paired")
  coAmp_color <- brewer.pal(n = 8, name = "OrRd")

  # Final plot
  plots_list <- list()
  for (NMD_method in c("ASE","endogenous")) {
    NMDeff_and_coCNV_genes[,paste0("NMDeff_",NMD_method,"_1")] <- NMDeff_and_coCNV_genes[,paste0("NMDeff_",NMD_method,"_1")] - abs(max(NMDeff_and_coCNV_genes[,paste0("NMDeff_",NMD_method,"_1")]))
    NMDeff_and_coCNV_genes$top_hits <- NA
    NMDeff_and_coCNV_genes <- NMDeff_and_coCNV_genes[order(-NMDeff_and_coCNV_genes[,paste0("NMDeff_",NMD_method,"_1")]),]
    NMDeff_and_coCNV_genes[1:30,"top_hits"] <- NMDeff_and_coCNV_genes[1:30,"gene_symbols"]
    png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/",NMD_method,"_NMDeff_and_freq_coAmp_",CNV_level_type,"_",CNV_type,"_",gene,"_chr_",chr,".png"), width = 4500, height = 3500, res = 300)
    p <- ggplot(data = NMDeff_and_coCNV_genes, aes(x = genome_location, y = eval(parse(text=paste0("NMDeff_",NMD_method,"_1"))), color = as.numeric(coCNV_freq_1))) +
      geom_point(size = 3) +
      xlab(paste0("Chr ",chr," Genome Location")) + ggtitle(paste0(gsub("endogenous","ETG",NMD_method),paste0(" iNMDeff across Chr ",chr," ",CNV_type))) +
      ylab(paste0("median NMD efficiency    CNA ",CNV_type," frequency")) +
      theme_classic() +
      geom_bar(aes(x = genome_location, y = CNV_freq_1), stat = "identity", color = CNV_color[1]) +
      geom_bar(aes(x = genome_location, y = CNV_freq_2), stat = "identity", color = CNV_color[2]) +
      geom_label_repel(data = subset(NMDeff_and_coCNV_genes, gene_symbols %in% genes_label), #%in% c(NMD_genes$gene_symbol)), 
                      aes(label = gene_symbols), color = "black", size = 8, arrow = arrow(length = unit(0.01, 'npc')),
                          point.padding = 0.2, nudge_x = .05, nudge_y = .05) +
      # geom_text_repel(data = subset(NMDeff_and_coCNV_genes, gene_symbols %in% c(NMDeff_and_coCNV_genes$top_hits)), 
      #                 aes(label = gene_symbols), color = "red", size = 3, max.overlaps = nrow(NMDeff_and_coCNV_genes),
      #                     point.padding = 0.2, nudge_x = -.05, nudge_y = -.05) +
      geom_hline(yintercept= 0, size = 2.5, color = "black") +
      scale_colour_gradient(name = paste0("Frequency of co-",CNV_type," with ",gene), 
              low = coAmp_color[3], high = coAmp_color[8],
              labels = scales::percent_format(scale = 100),
              guide = guide_colorbar(barheight = 1, barwidth = 17)) + #limits = c(0.1,0.4)) +
      theme(plot.title = element_text(hjust = 0.5, size = 45),
            axis.title.x = element_text(color="black", size=45),
            axis.text.x = element_blank(),
            axis.title.y = element_text(color="black", size=33,hjust=0.1),
            axis.text.y = element_text(color="black", size=35),
            legend.title=element_text(size=35),
            legend.text=element_text(size=20),
            legend.position = "bottom")
    plots_list[[NMD_method]] <- p
    print(p)
    dev.off()
    # Save results
    if (CNV_level_type == "focal" & CNV_type == "amp") {
      if (chr == "1") {
        if (NMD_method == "endogenous") {
          write.table(NMDeff_and_coCNV_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3D.txt", 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
          saveRDS(NMDeff_and_coCNV_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3D.RData")
        } else if (NMD_method == "ASE") {
          write.table(NMDeff_and_coCNV_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14A.txt", 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
          saveRDS(NMDeff_and_coCNV_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14A.RData")     
        }
      } else if (chr == "2") {
        if (NMD_method == "endogenous") {
          write.table(NMDeff_and_coCNV_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14B.txt", 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
          saveRDS(NMDeff_and_coCNV_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14B.RData")
        } else if (NMD_method == "ASE") {
          write.table(NMDeff_and_coCNV_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14C.txt", 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
          saveRDS(NMDeff_and_coCNV_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/SuppFig14C.RData")     
        }      
      }
    }
  }

  # Plot with NMD Codependency scores

  # Merge 
  NMD_codependency_and_coCNV_genes <- merge(gene_coCNV_freq,NMD_codependency_scores, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
  # Add Genome location and order
  NMD_codependency_and_coCNV_genes <- merge(NMD_codependency_and_coCNV_genes , TCGA_CNV_genes_updated[,-grep("TCGA",colnames(TCGA_CNV_genes_updated))],
                                  by.x ="Row.names", by.y = "ensembl_gene_id", all.x = TRUE)
  # Order by Genome Location
  NMD_codependency_and_coCNV_genes$genome_location <- as.numeric(gsub(".*:(.*)-.*","\\1",NMD_codependency_and_coCNV_genes$genome_location))
  NMD_codependency_and_coCNV_genes <- NMD_codependency_and_coCNV_genes[order(NMD_codependency_and_coCNV_genes$genome_location),]

  sig_band <- c("q21.1","q21.2","q21.3","q22","q23.1","q23.2","q23.3","q24.1","q24.2","q24.3","q25.1","q25.2",
                  "q25.3","q31.1","q31.2","q31.3")
  # sig_band <- c("q21.1","q21.2","q21.3","q22")

  var <- c( paste0(c("SMG1","UPF1"),c("_KO")) , paste0(c("SMG1","UPF1"),c("_KD")) )
  for (NMD_gene in var) {
    print(NMD_gene)
    # Remove negative values
    col <- paste0(NMD_gene,"_codependency_score")
    #NMD_codependency_and_coCNV_genes[which(NMD_codependency_and_coCNV_genes[,col] < 0),col] <- NA
    #NMD_codependency_and_coCNV_genes[,col] <- -NMD_codependency_and_coCNV_genes[,col]
    # Threshold
    threshold <- as.numeric(quantile(NMD_codependency_and_coCNV_genes[,col],seq(0,1,0.02),na.rm = TRUE)[2])
    # Min Ylim
    min_ylim <- round(min(NMD_codependency_and_coCNV_genes[,col],na.rm=TRUE)-0.02,2)
    NMD_codependency_and_coCNV_genes <- NMD_codependency_and_coCNV_genes[!is.na(NMD_codependency_and_coCNV_genes$CL_chr1q),]
    NMD_codependency_and_coCNV_genes <- NMD_codependency_and_coCNV_genes[grep("q",NMD_codependency_and_coCNV_genes$band),]
    # Top Hits
    NMD_codependency_and_coCNV_genes <- NMD_codependency_and_coCNV_genes %>%
          dplyr::group_by(CL_chr1q) %>%
          dplyr::arrange(CL_chr1q,desc(eval(parse(text=col))) ) %>%
          #dplyr::mutate(top_hits = ifelse(row_number() <= 20 & genome_location <= 175010789, gene_symbols, NA)) %>%
          dplyr::mutate(top_hits = ifelse(row_number() <= 20 & band %in% sig_band, gene_symbols, NA)) %>%
          ungroup()

    png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/KO_KD_DepMap/",NMD_gene,"_codependency_score_and_freq_co",CNV_type,"_",gene,".png"), 
                width = 5500, height = 3500, res = 300)
    p <- ggplot(data = NMD_codependency_and_coCNV_genes, aes(x = genome_location, y = eval(parse(text=paste0(NMD_gene,"_codependency_score"))), 
          color = as.numeric(coCNV_freq_1))) +
          facet_grid(. ~ CL_chr1q) +
          geom_point(size = 3) +
          xlab(paste0("Chr ",chr," Genome Location")) +
          ylab(paste0("Codependency score with ",NMD_gene," CNV amp frequency")) +
          theme_classic() +
          # geom_bar(aes(x = genome_location, y = CNV_freq_1), stat = "identity", color = "blue") +
          # geom_bar(aes(x = genome_location, y = CNV_freq_2), stat = "identity", color = "red") +
          geom_label_repel(data = subset(NMD_codependency_and_coCNV_genes, gene_symbols %in% c(NMD_genes$gene_symbol)), 
                          aes(label = gene_symbols), color = "black", size = 5, arrow = arrow(length = unit(0.01, 'npc')),
                              point.padding = 0.2, nudge_x = .05, nudge_y = .05) +
          geom_text_repel(data = subset(NMD_codependency_and_coCNV_genes, gene_symbols %in% c(NMD_codependency_and_coCNV_genes$top_hits)), 
                          aes(label = gene_symbols), color = "red", size = 3, max.overlaps = nrow(NMD_codependency_and_coCNV_genes),
                              point.padding = 0.2, nudge_x = -.05, nudge_y = -.05) +
          geom_hline(yintercept= 0, size = 2.5, color = "black") +
          scale_colour_gradient(name = paste0("Frequency of co-amp with ",gene), low = "yellow", high = "red") + #limits = c(0.1,0.4)) +
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                axis.title.x = element_text(color="black", size=20),
                axis.text.x = element_blank(),
                axis.title.y = element_text(color="black", size=20,hjust=0.1),
                axis.text.y = element_text(color="black", size=20),
                legend.title=element_text(size=12),
                strip.text = element_text(size = 20),
                legend.position = "bottom")
    print(p)
    dev.off()

    NMD_KD_codependency_score_Chr1q_CNV_amp <- NMD_codependency_and_coCNV_genes[NMD_codependency_and_coCNV_genes$CL_chr1q == 1,c("gene_symbol","genome_location",col)]
    NMD_KD_codependency_score_Chr1q_CNV_del <- NMD_codependency_and_coCNV_genes[NMD_codependency_and_coCNV_genes$CL_chr1q == -1,c("gene_symbol","genome_location",col)]
    NMD_KD_codependency_score_Chr1q_CNV_neutral <- NMD_codependency_and_coCNV_genes[NMD_codependency_and_coCNV_genes$CL_chr1q == 0,c("gene_symbol","genome_location",col)]
    colnames(NMD_KD_codependency_score_Chr1q_CNV_amp)[3] <- paste0(col,"_amp")
    colnames(NMD_KD_codependency_score_Chr1q_CNV_del)[3] <- paste0(col,"_del")
    colnames(NMD_KD_codependency_score_Chr1q_CNV_neutral)[3] <- paste0(col,"_neutral")
    NMD_KD_codependency_score_Chr1q_CNV <- merge(NMD_KD_codependency_score_Chr1q_CNV_amp,NMD_KD_codependency_score_Chr1q_CNV_del)
    NMD_KD_codependency_score_Chr1q_CNV <- merge(NMD_KD_codependency_score_Chr1q_CNV,NMD_KD_codependency_score_Chr1q_CNV_neutral)
    # Delta difference
    NMD_KD_codependency_score_Chr1q_CNV$delta <- NMD_KD_codependency_score_Chr1q_CNV[,paste0(col,"_amp")] - NMD_KD_codependency_score_Chr1q_CNV[,paste0(col,"_neutral")]
    NMD_KD_codependency_score_Chr1q_CNV$top_hits <- NA
    NMD_KD_codependency_score_Chr1q_CNV <- NMD_KD_codependency_score_Chr1q_CNV[order(NMD_KD_codependency_score_Chr1q_CNV$delta),]
    NMD_KD_codependency_score_Chr1q_CNV[1:30,"top_hits"] <- as.character(NMD_KD_codependency_score_Chr1q_CNV[1:30,"gene_symbol"])

    p <- ggplot(data = NMD_KD_codependency_score_Chr1q_CNV, aes(x = eval(parse(text=paste0(col,"_amp"))), y = eval(parse(text=paste0(col,"_neutral")))) ) +
            geom_point(size = 3) + 
            scale_x_continuous(limits = c(-0.5, 0.5), breaks= seq(-0.5, 0.5, by=0.1)) + 
            scale_y_continuous(limits = c(-0.5, 0.5), breaks= seq(-0.5, 0.5, by=0.1)) +
            xlab("Chr 1q - CNA amp") + ylab("Chr 1q - CNA neutral") +
            geom_label_repel(data = subset(NMD_KD_codependency_score_Chr1q_CNV, gene_symbol %in% c(NMD_KD_codependency_score_Chr1q_CNV$top_hits)), 
                            aes(label = gene_symbol), color = "black", size = 5, arrow = arrow(length = unit(0.01, 'npc')),
                                point.padding = 0.2, nudge_x = .05, nudge_y = .01) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.title.x = element_text(color="black", size=20, face="bold"),
              axis.text.x = element_text(color="black", size=20),
              axis.title.y = element_text(color="black", size=20, face="bold"),
              axis.text.y = element_text(color="black", size=20),
              legend.title=element_text(size=12),
              strip.text = element_text(size = 20),
              legend.position = "bottom")
    png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/KO_KD_DepMap/",NMD_gene,"_codependency_scores_scatterplot.png"), 
                width = 5500, height = 3500, res = 300)
    print(p)
    dev.off() 

  }
  # plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/NMDeff_and_freq_coAmp_",CNV_level_type,"_",gene,".png")

  # p1 <- plots_list$ASE
  # p2 <- plots_list$endogenous
  
  # png(plot, width = 6000, height = 3500, res = 300)
  # p_final <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
  # print(p_final)
  # dev.off()

}
