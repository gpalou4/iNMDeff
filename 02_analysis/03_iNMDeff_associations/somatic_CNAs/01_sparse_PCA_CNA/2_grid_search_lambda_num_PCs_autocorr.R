library("ggplot2")
library("dplyr")
library("cowplot")

calc_autocorr3 <- function(df, lag) {
  re = unlist(apply(df, 2, function(d) {
    r = acf(d, lag.max = lag, plot = FALSE)
    r = r$acf[2]
    # r_ran = unlist(
    #   lapply(1:10, function(i){
    #     d_ran = sample(d)
    #     r_ran = acf(d_ran, lag.max = 1, plot = FALSE)
    #     return(r_ran$acf[2])
    #   })
    # )
    # zsco = (r - mean(r_ran)) / sd(r_ran)
    # return(zsco)
    r
  })
  )
  return(re)
} 
# Calculate p-val for the Z-score
# df_re = data.frame(sig = names(a), zsco= a, stringsAsFactors = F)
# df_re$pvalue = 2*pnorm(-abs(df_re$zsco))
# df_re$padj = p.adjust(df_re$pvalue, n = dim(gene_symbol_chr)[1])

add_genome_location <- function(TCGA_CNV_PCA) {
  # Add gene symbols / ensembl gene id to TCGA_CNV_genes
  genes_symbols <- unlist(lapply(strsplit(rownames(TCGA_CNV_PCA),"\\|"), function(x){
    x[1]
  }))
  ensembl_gene_id <- unlist(lapply(strsplit(rownames(TCGA_CNV_PCA),"\\|"), function(x){
    x[2]
  }))
  ensembl_gene_id <- gsub("(.*)\\..*","\\1",ensembl_gene_id)
  TCGA_CNV_PCA$gene_symbols <- gsub("_amp|_del","",genes_symbols)
  TCGA_CNV_PCA$ensembl_gene_id <- ensembl_gene_id
  # Update TCGA_CNV_genes ENSEMBLE GENE IDs (90%)
  # Genes with no ENSEMBL gene ID are UNIQUE, so let's add it from our table
  missing_rows <- which(is.na(TCGA_CNV_PCA$ensembl_gene_id))
  matching_genes <- ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_name %in% TCGA_CNV_PCA[,"gene_symbols"],c("gene_id","gene_name")]
  matching_genes <- matching_genes[!duplicated(matching_genes),]
  TCGA_CNV_genes_1 <- TCGA_CNV_PCA[missing_rows,]
  TCGA_CNV_genes_1$gene_symbol_CNV <- rownames(TCGA_CNV_genes_1)
  TCGA_CNV_genes_1 <- merge(TCGA_CNV_genes_1,matching_genes, by.x = "gene_symbols", by.y = "gene_name", all.x = TRUE)
  TCGA_CNV_genes_1$ensembl_gene_id <- TCGA_CNV_genes_1$gene_id
  TCGA_CNV_genes_1$gene_id <- NULL
  TCGA_CNV_genes_1 <- TCGA_CNV_genes_1[-which(duplicated(TCGA_CNV_genes_1$gene_symbol_CNV)),]
  TCGA_CNV_genes_2 <- TCGA_CNV_PCA[-missing_rows,]
  TCGA_CNV_genes_2$gene_symbol_CNV <- rownames(TCGA_CNV_genes_2)
  TCGA_CNV_genes_2 <- merge(TCGA_CNV_genes_2,matching_genes, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_2$gene_name <- NULL
  TCGA_CNV_genes_updated <- rbind(TCGA_CNV_genes_1,TCGA_CNV_genes_2)
  # Add Chromosome location
  TCGA_CNV_genes_filt <- merge(TCGA_CNV_genes_updated,ensembl_v88_gtf, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = TRUE)
  TCGA_CNV_genes_filt <- TCGA_CNV_genes_filt[!is.na(TCGA_CNV_genes_filt$genome_location),]
  TCGA_CNV_genes_filt$genome_location <- as.numeric(gsub(".*:(.*)-.*","\\1",TCGA_CNV_genes_filt$genome_location))
  # TCGA_CNV_genes_filt <- merge(TCGA_CNV_genes_filt, anno_gene[,c("ensembl_gene_id","band")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  # TCGA_CNV_genes_filt$chr_arm <- factor(substr(TCGA_CNV_genes_filt$band,1,1))
  # TCGA_CNV_genes_filt <- TCGA_CNV_genes_filt[order(TCGA_CNV_genes_filt$genome_location),]
  # Order
  TCGA_CNV_genes_filt <- TCGA_CNV_genes_filt[order(TCGA_CNV_genes_filt$seqnames,TCGA_CNV_genes_filt$genome_location),]
  return(TCGA_CNV_genes_filt)
}

autocorr_PCs_CNV <- function(PCA_type, view = "var", alpha = NULL, scale_PCs = NULL, num_PCs = NULL) {

  print("PARAMETERS --> ")
  print(paste(PCA_type,view,alpha,scale_PCs, sep = "____"))

  # 1) Open sparse PCA CNV
  if (PCA_type == "normal") {
    input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_PCA_",view,".txt")
    alpha_char <- ""
  } else if (PCA_type == "sparse") {
    input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_",view,"_",alpha,"_robust_no_num_PCs_",num_PCs,".txt")
    alpha_char <- paste0(alpha,"_")
  }
  TCGA_CNV_PCA <- read.table(file = input_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  # Remove PCs with 0
  cols <- colnames(TCGA_CNV_PCA)[which( colSums(TCGA_CNV_PCA) != 0 )]
  TCGA_CNV_PCA <- TCGA_CNV_PCA[,cols]
  print("PCA dimensions --> ")
  print(dim(TCGA_CNV_PCA))
  nPCs <- ncol(TCGA_CNV_PCA)
  # Scale PCs
  if ( !is.null(scale_PCs) ) {
    TCGA_CNV_PCA <- data.frame(scale(TCGA_CNV_PCA, scale = TRUE, center = TRUE))
    scale_char <- "scaled"
  } else {
    scale_char <- "not_scaled"
  }
  # if ( !is.null(num_PCs) ) {
  #     if (nPCs < num_PCs) {
  #       return(NULL)
  #     } else {
  #       TCGA_CNV_PCA <- TCGA_CNV_PCA[,1:num_PCs]
  #       nPCs <- num_PCs
  #     }
  # }
  # print("PCA dimensions after --> ")
  # print(dim(TCGA_CNV_PCA))
  # if (ncol(TCGA_CNV_PCA) == 0) {
  #   return(NULL)
  # }

  # 2) Obtain Chromosome location and order the PCs by position
  TCGA_CNV_genes <- add_genome_location(TCGA_CNV_PCA)
  cols <- grep("Dim",colnames(TCGA_CNV_genes))

  # 3) Split by AMP and DEL 
  TCGA_CNV_genes_amp <- TCGA_CNV_genes[grep("_amp",TCGA_CNV_genes$gene_symbol_CNV),]
  TCGA_CNV_genes_del <- TCGA_CNV_genes[grep("_del",TCGA_CNV_genes$gene_symbol_CNV),]
  
  # 4) Test autocorrelations for each PC
  PCs_autocorr_amp <- calc_autocorr3(TCGA_CNV_genes_amp[,cols], lag = 1)
  PCs_autocorr_del <- calc_autocorr3(TCGA_CNV_genes_del[,cols], lag = 1)
  # Merge
  PCs_autocorr_amp_df <- data.frame(CNV_amp = PCs_autocorr_amp)
  PCs_autocorr_del_df <- data.frame(CNV_del = PCs_autocorr_del)
  PCs_autocorr_df <- merge(PCs_autocorr_amp_df,PCs_autocorr_del_df, by = "row.names", all.x = TRUE)
  colnames(PCs_autocorr_df)[1] <- "PC"
  PCs_autocorr_df$num_PC <- as.numeric(gsub("Dim\\.(.*)","\\1",PCs_autocorr_df$PC))
  PCs_autocorr_df <- PCs_autocorr_df[order(PCs_autocorr_df$num_PC),]
  # 5) Histogram of correlations
  PCs_autocorr_stack <- stack(PCs_autocorr_df[,-4])
  PCs_autocorr_stack$PC <- rep(PCs_autocorr_df$num_PC,2)

  # Barplot of correlations
  # output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PCA_CNV/",PCA_type,"_PCA_",view,"_",scale_char,"_",alpha_char,"autocorrelations.png")
  # png(output_path, width = 5000, height = 3000, res = 300)
  # p <- ggplot(data = PCs_autocorr_stack, aes(x = as.factor(PC), y = values, fill = ind)) +
  #   geom_bar(stat = "identity", position = "dodge") +
  #   #geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
  #   ylab("Autocorrelation") + xlab("PCs") + ggtitle(paste0("Autocorrelation of ",PCA_type," PCA for - ",view)) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 20),
  #         axis.text.x = element_text(color="black", size=10, angle = 75),
  #         axis.title.y = element_text(color="black", size=13, face="bold"),
  #         axis.title.x = element_text(color="black", size=13, face="bold"))
  # print(p)
  # dev.off()

  # Histogram
  p <- ggplot(data = PCs_autocorr_stack, aes(x = values, fill = ind)) +
    geom_histogram(col = "black", binwidth = 0.02, bins = 100) +
    #geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
    ylab("Frequency") + xlab("Autocorrelation") + ggtitle(paste0("Autocorrelation of ",PCA_type," PCA \n Alpha --> ",alpha, " - #PCs --> ", nPCs)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.x = element_text(color="black", size=10, angle = 75),
          axis.title.y = element_text(color="black", size=13, face="bold"),
          axis.title.x = element_text(color="black", size=13, face="bold"))
  autocorrelation_PCs <- na.omit(PCs_autocorr_stack$values)
  percentiles <- quantile(autocorrelation_PCs,probs = seq(0, 1, 0.01))
  #sum(autocorrelation_PCs >= percentiles["10%"]) / length(autocorrelation_PCs)
  return(list(plot = p, percentiles = percentiles, eff_num_PCs = nPCs))
}

# 1) ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
ensembl_v88_gtf <- ensembl_v88_gtf[ensembl_v88_gtf$type == "gene",]
ensembl_v88_gtf$genome_location <- gsub("chr","",paste0(ensembl_v88_gtf$seqnames,":",ensembl_v88_gtf$start,"-",ensembl_v88_gtf$end))
ensembl_v88_gtf <- ensembl_v88_gtf[,c("gene_name","gene_id","genome_location","seqnames","start")]
ensembl_v88_gtf$gene_id <- gsub("(.*)\\..*","\\1",ensembl_v88_gtf$gene_id)

# 1.3) ENSEMBL gene id and Gene Symbol conversion
ensembl_v88_gene_transcript_genesymbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                                      sep = "\t", header = TRUE)

# 2) Autocorrelations

percentile <- "1"
PCs <- c(50,100,150,200,250,300,350,400,450,500)
alphas <- c("1e-04","0.00015","2e-04","3e-04","5e-04","7e-04")
grid_search <- matrix(nrow = length(PCs), ncol = length(alphas))
rownames(grid_search) <- PCs
colnames(grid_search) <- paste0("a_",alphas)
grid_search <- data.frame(grid_search)
effective_num_PCs <- grid_search

for (num_PCs in PCs) {
  print(paste0("Num of PCs --> ",num_PCs))
  row <- rownames(grid_search) %in% num_PCs
  plots_list <- NULL
  plots_list <- lapply(alphas, function(alpha) {
    res <- autocorr_PCs_CNV(PCA_type = "sparse", view = "var", alpha = alpha, scale_PCs = NULL, num_PCs = num_PCs)
    if (! is.null(res) ) {
      col <- colnames(grid_search) %in% gsub("\\-","\\.",paste0("a_",alpha))
      grid_search[row,col] <<- as.numeric(res$percentiles[paste0(percentile,"%")])
      effective_num_PCs[row,col] <<- as.numeric(res$eff_num_PCs)
      res$plot
    }
  })
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PCA_CNV/All_autocorrelations_num_PCs_",num_PCs,"_percentile_",percentile,".png")
  png(output_path, width = 7000, height = 3000, res = 300)
  p <- cowplot::plot_grid(plotlist=plots_list, labels = "AUTO", align = "v", ncol = 3, nrow = 2)
  print(p)
  dev.off()
}

# 3) Plot of grid search
grid_search_stacked <- stack(grid_search)
grid_search_stacked$PCs <- rep(rownames(grid_search),length(alphas))
grid_search_stacked$PCs <- factor(grid_search_stacked$PCs, levels = rownames(grid_search))
effective_num_PCs_stacked <- stack(effective_num_PCs)
grid_search_stacked$eff_num_PCs <- effective_num_PCs_stacked$values

# Grid Search Barlot
library("RColorBrewer")
display.brewer.all()

p <- ggplot(data = grid_search_stacked, aes(x = as.factor(PCs), y = values, fill = ind)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "YlOrBr", direction = -1) +
  labs(fill = "Alpha") +
  ylab(paste0(percentile,"% percentile of autocorrelations")) + xlab("Number of PCs") + 
  ggtitle(paste0("Sparse-PCA Grid Search of Best Parameters")) +
  geom_text(aes(label=eff_num_PCs), position=position_dodge(width=0.9), vjust=-0.25, size = 4, angle=45) +
  theme_classic(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 35),
      axis.text.x = element_text(color="black", size=30),
      axis.text.y = element_text(color="black", size=30),
      axis.title.y = element_text(color="black", size=35),
      axis.title.x = element_text(color="black", size=35))

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PCA_CNV/Grid_search_PCs_alpha_percentile_",percentile,".png")
png(output_path, width = 5000, height = 3000, res = 300)
print(p)
dev.off()

colnames(grid_search_stacked) <- c("correlation","alpha","PCs","eff_num_PCs")
# Save results
write.table(grid_search_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(grid_search_stacked, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/SuppFig9A.RData")


# p1 <- ggplot(data = grid_search_stacked, aes(x = as.factor(PCs), y = values, fill = ind)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   #geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
#   ylab("Percentile 1% autocorrelation") + xlab("PCs") + ggtitle(paste0("Grid Search for #PCs-Lambda")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         axis.text.x = element_text(color="black", size=10, angle = 75),
#         axis.title.y = element_text(color="black", size=13, face="bold"),
#         axis.title.x = element_text(color="black", size=13, face="bold"))

# alphas_PCs <- data.frame(alpha=alphas,max_PCs = c(1000,350,199,108,60,41))
# alphas_PCs$alpha <- factor(alphas_PCs$alpha, levels = alphas_PCs$alpha)

# # Effective number of PCs
# p2 <- ggplot(data = alphas_PCs, aes(x = as.factor(alpha), y = max_PCs)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_text(aes(label=max_PCs)) +
#   #geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
#   ylab("Number of non-zero PCs") + xlab("Alpha") + ggtitle(paste0("Effective num of PCs")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         axis.text.x = element_text(color="black", size=10, angle = 75),
#         axis.title.y = element_text(color="black", size=13, face="bold"),
#         axis.title.x = element_text(color="black", size=13, face="bold"))
# # Plot
# output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/PCA_CNV/Grid_search_PCs_alpha.png")
# png(output_path, width = 5000, height = 3000, res = 300)
# p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1,rel_widths = c(3,1.5))
# print(p)
# dev.off()

# Hola, te quiero amor de mi vida.



#######################


# CNV_100 <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_var_1e-04_robust_no_num_PCs_100.txt")
# CNV_500 <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_var_1e-04_robust_no_num_PCs_1000.txt"
# TCGA_CNV_PCA_100 <- read.table(file = CNV_100, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# # Remove PCs with 0
# cols <- colnames(TCGA_CNV_PCA_100)[which( colSums(TCGA_CNV_PCA_100) != 0 )]
# TCGA_CNV_PCA_100 <- TCGA_CNV_PCA_100[,cols]
# TCGA_CNV_PCA_500 <- read.table(file = CNV_500, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# # Remove PCs with 0
# cols <- colnames(TCGA_CNV_PCA_500)[which( colSums(TCGA_CNV_PCA_500) != 0 )]
# TCGA_CNV_PCA_500 <- TCGA_CNV_PCA_500[,cols]

# dim(TCGA_CNV_PCA_100)
# dim(TCGA_CNV_PCA_500)

# head(TCGA_CNV_PCA_100)
# head(TCGA_CNV_PCA_500)


# colnames(TCGA_CNV_PCA_100)[!colnames(TCGA_CNV_PCA_100) %in% colnames(TCGA_CNV_PCA_500)]

# TCGA_CNV_PCA_final <- merge(TCGA_CNV_PCA_500,TCGA_CNV_PCA_100[,"Dim.75", drop = FALSE], by ="row.names",all.x =TRUE)
# rownames(TCGA_CNV_PCA_final) <- TCGA_CNV_PCA_final$Row.names
# TCGA_CNV_PCA_final$Row.names <- NULL
# df <- cor(TCGA_CNV_PCA_final, method = "pearson", use = "pairwise.complete.obs")
# df



# colnames(TCGA_CNV_PCA_500)[!colnames(TCGA_CNV_PCA_500) %in% colnames(TCGA_CNV_PCA_100)]

# TCGA_CNV_PCA_final <- merge(TCGA_CNV_PCA_100,TCGA_CNV_PCA_500[,"Dim.215", drop = FALSE], by ="row.names",all.x =TRUE)
# rownames(TCGA_CNV_PCA_final) <- TCGA_CNV_PCA_final$Row.names
# TCGA_CNV_PCA_final$Row.names <- NULL
# df <- cor(TCGA_CNV_PCA_final, method = "pearson", use = "pairwise.complete.obs")
# df

