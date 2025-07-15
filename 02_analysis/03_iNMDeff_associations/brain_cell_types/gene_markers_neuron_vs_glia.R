library(ggplot2)
library(dplyr)
library(cowplot)

# 1) Data

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

# 1.1) sample NMD efficiencies TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

# 1.2) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)

# 1.3) Gene-level TPM RNAseq data
# TCGA
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt"
RNAseq_TCGA_TPM_all <- read.table(file = output_path, header = TRUE, sep = "\t", row.names = 1)
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))
# Check Samples with NAs
na_counts <- colSums(is.na(RNAseq_TCGA_TPM_all))
cols_na <- names(na_counts[na_counts > 0])
RNAseq_TCGA_TPM_all <- RNAseq_TCGA_TPM_all[, !colnames(RNAseq_TCGA_TPM_all) %in% cols_na]
# Check Genes with NAs
na_counts <- rowSums(is.na(RNAseq_TCGA_TPM_all))
rows_na <- names(na_counts[na_counts > 0])
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))

# GTEx
GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
RNAseq_GTEx_TPM_all <- read.table(file = GTEx_path, header = TRUE, sep = "\t", row.names = 1, skip = 2)
print("Dimensions -->")
print(dim(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all <- RNAseq_GTEx_TPM_all[-grep("PAR", rownames(RNAseq_GTEx_TPM_all)), ]
rownames(RNAseq_GTEx_TPM_all) <- gsub("(.*)\\..*", "\\1", rownames(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all$Description <- NULL

# 1.4) Ensembl gene ID vs Gene Symbol conversion
ensembl_gene_symbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_gene_symbol <- ensembl_gene_symbol[!duplicated(ensembl_gene_symbol[,1:2]),1:2]

# 1.5) Gene markers for neuron and glia
gene_name <- c("Trpm3","RBFOX3","MAP2","AQP4","SLC1A2","S100B","GFAP","OLIG1","OLIG2","MAG","MOG")
ensembl_gene_id <- c("ENSG00000083067","ENSG00000167281","ENSG00000078018","ENSG00000171885","ENSG00000110436",
                    "ENSG00000160307","ENSG00000131095","ENSG00000184221","ENSG00000205927","ENSG00000105695","ENSG00000204655")
cell_type <- c("neuron","neuron","neuron","astrocyte","astrocyte","astrocyte","astrocyte","oligodendrocyte","oligodendrocyte","oligodendrocyte","oligodendrocyte")
marker_genes <- data.frame(gene_name = gene_name, ensembl_gene_id = ensembl_gene_id, cell_type = cell_type)

# 2) Brain

# 2.1) Stratify TCGA Brain subtypes based on molecular/histological data
### TCGA ###
sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA[sample_NMD_efficiencies_TCGA$cancer_type %in% c("GBM","LGG"),]

### Subtype_Selected ###
tags <- as.character(sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"Subtype_Selected"])
tags <- gsub(".*\\.","",tags)
add <- sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"cancer_type"]
tags <- paste0(add,"_",tags)
sample_NMD_efficiencies_TCGA_filt$cancer_type_strat <- as.character(sample_NMD_efficiencies_TCGA_filt$cancer_type_strat)
sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"cancer_type_strat"] <- tags

### GTEx ###
sample_NMD_efficiencies_GTEx_filt <- sample_NMD_efficiencies_GTEx[grep("Brain|Nerve",sample_NMD_efficiencies_GTEx$tissue),]

# 3) Merge gene expression data

# 3.1) TCGA
RNAseq_TCGA_TPM_all_filt <- RNAseq_TCGA_TPM_all[rownames(RNAseq_TCGA_TPM_all) %in% marker_genes$ensembl_gene_id,]
RNAseq_TCGA_TPM_all_filt <- merge(RNAseq_TCGA_TPM_all_filt,marker_genes, by.x = "row.names", by.y = "ensembl_gene_id")
rownames(RNAseq_TCGA_TPM_all_filt) <- RNAseq_TCGA_TPM_all_filt$gene_name
RNAseq_TCGA_TPM_all_filt[,c("Row.names","gene_name","cell_type")] <- NULL
RNAseq_TCGA_TPM_all_filt <- data.frame(t(RNAseq_TCGA_TPM_all_filt))
rownames(RNAseq_TCGA_TPM_all_filt) <- gsub("\\.","-",rownames(RNAseq_TCGA_TPM_all_filt))
sample_NMD_efficiencies_TCGA_filt <- merge(sample_NMD_efficiencies_TCGA_filt, RNAseq_TCGA_TPM_all_filt, by.x = "sample", by.y = "row.names", all.x = TRUE)
# 3.2) GTEx
RNAseq_GTEx_TPM_all_filt <- RNAseq_GTEx_TPM_all[rownames(RNAseq_GTEx_TPM_all) %in% marker_genes$ensembl_gene_id,]
RNAseq_GTEx_TPM_all_filt <- merge(RNAseq_GTEx_TPM_all_filt,marker_genes, by.x = "row.names", by.y = "ensembl_gene_id")
rownames(RNAseq_GTEx_TPM_all_filt) <- RNAseq_GTEx_TPM_all_filt$gene_name
RNAseq_GTEx_TPM_all_filt[,c("Row.names","gene_name","cell_type")] <- NULL
RNAseq_GTEx_TPM_all_filt <- data.frame(t(RNAseq_GTEx_TPM_all_filt))
rownames(RNAseq_GTEx_TPM_all_filt) <- gsub("\\.","-",rownames(RNAseq_GTEx_TPM_all_filt))
sample_NMD_efficiencies_GTEx_filt <- merge(sample_NMD_efficiencies_GTEx_filt, RNAseq_GTEx_TPM_all_filt, by.x = "sample_full_barcode", by.y = "row.names", all.x = TRUE)

NMDeff_vs_neuron_and_glia_markers <- function(dataset, var_to_plot) {

    # 1) NMDeff and neuron proportion
    # Dataset
    NMD_method_order <- "ETG"
    if (var_to_plot == "NMD Consensus") {
        var_to_plot_char <- "ETG"
    } else if (var_to_plot == "PTC NMD-triggering 0.2") {
        var_to_plot_char <- "ASE" 
    } else if (var_to_plot == "neuron") {
        var_to_plot_char <- "Neuron_Cell_Type"
    }
    
    if (dataset == "TCGA") { 
        #TCGA
        sample_NMD_efficiencies <- sample_NMD_efficiencies_TCGA_filt 
        #RNAseq_TPM_all <- RNAseq_TCGA_TPM_all
        tissue_var <- "cancer_type"
        tissue_var_strat <- "cancer_type_strat"
        cell_types_char <- c("neuron")
    } else if (dataset == "GTEx") {
        #GTEx
        sample_NMD_efficiencies <- sample_NMD_efficiencies_GTEx_filt
        tissue_var <- "tissue"
        tissue_var_strat <- "acronyms"
        cell_types_char <- c("neuron")
    }
    colnames(sample_NMD_efficiencies)[colnames(sample_NMD_efficiencies) %in% c("NMD Consensus","PTC NMD-triggering 0.2")] <- c("ETG","ASE")

    # Order by Endogenous NMDeff median
    if (dataset == "GTEx") {
        NMDeff_medians <- aggregate(eval(parse(text=NMD_method_order)) ~ eval(parse(text=tissue_var_strat)), data = sample_NMD_efficiencies, median)
        NMDeff_medians <- NMDeff_medians[order(NMDeff_medians[,2]),]
        sample_NMD_efficiencies[,tissue_var_strat] <- factor(sample_NMD_efficiencies[,tissue_var_strat], levels = unique(NMDeff_medians[,1]))
    } else if (dataset == "TCGA") {
        NMDeff_medians <- aggregate(eval(parse(text=NMD_method_order)) ~ eval(parse(text=tissue_var_strat)), data = sample_NMD_efficiencies, median)
        NMDeff_medians <- NMDeff_medians[order(NMDeff_medians[,2]),]
        sample_NMD_efficiencies[,tissue_var_strat] <- factor(sample_NMD_efficiencies[,tissue_var_strat], levels = unique(NMDeff_medians[,1]))
    }

    # return(sample_NMD_efficiencies)

    # p <- ggplot(data = sample_NMD_efficiencies, aes(x = eval(parse(text=tissue_var_strat)), 
    #                     y = eval(parse(text=var_to_plot)), fill = factor(eval(parse(text=tissue_var))))) +
    #         geom_boxplot() + coord_flip(ylim = c(0,1)) + #coord_cartesian(ylim = c(-3,3)) +
    #         #facet_grid(. ~ NMD_method , scale = "free") +
    #         ylab("") + ggtitle("") + xlab("") + 
    #         theme_bw() + scale_fill_brewer(palette = "Set3") +
    #         theme(plot.title = element_text(hjust = 0.5, size = 20),
    #                 axis.title.x = element_text(color="black", size=15),
    #                 strip.text = element_text(size = 9),
    #                 axis.title.y = element_text(color="black", size=15),
    #                 axis.text.x = element_text(color="black", size=12),# angle = 90, hjust = 1),
    #                 axis.text.y = element_text(color="black", size=12),
    #                 legend.position='top', legend.text = element_text(size = 10)) +
    #                 guides(fill = guide_legend(title = "", override.aes = list(size = 10)))
    # plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/",dataset,"_",var_to_plot,"_brains.png")  
    # png(plot_path, width = 2500, height = 2500, res = 300)
    # print(p)
    # dev.off()

    # 2) Add NMDeff vs marker gene correlations
    cols <- c(as.character(marker_genes$gene_name),"ETG","ASE")
    #cols <- c("ETG","ASE","neuron")
    # Split the data by cancer_type
    split_data <- split(sample_NMD_efficiencies, sample_NMD_efficiencies[,tissue_var_strat])
    # Compute the correlation for each cancer type
    cor_res_list <- lapply(split_data, function(df) {
        cor(as.matrix(df[,cols]), use = "pairwise.complete.obs", method = "pearson")
    })
    # Convert the list to a named list for clarity
    #names(cor_res_list) <- names(split_data)
    # Extract correlations of interest from each matrix and create a data frame
    cor_dfs <- lapply(names(cor_res_list), function(tissue_type) {
        cor_matrix <- cor_res_list[[tissue_type]]
        
        data.frame(
            tissue_type = tissue_type,
            corr_ETG_vs_neuron = cor_matrix["RBFOX3", "ETG"],
            corr_ASE_vs_neuron = cor_matrix["RBFOX3", "ASE"],
            corr_ETG_vs_glia = cor_matrix["AQP4", "ETG"],
            corr_ASE_vs_glia = cor_matrix["AQP4", "ASE"]
        )
    })
    # Bind all data frames together
    corr_final_df <- do.call(rbind, cor_dfs)
    # Merge
    sample_NMD_efficiencies <- merge(sample_NMD_efficiencies,corr_final_df, by.x = tissue_var_strat, by.y = "tissue_type", all.x = TRUE)
    return(sample_NMD_efficiencies)
}
# TCGA

TCGA_NMDeff_gene_markers <- NMDeff_vs_neuron_and_glia_markers(dataset = "TCGA", var_to_plot = "NMD Consensus") 
GTEx_NMDeff_gene_markers <- NMDeff_vs_neuron_and_glia_markers(dataset = "GTEx", var_to_plot = "NMD Consensus") 

# Remove low sample sizes
TCGA_NMDeff_gene_markers <- TCGA_NMDeff_gene_markers[!TCGA_NMDeff_gene_markers$cancer_type_strat %in%
                            c("LGG_NA","GBM_G-CIMP-high","GBM_G-CIMP-low","LGG_G-CIMP-low"),]

write.table(TCGA_NMDeff_gene_markers, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_NMDeff_gene_markers, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6A.RData")
write.table(GTEx_NMDeff_gene_markers, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_NMDeff_gene_markers, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/SuppFig6B.RData")

# RBFOX3, AQP4 --> Are the best candidates

for (gene in as.character(marker_genes$gene_name)) {

    # gene <- as.character(marker_genes$gene_name[1])
    # gene
    # cor.test(TCGA_NMDeff_gene_markers[,gene],TCGA_NMDeff_gene_markers$ETG)
    # cor.test(GTEx_NMDeff_gene_markers[,gene],GTEx_NMDeff_gene_markers$ETG)

    p <- ggplot(data = GTEx_NMDeff_gene_markers, aes(x = acronyms, 
                        y = log2(eval(parse(text=gene))), fill = factor(acronyms))) +
            geom_boxplot() + coord_flip() + #coord_cartesian(ylim = c(-3,3)) +
            #facet_grid(. ~ NMD_method , scale = "free") +
            ylab("") + ggtitle("") + xlab("") + 
            theme_bw() + scale_fill_brewer(palette = "Set3") +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
                    axis.title.x = element_text(color="black", size=15),
                    strip.text = element_text(size = 9),
                    axis.title.y = element_text(color="black", size=15),
                    axis.text.x = element_text(color="black", size=12),# angle = 90, hjust = 1),
                    axis.text.y = element_text(color="black", size=12),
                    legend.position='top', legend.text = element_text(size = 10)) +
                    guides(fill = guide_legend(title = "", override.aes = list(size = 10)))
    plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/marker_genes/GTEx_NMDeff_gene_",gene,"_brains.png")  
    png(plot_path, width = 2500, height = 2500, res = 300)
    print(p)
    dev.off()

    p <- ggplot(data = TCGA_NMDeff_gene_markers, aes(x = cancer_type_strat, 
                        y = log2(eval(parse(text=gene))), fill = factor(cancer_type))) +
            geom_boxplot() + coord_flip() + #coord_cartesian(ylim = c(-3,3)) +
            #facet_grid(. ~ NMD_method , scale = "free") +
            ylab("") + ggtitle("") + xlab("") + 
            theme_bw() + scale_fill_brewer(palette = "Set3") +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
                    axis.title.x = element_text(color="black", size=15),
                    strip.text = element_text(size = 9),
                    axis.title.y = element_text(color="black", size=15),
                    axis.text.x = element_text(color="black", size=12),# angle = 90, hjust = 1),
                    axis.text.y = element_text(color="black", size=12),
                    legend.position='top', legend.text = element_text(size = 10)) +
                    guides(fill = guide_legend(title = "", override.aes = list(size = 10)))
    plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/marker_genes/TCGA_NMDeff_gene_",gene,"_brains.png")  
    png(plot_path, width = 2500, height = 2500, res = 300)
    print(p)
    dev.off()

}
