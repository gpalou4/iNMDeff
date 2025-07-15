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

# 1.5) GTEx Cellular Deconvolution
GTEx_cell_deconv <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/cellular_deconvolution/Donovan_et_al_2020_Nat_Comm.csv",
                        header = TRUE)
sample_NMD_efficiencies_GTEx <- merge(sample_NMD_efficiencies_GTEx,GTEx_cell_deconv, by.x = "sample_full_barcode", by.y = "Input.Sample", all.x = TRUE)

# 1.6) TCGA Cellular Deconvolution
# Done by J.Lanillos
adata_TCGA_pancancer <- read.csv(file = "/g/strcombio/fsupek_cancer1/jlanillos/UCDeconvolve/tcga/adata_tcga_pancancer.tsv",
                                    header=TRUE, sep = "\t")
# updated by Liza
adata_TCGA_pancancer <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_cell_type_fraction/UCDeconvolve_TCGA_pancancer.tsv",
                                    header=TRUE, sep = "\t")
#"/g/strcombio/fsupek_cancer1/jlanillos/TCIA/TCIA-CellTypeFractionsData.withheader.tsv" CYBERSORT
# Filter Brain
adata_TCGA_pancancer <- adata_TCGA_pancancer[which(adata_TCGA_pancancer$cancer_type %in% c("LGG","GBM","PCPG")),]
# cols <- grep("neur|sample",colnames(adata_TCGA_pancancer))
# adata_TCGA_pancancer_filt <- adata_TCGA_pancancer[,cols]
adata_TCGA_pancancer$cancer_type <- NULL
adata_TCGA_pancancer_filt <- adata_TCGA_pancancer
# Fix sample ID
sample_NMDeff_TCGA_filt <- sample_NMD_efficiencies_TCGA
sample_NMDeff_TCGA_filt$sample <- substr(sample_NMD_efficiencies_TCGA$TCGA_full_barcode,1,15)
# Merge
sample_NMD_efficiencies_TCGA_filt <- merge(sample_NMDeff_TCGA_filt,adata_TCGA_pancancer_filt, by.x = "sample", by.y = "sample", all.x = TRUE)

# 2) Brain

# 2.1) Stratify TCGA Brain subtypes based on molecular/histological data
### TCGA ###
sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA_filt[sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG"),]

### Subtype_Selected ###
tags <- as.character(sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"Subtype_Selected"])
tags <- gsub(".*\\.","",tags)
add <- sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"cancer_type"]
tags <- paste0(add,"_",tags)
sample_NMD_efficiencies_TCGA_filt$cancer_type_strat <- as.character(sample_NMD_efficiencies_TCGA_filt$cancer_type_strat)
sample_NMD_efficiencies_TCGA_filt[which(sample_NMD_efficiencies_TCGA_filt$cancer_type %in% c("GBM","LGG")),"cancer_type_strat"] <- tags

### GBM_MGMT_status,GBM_expression_subclass,LGG_histological_type ###
# sample_NMD_efficiencies_TCGA_stacked <- stack(sample_NMD_efficiencies_TCGA_filt[,c("NMD Consensus","PTC NMD-triggering 0.2")])
# sample_NMD_efficiencies_TCGA_stacked$cancer_type <- rep(sample_NMD_efficiencies_TCGA_filt$cancer_type_strat,2)
# colnames(sample_NMD_efficiencies_TCGA_stacked) <- c("NMD_efficiency","NMD_method","cancer_type")

### GTEx ###
sample_NMD_efficiencies_GTEx_filt <- sample_NMD_efficiencies_GTEx[grep("Brain|Nerve",sample_NMD_efficiencies_GTEx$tissue),]

# 2.2) NMD factors gene expression
NMDeff_vs_cell_type_or_NMD_factors_exp <- function(dataset, cancer_type_char, NMD_methods, selected_genes, cell_types, cell_types_char = NULL) {

    # Dataset
    if (dataset == "TCGA") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_TCGA_filt 
        RNAseq_TPM_all <- RNAseq_TCGA_TPM_all
        tissue_var <- "cancer_type"
        tissue_var_strat <- "cancer_type_strat"
    } else if (dataset == "GTEx") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_GTEx_filt
        RNAseq_TPM_all <- RNAseq_GTEx_TPM_all
        sample_NMD_efficiencies$sample <- NULL
        sample_NMD_efficiencies$sample <- sample_NMD_efficiencies$sample_full_barcode
        tissue_var <- "tissue"
        tissue_var_strat <- "acronyms"
        #sample_NMD_efficiencies$tissue <- as.character("Brain")
    }
    # Ensembl gene id to Gene Symbol
    tmp <- merge(RNAseq_TPM_all,ensembl_gene_symbol, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
    tmp <- tmp[!duplicated(tmp$gene_name),]
    rownames(tmp) <- tmp$gene_name
    tmp[,c("Row.names","gene_name")] <- NULL
    tmp <- as.matrix(tmp)
    #  Use NMD targets
    RNAseq_TPM_all_filt <- t(tmp[rownames(tmp) %in% selected_genes,])
    rownames(RNAseq_TPM_all_filt) <- gsub("\\.","-",rownames(RNAseq_TPM_all_filt))
    # Cancer_type
    if (cancer_type_char != "all") {
        if (dataset == "TCGA") {
            sample_NMD_efficiencies_filt <- sample_NMD_efficiencies %>%
                filter(cancer_type == cancer_type_char)
        } else if (dataset == "GTEx") {
            sample_NMD_efficiencies_filt <- sample_NMD_efficiencies %>%
                filter(acronyms == cancer_type_char)         
        }
    } else {
        sample_NMD_efficiencies_filt <- sample_NMD_efficiencies
    }
    # Merge RNA-seq
    NMDeff_samples <- merge(sample_NMD_efficiencies_filt,RNAseq_TPM_all_filt, by.x = "sample", by.y = "row.names", all.x = TRUE)
    # Sort samples by NMDeff
    NMDeff_samples <- NMDeff_samples[order(NMDeff_samples[,"NMD Consensus"]),]
    NMDeff_median <- median(NMDeff_samples[,"NMD Consensus"],na.rm = TRUE)
    NMDeff_samples[which(NMDeff_samples[,"NMD Consensus"] < NMDeff_median),"NMDeff_type"] <- "Low"
    NMDeff_samples[which(NMDeff_samples[,"NMD Consensus"] >= NMDeff_median),"NMDeff_type"] <- "High"
    NMDeff_samples$sample <- factor(NMDeff_samples$sample, levels = NMDeff_samples$sample)
    # Stack
    if (!is.null(cell_types_char)) {
        cols <- c(NMD_methods,cell_types,selected_genes)
    } else {
        cols <- c(NMD_methods,selected_genes)
    }
    #Correct by sample library size?
    NMDeff_samples[,selected_genes] <- NMDeff_samples[,selected_genes] / NMDeff_samples$sample_lib_size
    n <- length(cols)
    NMDeff_samples_stack <- stack(NMDeff_samples[,cols])
    NMDeff_samples_stack[NMDeff_samples_stack$ind %in% selected_genes,"values"] <- log2(NMDeff_samples_stack[NMDeff_samples_stack$ind %in% selected_genes,"values"])
    NMDeff_samples_stack$sample <- rep(NMDeff_samples$sample,n)
    NMDeff_samples_stack$NMDeff_type <- rep(NMDeff_samples$NMDeff_type,n)
    NMDeff_samples_stack[,tissue_var_strat] <- rep(NMDeff_samples[,tissue_var_strat],n)
    NMDeff_samples_stack[,tissue_var_strat] <- gsub(".*\\_","",NMDeff_samples_stack[,tissue_var_strat])
    NMDeff_samples_stack[,tissue_var] <- rep(NMDeff_samples[,tissue_var],n)
    NMDeff_samples_stack <- na.omit(NMDeff_samples_stack)
    tmp <- names(table(NMDeff_samples_stack$ind))
    NMD_method_char <- tmp[!tmp %in% selected_genes]
    NMDeff_samples_stack$ind <- factor(NMDeff_samples_stack$ind, levels = c(NMD_method_char,selected_genes))
    colnames(NMDeff_samples_stack) <- c("NMD_efficiency","NMD_method","sample","NMDeff_type",tissue_var_strat,tissue_var)
    # Order by median
    if (dataset == "GTEx") {
        NMDeff_medians <- aggregate(NMD_efficiency ~ eval(parse(text=tissue_var)) + NMD_method, data = NMDeff_samples_stack, median)
        NMDeff_medians <- NMDeff_medians[grep("NMD Consensus", NMDeff_medians$NMD_method),]
        NMDeff_medians <- NMDeff_medians[order(NMDeff_medians$NMD_efficiency),]
        NMDeff_samples_stack[,tissue_var] <- factor(NMDeff_samples_stack[,tissue_var], levels = unique(NMDeff_medians[,1]))
        if (!is.null(cell_types_char)) {
            panels_to_show <- c(NMD_methods,cell_types)
            char <- "cell_type_deconvolution"
        } else {
            panels_to_show <- c(NMD_methods,"UPF2","UPF3A","UPF3B","SMG1","SMG5","SMG7","SMG6","SMG8","MIR128-2","MIR128-1")
            char <- "NMDfactors_gene_expression"
        }
    } else if (dataset == "TCGA") {
        NMDeff_medians <- aggregate(NMD_efficiency ~ eval(parse(text=tissue_var_strat)) + NMD_method, data = NMDeff_samples_stack, median)
        NMDeff_medians <- NMDeff_medians[grep("NMD Consensus", NMDeff_medians$NMD_method),]
        NMDeff_medians <- NMDeff_medians[order(NMDeff_medians$NMD_efficiency),]
        NMDeff_samples_stack[,tissue_var_strat] <- factor(NMDeff_samples_stack[,tissue_var_strat], levels = unique(NMDeff_medians[,1]))
        panels_to_show <- c(NMD_methods,"UPF2","UPF3A","UPF3B","SMG1")
        char <- "NMDfactors_gene_expression"
    }
    
    # Plot
    # NMDeff_samples_stack_filt <- NMDeff_samples_stack %>% filter(NMD_method %in% c("UPF1","UPF2","UPF3A","UPF3B","SMG1","SMG5","SMG7","SMG6","SMG8"))
    # NMDeff_samples_stack_filt <- NMDeff_samples_stack %>% filter(NMD_method %in% c("SMG9","RBM8A","EIF3A","CASC3","EIF4A3","MAGOH"))
    NMDeff_samples_stack_filt <- NMDeff_samples_stack %>% filter(NMD_method %in% panels_to_show)

    p <- ggplot(data = NMDeff_samples_stack_filt, aes(x = eval(parse(text=tissue_var)), y = NMD_efficiency, fill = factor(eval(parse(text=tissue_var_strat))))) +
            geom_boxplot() + coord_flip() +
            facet_grid(. ~ NMD_method , scale = "free") +
            #geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
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
    plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/",dataset,"_NMDeff_and_",char,".png")  
    png(plot_path, width = 6500, height = 2500, res = 300)
    print(p)
    dev.off()

    return(NMDeff_samples_stack)    

}

selected_genes <- c("UPF1","UPF2","UPF3A","UPF3B","SMG1","SMG5","SMG7","SMG6","SMG8","SMG9","RBM8A","EIF3A","CASC3","EIF4A3","MAGOH","MIR128-2","MIR128-1")
cell_types <- c("astrocyte_of_the_cerebral_cortex","Bergmann_glial_cell","brain_pericyte","endothelial_cell","neuron","oligodendrocyte","oligodendrocyte_precursor_cell")
NMD_methods <- c("NMD Consensus", "PTC NMD-triggering 0.2")

# NMD factors gene expression
NMDeff_vs_cell_type_or_NMD_factors_exp(dataset = "TCGA", cancer_type_char = "all", NMD_methods = NMD_methods, 
                                        selected_genes = selected_genes, cell_types = cell_types, cell_types_char = NULL)
NMDeff_vs_cell_type_or_NMD_factors_exp(dataset = "GTEx", cancer_type_char = "all", NMD_methods = NMD_methods, 
                                        selected_genes = selected_genes, cell_types = cell_types, cell_types_char = NULL)
# Cell type deconvolution only in GTEx
NMDeff_vs_cell_type_or_NMD_factors_exp(dataset = "GTEx", cancer_type_char = "all", NMD_methods = NMD_methods, 
                                        selected_genes = selected_genes, cell_types = cell_types, cell_types_char = "yes")


NMDeff_vs_neurons <- function(dataset, var_to_plot) {

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

    if (var_to_plot != "neuron") {
        p <- ggplot(data = sample_NMD_efficiencies, aes(x = eval(parse(text=tissue_var_strat)), 
                            y = eval(parse(text=var_to_plot_char)), fill = factor(eval(parse(text=tissue_var))))) +
                geom_boxplot() + coord_flip(ylim = c(-3,3)) + #coord_cartesian(ylim = c(-3,3)) +
                #facet_grid(. ~ NMD_method , scale = "free") +
                geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
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
        plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/",dataset,"_iNMDeff_",var_to_plot_char,"_brains.png")  
        png(plot_path, width = 2500, height = 2500, res = 300)
        print(p)
        dev.off()
    } else if (var_to_plot == "neuron") {
        p <- ggplot(data = sample_NMD_efficiencies, aes(x = eval(parse(text=tissue_var_strat)), 
                            y = eval(parse(text=var_to_plot)), fill = factor(eval(parse(text=tissue_var))))) +
                geom_boxplot() + coord_flip(ylim = c(0,1)) + #coord_cartesian(ylim = c(-3,3)) +
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
        plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Brain/",dataset,"_",var_to_plot,"_brains.png")  
        png(plot_path, width = 2500, height = 2500, res = 300)
        print(p)
        dev.off()
    }
    # 2) Add NMDeff vs Neuron correlations
    cols <- c("ETG","ASE","neuron")
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
            corr_ETG_vs_neuron = cor_matrix["neuron", "ETG"],
            corr_ASE_vs_neuron = cor_matrix["neuron", "ASE"]
        )
    })
    # Bind all data frames together
    corr_final_df <- do.call(rbind, cor_dfs)
    # Merge
    sample_NMD_efficiencies <- merge(sample_NMD_efficiencies,corr_final_df, by.x = tissue_var_strat, by.y = "tissue_type", all.x = TRUE)
    return(sample_NMD_efficiencies)
}
# TCGA
tmp <- NMDeff_vs_neurons(dataset = "TCGA", var_to_plot = "NMD Consensus") 
tmp <- NMDeff_vs_neurons(dataset = "TCGA", var_to_plot = "PTC NMD-triggering 0.2") 
TCGA_sample_NMD_eff <- NMDeff_vs_neurons(dataset = "TCGA", var_to_plot = "neuron") 
# Remove low sample sizes
TCGA_sample_NMD_eff <- TCGA_sample_NMD_eff[!TCGA_sample_NMD_eff$cancer_type_strat %in%
                            c("LGG_NA","GBM_G-CIMP-high","GBM_G-CIMP-low","LGG_G-CIMP-low"),]
# Save
write.table(TCGA_sample_NMD_eff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(TCGA_sample_NMD_eff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2E.RData")

# GTEx
GTEx_sample_NMD_eff <- NMDeff_vs_neurons(dataset = "GTEx", var_to_plot = "neuron") 
# Remove low sample sizes
# GTEx_sample_NMD_eff <- GTEx_sample_NMD_eff[!GTEx_sample_NMD_eff$cancer_type_strat %in%
#                             c("LGG_NA","GBM_G-CIMP-high","GBM_G-CIMP-low","LGG_G-CIMP-low"),]
# Save
write.table(GTEx_sample_NMD_eff, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2F.txt", 
         vscode-terminal:/4a2942bc79609de7273ae7a93f6cb14e/10       sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_sample_NMD_eff, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2F.RData")
