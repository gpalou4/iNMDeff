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

library(dplyr)
library(RColorBrewer)
library(ggplot2)

# 1) Data

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

# 2) Check gene expression of well-known NMD-targets (controls) for every NMDeff bin

NMD_targets_exp_vs_NMDeff <- function(dataset, cancer_type_char, NMD_method, selected_genes) {

    # Dataset
    if (dataset == "TCGA") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_TCGA 
        RNAseq_TPM_all <- RNAseq_TCGA_TPM_all
    } else if (dataset == "GTEx") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_GTEx
        RNAseq_TPM_all <- RNAseq_GTEx_TPM_all
        sample_NMD_efficiencies$sample <- NULL
        sample_NMD_efficiencies$sample <- sample_NMD_efficiencies$sample_full_barcode
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
    NMDeff_samples <- NMDeff_samples[order(NMDeff_samples[,NMD_method]),]
    NMDeff_median <- median(NMDeff_samples[,NMD_method],na.rm = TRUE)
    NMDeff_samples[which(NMDeff_samples[,NMD_method] < NMDeff_median),"NMDeff_type"] <- "Low"
    NMDeff_samples[which(NMDeff_samples[,NMD_method] >= NMDeff_median),"NMDeff_type"] <- "High"
    NMDeff_samples$sample <- factor(NMDeff_samples$sample, levels = NMDeff_samples$sample)

    # for (selected_gene in selected_genes) {
    #     # Stack
    #     NMDeff_samples_stack <- stack(NMDeff_samples[,c(NMD_method,selected_gene)])
    #     NMDeff_samples_stack[NMDeff_samples_stack$ind == selected_gene,"values"] <- sqrt(NMDeff_samples_stack[NMDeff_samples_stack$ind == selected_gene,"values"])
    #     NMDeff_samples_stack$ind <- ifelse(as.character(NMDeff_samples_stack$ind) == "endogenous_NMD_Consensus","NMD Consensus",as.character(NMDeff_samples_stack$ind))
    #     NMDeff_samples_stack$ind <- ifelse(as.character(NMDeff_samples_stack$ind) == "ASE_PTC_NMD_triggering_0.2","PTC NMD-triggering",as.character(NMDeff_samples_stack$ind))
    #     NMDeff_samples_stack$sample <- rep(NMDeff_samples$sample,2)
    #     NMDeff_samples_stack$NMDeff_type <- rep(NMDeff_samples$NMDeff_type,2)
    #     NMDeff_samples_stack <- na.omit(NMDeff_samples_stack)
    #     tmp <- names(table(NMDeff_samples_stack$ind))
    #     NMD_method_char <- tmp[tmp != selected_gene]
    #     NMDeff_samples_stack$ind <- factor(NMDeff_samples_stack$ind, levels = c(NMD_method_char,selected_gene))

    #     #p <- ggplot(data = NMDeff_samples, aes(x = factor(sample, levels = sample), y = sqrt(eval(parse(text=selected_gene))), color = factor(NMDeff_type, levels = c("Low","High")))) +
    #     p <- ggplot(data = NMDeff_samples_stack, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
    #         geom_point(data = subset(NMDeff_samples_stack, ind != selected_gene),alpha=0.5, size = 2) + #ggtitle(selected_gene) +
    #         # geom_point(data = subset(NMDeff_samples_stack, ind == selected_gene),alpha=0.1, size = 1) + #ggtitle(selected_gene) +
    #         geom_smooth(data = subset(NMDeff_samples_stack, ind == selected_gene),
    #                     aes(group = factor(NMDeff_type, levels = c("Low","High"))), 
    #                     method = "gam", se = TRUE, linewidth = 2.5) +
    #         # geom_line(aes(group = 1)) +
    #         facet_grid(ind ~ ., scales = "free_y") +#, space = "free_y") +
    #         xlab("Individuals sorted by iNMDeff") + ylab("TPM       iNMDeff") +
    #         theme_classic(base_size = 20) + 
    #         theme(plot.title = element_text(hjust = 0.5, size = 35),
    #                 axis.title.x = element_text(color="black", size=35),
    #                 axis.title.y = element_text(color="black", size=35),
    #                 panel.spacing = grid::unit(2, "lines"),
    #                 #axis.text.x = element_text(color="black", size=15, angle = 90, hjust = 0.5, vjust = 0.5),
    #                 axis.text.x = element_blank(),
    #                 axis.text.y = element_text(color="black", size=30),
    #                 strip.text = element_text(size = 22),
    #                 legend.position = "top",
    #                 legend.text = element_text(color="black", size=30),
    #                 legend.title = element_text(color="black", size=35)) +
    #         #scale_color_brewer(palette = "Accent", labels = c("Low","High"), direction = -1) +
    #         guides(color = guide_legend(override.aes = list(size = 13), title = "iNMDeff")) +
    #         scale_color_manual(  
    #                 values = c("Low" = "#2D3263", "High" = "#F07626"), 
    #                 labels = c("Low", "High")
    #                 )
    #     plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/sample_level/NMD_targets_controls/",dataset,"_",cancer_type_char,"_",NMD_method,"_",selected_gene,".png")
    #     png(plot_path, width = 3500, height = 2500, res = 300)
    #     print(p)
    #     dev.off()
    # }
    # Stack
    n <- length(selected_genes) + 1
    NMDeff_samples_stack <- stack(NMDeff_samples[,c(NMD_method,selected_genes)])
    NMDeff_samples_stack[NMDeff_samples_stack$ind %in% selected_genes,"values"] <- sqrt(NMDeff_samples_stack[NMDeff_samples_stack$ind %in% selected_genes,"values"])
    NMDeff_samples_stack$ind <- ifelse(as.character(NMDeff_samples_stack$ind) == "endogenous_NMD_Consensus","NMD Consensus",as.character(NMDeff_samples_stack$ind))
    NMDeff_samples_stack$ind <- ifelse(as.character(NMDeff_samples_stack$ind) == "ASE_PTC_NMD_triggering_0.2","PTC NMD-triggering",as.character(NMDeff_samples_stack$ind))
    NMDeff_samples_stack$sample <- rep(NMDeff_samples$sample,n)
    NMDeff_samples_stack$NMDeff_type <- rep(NMDeff_samples$NMDeff_type,n)
    NMDeff_samples_stack <- na.omit(NMDeff_samples_stack)
    tmp <- names(table(NMDeff_samples_stack$ind))
    NMD_method_char <- tmp[!tmp %in% selected_genes]
    NMDeff_samples_stack$ind <- factor(NMDeff_samples_stack$ind, levels = c(NMD_method_char,selected_genes))

    #p <- ggplot(data = NMDeff_samples, aes(x = factor(sample, levels = sample), y = sqrt(eval(parse(text=selected_gene))), color = factor(NMDeff_type, levels = c("Low","High")))) +
    p <- ggplot(data = NMDeff_samples_stack, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
        geom_point(data = subset(NMDeff_samples_stack, !ind %in% selected_genes),alpha=0.5, size = 2) + #ggtitle(selected_gene) +
        # geom_point(data = subset(NMDeff_samples_stack, ind == selected_gene),alpha=0.1, size = 1) + #ggtitle(selected_gene) +
        geom_smooth(data = subset(NMDeff_samples_stack, ind %in% selected_genes),
                    aes(group = factor(NMDeff_type, levels = c("Low","High"))), 
                    method = "gam", se = TRUE, linewidth = 2.5) +
        # geom_line(aes(group = 1)) +
        facet_grid(ind ~ ., scales = "free_y") +#, space = "free_y") +
        xlab("Individuals sorted by iNMDeff") + ylab("TPM       iNMDeff") +
        theme_classic(base_size = 20) + 
        theme(plot.title = element_text(hjust = 0.5, size = 35),
                axis.title.x = element_text(color="black", size=35),
                axis.title.y = element_text(color="black", size=35),
                panel.spacing = grid::unit(2, "lines"),
                #axis.text.x = element_text(color="black", size=15, angle = 90, hjust = 0.5, vjust = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_text(color="black", size=30),
                strip.text = element_text(size = 22),
                legend.position = "top",
                legend.text = element_text(color="black", size=30),
                legend.title = element_text(color="black", size=35)) +
        #scale_color_brewer(palette = "Accent", labels = c("Low","High"), direction = -1) +
        guides(color = guide_legend(override.aes = list(size = 13), title = "iNMDeff")) +
        scale_color_manual(  
                values = c("Low" = "#2D3263", "High" = "#F07626"), 
                labels = c("Low", "High")
                )
                
    plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/sample_level/NMD_targets_controls/",dataset,"_",cancer_type_char,"_",NMD_method,"_",paste0(selected_genes,collapse="_"),".png")
    png(plot_path, width = 3500, height = 2500, res = 300)
    print(p)
    dev.off()

    return(NMDeff_samples_stack)    

}



neg_controls <- c("CASC11","TSPAN6","DPM1")


NMD_targets_UPR_response <- c("TRAF2","FSD1L","CNPY3","ERN1","ATF4","ATF3")
NMD_targets <- c("TRA2B","PTBP2","GAS5","RP9P")
NMD_genes <- c("UPF1","UPF2","SMG7","SMG5","SMG6","SMG8","SMG9","SMG1","UPF3B","UPF3A","RBM8A","EIF3A","CASC3","EIF4A3","MAGOH")
selected_genes <- c(NMD_targets,NMD_genes,NMD_targets_UPR_response)

selected_genes <- c("RP9P","GAS5","SMG5")

# TCGA STAD
# NMD_targets_exp_vs_NMDeff(dataset = "TCGA", NMD_method = "endogenous_NMD_Consensus", selected_genes = selected_genes, cancer_type_char = "STAD")

# All TCGA
NMDeff_gene_exp_stacked <- NMD_targets_exp_vs_NMDeff(dataset = "TCGA", NMD_method = "endogenous_NMD_Consensus", selected_genes = selected_genes, cancer_type_char = "all")

write.table(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1D.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
write.table(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2D.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# save R object
saveRDS(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1D.RData")
saveRDS(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2D.RData")
NMD_targets_exp_vs_NMDeff(dataset = "TCGA", NMD_method = "ASE_PTC_NMD_triggering_0.2", selected_genes = selected_genes, cancer_type_char = "all")
# write.table(NMDeff_ASE_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Supp.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# All GTEx
NMDeff_gene_exp_stacked <- NMD_targets_exp_vs_NMDeff(dataset = "GTEx", NMD_method = "endogenous_NMD_Consensus", selected_genes = selected_genes, cancer_type_char = "all")
write.table(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1E.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
write.table(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2E.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# save R object
saveRDS(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1E.RData")
saveRDS(NMDeff_gene_exp_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/SuppFig2E.RData")
NMD_targets_exp_vs_NMDeff(dataset = "GTEx", NMD_method = "ASE_PTC_NMD_triggering_0.2", selected_genes = selected_genes, cancer_type_char = "all")
# write.table(NMDeff_ASE_stacked, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Supp.txt", 
#             sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)


######################################################################################################## 

# Mirar transcript level
# Mirar proteomic level?

# NMD_Consensus <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final.txt",
#             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# NMD_target_genes <- unique(NMD_Consensus$gene_symbol)


# df <- round(cor(NMDeff_samples[,selected_genes]),2)
# df <- df[,"SMG1", drop = FALSE]
# df <- df[order(df[,"SMG1"], decreasing = TRUE),,drop = FALSE]
# summary(df)


# lm_res <- lm(endogenous_NMD_Consensus ~ SMG1 + SMG5 + SMG6 + SMG7 + SMG8 + UPF1 + UPF3B + UPF3A + tissue + factor(sex) + age + factor(death_group)
#                 + sample_lib_size + sample, data = NMDeff_samples)   
# summary(lm_res)

########################################################################################################