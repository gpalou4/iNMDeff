# conda activate /home/dnaro/.conda/envs/NMD_regression_3
# library("edgeR")
.libPaths( rev( .libPaths() ) )
# Bayesian for NB regression
library("V8")
# library("shinystan")
library("rstan")
library("rstanarm")
library("readxl")
library("dplyr")
library("ggpubr")
library("ggsignif")

################################ SCRIPT ################################
################################ SCRIPT ################################
################################ SCRIPT ################################

# 1.1) NMD targets Consensus
path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_old.txt"
NMD_targets <- read.table(file = path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 1.2) cNMDeff estimates excluding 1 gene at a time

NMD_genes <- c("all",unique(NMD_targets$ensembl_gene_id))
NMD_genes <- unique(NMD_targets$ensembl_gene_id)

# 1.3) 
CL_metadata_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/CL_metadata_subset_and_controls.xlsx"
CL_metadata <- data.frame(read_excel(CL_metadata_path))
# Same order as CL data
CL_metadata$ID <- factor(paste0(CL_metadata$GSE_project,"_R",CL_metadata$replicate))
CL_metadata$Type <- as.factor(CL_metadata$Type)
CL_metadata <- CL_metadata[order(CL_metadata$ID),]

cNMDeff_final_df <- c()
for (NMD_gene in NMD_genes) {

    path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/excluding_genes/cNMDeff_all.txt")
    cNMDeff_all <- read.table(file = path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/excluding_genes/cNMDeff_",NMD_gene,".txt")
    cNMDeff_gene <- read.table(file = path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Difference between original and the new
    # df_all <- rbind(df_all,data.frame(delta = cNMDeff_all$nb_coeff - cNMDeff_gene$nb_coeff, 
    #         excluded_NMD_gene = NMD_gene,
    #         cell_line = cNMDeff_gene$cell_line))
    cNMDeff_final_df <- rbind(cNMDeff_final_df,cNMDeff_gene)

}
cNMDeff_final_df <- rbind(cNMDeff_final_df,cNMDeff_all)

df <- merge(cNMDeff_final_df,CL_metadata, 
            by.x = "cell_line", by.y = "ID", all.x = TRUE)
df$nb_coeff <- -df$nb_coeff
df$CI_2.5 <- -df$CI_2.5
df$CI_97.5 <- -df$CI_97.5
df$Type <- ifelse(df$Type == "control","WT","UPF1 KD")

path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/ETG_cNMDeff_cell_types.txt")
write.table(df, file = path,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# Fig 1C
path <- paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig1/Fig1C.txt")
write.table(df, file = path,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# Supp. Fig. 3
path <- paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3_new/SuppFig3A.txt")
write.table(df, file = path,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# 4) Analysis

df <- merge(cNMDeff_final_df,CL_metadata[,c("ID","Type","CL")], 
            by.x = "cell_line", by.y = "ID", all.x = TRUE)

df$nb_coeff <- -df$nb_coeff
df$CI_2.5 <- -df$CI_2.5
df$CI_97.5 <- -df$CI_97.5
df$Type <- ifelse(df$Type == "control","WT","UPF1 KD")

# Calculate mean and CI for each group (Type and CL)
summary_df <- df %>%
  group_by(Type, CL) %>%
  summarize(
    nb_coeff = nb_coeff,
    mean_nb_coeff = mean(nb_coeff),
    sd_nb_coeff = sd(nb_coeff),
    n = n(),
    se_nb_coeff = sd_nb_coeff / sqrt(n),   # Standard error
    # CI_lower = mean_nb_coeff - qt(0.975, df = n - 1) * se_nb_coeff, # 95% CI lower
    # CI_upper = mean_nb_coeff + qt(0.975, df = n - 1) * se_nb_coeff  # 95% CI upper
    CI_lower = mean(CI_2.5),
    CI_upper = mean(CI_97.5)
  )

plot <- ggplot(df, aes(x = CL, y = nb_coeff, fill = Type)) +
    geom_bar(
    data = df %>%
        group_by(CL, Type) %>%
        summarise(mean_nb_coeff = mean(nb_coeff), .groups = "drop"),
    aes(y = mean_nb_coeff),
        stat = "identity", position = position_dodge(width = 0.9), width = 0.7
    ) +
    geom_errorbar(data = summary_df,
                aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.9), width = 0.2) +
    geom_point(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
        size = 3, shape = 21, color = "black"
    ) +
    # Add comparison bars
    annotate("rect", xmin = 0.78, xmax = 1.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 0.78, xmax = 0.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 1.22, xmax = 1.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    annotate("rect", xmin = 1.78, xmax = 2.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 1.78, xmax = 1.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 2.22, xmax = 2.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    annotate("rect", xmin = 2.78, xmax = 3.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 2.78, xmax = 2.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 3.22, xmax = 3.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    coord_cartesian(ylim = c(0,4)) +
    stat_compare_means(
        method = "t.test",
        label = "p.format",
        paired = TRUE,
        label.y = max(df$nb_coeff) + 0.5  # Adjust for p-value placement
    ) +
    labs(
        title = "",
        x = "Cell Lines",
        y = "ETG cNMDeff"
    ) +
    theme_classic() +
    scale_fill_manual(values = c("WT" = "#56B4E9", "UPF1 KD" = "#E69F00")) +
    theme(
        legend.title = element_text(size = 10),
        legend.position = "top"
    )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/cNMDeff_UPF1_KD_with_excluded_genes.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

df <- merge(cNMDeff_final_df,CL_metadata[,c("ID","Type","CL")], 
            by.x = "cell_line", by.y = "ID", all.x = TRUE) %>%
            filter(NMD_gene_excluded == "all")

df$nb_coeff <- -df$nb_coeff
df$CI_2.5 <- -df$CI_2.5
df$CI_97.5 <- -df$CI_97.5
df$Type <- ifelse(df$Type == "control","WT","UPF1 KD")

# aggregate(nb_coeff ~ Type, data = df , median)

# Calculate mean and CI for each group (Type and CL)
summary_df <- df %>%
  group_by(Type, CL) %>%
  summarize(
    nb_coeff = nb_coeff,
    mean_nb_coeff = mean(nb_coeff),
    sd_nb_coeff = sd(nb_coeff),
    n = n(),
    se_nb_coeff = sd_nb_coeff / sqrt(n),   # Standard error
    # CI_lower = mean_nb_coeff - qt(0.975, df = n - 1) * se_nb_coeff, # 95% CI lower
    # CI_upper = mean_nb_coeff + qt(0.975, df = n - 1) * se_nb_coeff  # 95% CI upper
    CI_lower = mean(CI_2.5),
    CI_upper = mean(CI_97.5)
  )

plot <- ggplot(df, aes(x = CL, y = nb_coeff, fill = Type)) +
    geom_bar(
    data = df %>%
        group_by(CL, Type) %>%
        summarise(mean_nb_coeff = mean(nb_coeff), .groups = "drop"),
    aes(y = mean_nb_coeff),
        stat = "identity", position = position_dodge(width = 0.9), width = 0.7
    ) +
    geom_errorbar(data = summary_df,
                aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.9), width = 0.2) +
    geom_point(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
        size = 3, shape = 21, color = "black"
    ) +
    # Add comparison bars
    annotate("rect", xmin = 0.78, xmax = 1.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 0.78, xmax = 0.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 1.22, xmax = 1.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    annotate("rect", xmin = 1.78, xmax = 2.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 1.78, xmax = 1.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 2.22, xmax = 2.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    annotate("rect", xmin = 2.78, xmax = 3.22, ymin = 3.15, ymax = 3.15, alpha=1,colour = "black")+
    annotate("rect", xmin = 2.78, xmax = 2.78, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 3.22, xmax = 3.22, ymin = 3.05, ymax = 3.15, alpha=1, colour = "black")+

    coord_cartesian(ylim = c(0,4)) +
    stat_compare_means(
        method = "t.test",
        label = "p.format",
        paired = TRUE,
        label.y = max(df$nb_coeff) + 0.5  # Adjust for p-value placement
    ) +
    labs(
        title = "",
        x = "Cell Lines",
        y = "ETG cNMDeff"
    ) +
    theme_classic() +
    scale_fill_manual(values = c("WT" = "#56B4E9", "UPF1 KD" = "#E69F00")) +
    theme(
        legend.title = element_text(size = 10),
        legend.position = "top"
    )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/cNMDeff_UPF1_KD.png"
ggsave(final_figure_path, plot, width = 175, height = 110, units = "mm") 

 




