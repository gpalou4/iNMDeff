source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

# Function to save ggsurv plots with "ggsave"
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)} 

ggplot_theme_bw <- function() {

  theme_bw() +
  theme(

    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.title.align = 0.5,
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,0,0,0),
    legend.text = element_text(colour = "black", size = 8),

    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 12, hjust = 0.5),

    plot.title = element_text(size = 14, hjust = 0.5),
     strip.background = element_rect(fill = "#f4f1f8e4", colour = "#000000"),

    strip.text = element_text(colour = "black", size = 12)
  )
}

################################
############ FIG 6A ############
################################

# Data
input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/panel_A.RData")
# input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6A_subset.RData")
# Change
input_figureA$mutation_type <- as.character(input_figureA$mutation_type)
input_figureA$mutation_type <- ifelse(input_figureA$mutation_type == "wmis_cv", "Missense",input_figureA$mutation_type)
input_figureA$mutation_type <- ifelse(input_figureA$mutation_type == "wnon_cv", "Nonsense",input_figureA$mutation_type)
input_figureA$mutation_type <- ifelse(input_figureA$mutation_type == "wind_cv", "Indel",input_figureA$mutation_type)
table(input_figureA$mutation_type)
# Median differences manually

# dndscv_wilcox_test <- input_figureA %>%
#   group_by(mutation_type,gene_type,NMD_variant_type,percentile,NMD_method) %>%
#   filter(length(unique(NMDeff)) == 2) %>%  # Ensure exactly two levels
#   summarise(wilcox_p_value = wilcox.test(dNdScv_ratio ~ NMDeff,na.rm =TRUE)$p.value)
dndscv_medians <- input_figureA %>%
  group_by(mutation_type,NMDeff,gene_type,NMD_variant_type,percentile,NMD_method) %>%
  summarise(dNdScv_median = median(dNdScv_ratio, na.rm = TRUE))
# dndscv_wilcox_test %>%
#   filter(mutation_type == "Missense" & gene_type == "TSG" & NMD_method == "END" &variants_NMDeff_TCGA_df$ percentile == 50) %>%
#   as.data.frame()
dndscv_medians %>%
  filter(mutation_type == "Nonsense" & gene_type == "oncogene" & NMD_method == "END" & percentile == 50) %>%
  as.data.frame()

list_plots <- list() 
combinations <- combn(names(table(input_figureA$NMDeff)), 2, simplify = FALSE)
for (NMD_method_char in c("END","ASE")) {
  for (mutation_type_char in c("Nonsense","Missense")) {
    if ( mutation_type_char == "Nonsense" ) {
      if (NMD_method_char == "ASE") { 
            ylim <- c(0,6)
      } else if (NMD_method_char == "END") {
            ylim <- c(0,6)
      }
    } else if ( mutation_type_char == "Missense" ) {
      if (NMD_method_char == "ASE") { 
        ylim <- c(0,3.5)
      } else if (NMD_method_char == "END") {
        ylim <- c(0,3)
      }
    }
    p <- input_figureA %>%
        mutate(NMD_variant_type = factor(NMD_variant_type, levels = c("NMD-triggering","NMD-evading"))) %>%
        filter(gene_type == "TSG") %>%
        # filter(percentile %in% c(20,30,40,50)) %>%
        filter(percentile %in% c(50)) %>%
        # filter(!gene_name %in% low_higher_genes$gene_name) %>%
        # filter(n_non >= 3 & n_syn >= 1) %>%
        # filter(qtrunc_cv < 0.50) %>%
        filter(NMD_method %in% NMD_method_char) %>%
        filter(mutation_type %in% mutation_type_char) %>%
          ggplot(aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                  geom_violin() + scale_fill_brewer(palette = "Pastel1") + xlab("") + ylab("dNdS ratio") +
                  geom_jitter(aes(fill = factor(NMDeff)), 
                    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                    alpha = 0.1, size = 2) + xlab("") +
                  labs(fill = "ETG iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                  ggtitle(paste0(mutation_type_char," in TSG")) +
                  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                  facet_grid(. ~ NMD_variant_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                  geom_boxplot(width=0.3, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                  ggplot_theme_bw() +
                  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
                  # theme(legend.key.size = unit(0.1, "lines")) +
                  theme(strip.text = element_text(size = 9))+
                stat_compare_means(
                                  # data=paired_data,
                                  # aes(x=dNdScv_ratio_high,y=dNdScv_ratio_low),
                                  aes(group = NMDeff),
                                  method.args = list(alternative = "greater"),
                                  comparisons = combinations, 
                                  size = 5, paired = FALSE,
                                  label.y = c(1.25),
                                  label = "p.format", method = "wilcox.test", hide.ns = FALSE)
      list_plots[[length(list_plots) + 1]] <- p
  }
}

# input_figureA_filt <- input_figureA %>% 
#           filter(mutation_type == "Nonsense" & gene_type == "TSG" &
#                 percentile == 20 & NMD_method == "END" & 
#                 NMD_variant_type == "NMD-triggering") %>%
#           filter(n_non >= 5)

# table(input_figureA_filt$ptrunc_cv < 0.25,input_figureA_filt$NMDeff)
# 
wilcox_test <- function(input_figureA, mutation_type_char, gene_type_char, NMD_variant_type_char,
                        NMD_method_char, percentile_char) {
  input_figureA_filt <- input_figureA %>% 
            filter(mutation_type == mutation_type_char & gene_type == gene_type_char &
                  percentile == percentile_char & NMD_method == NMD_method_char & 
                  NMD_variant_type == NMD_variant_type_char)
  # input_figureA_filt <- input_figureA_filt %>%
  #               filter(n_non >= 10)
  # Subset data for 'High' and 'Low' NMDeff
  high_group <- subset(input_figureA_filt, NMDeff == "High")
  colnames(high_group)[colnames(high_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_high"
  low_group <- subset(input_figureA_filt, NMDeff == "Low")
  colnames(low_group)[colnames(low_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_low"
  # Merge subsets by 'gene_name'
  paired_data <- merge(high_group[,c("gene_name","dNdScv_ratio_high")], low_group[,c("gene_name","dNdScv_ratio_low")], by = c("gene_name"))       
  wilcox_res <- wilcox.test(paired_data$dNdScv_ratio_high, paired_data$dNdScv_ratio_low, paired = TRUE, 
                            method = "wilcox.test", alternative = "greater")
  # t.test(log2(dNdScv_ratio) ~ NMDeff, data = input_figureA_filt, paired = FALSE, alternative = "greater")
  return(wilcox_res$p.value)
}

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, list_plots[[1]], width = 250, height = 200, units = "mm") 

# Wilcoxon tests with pairing data
pvalue1 <- wilcox_test(input_figureA, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
# pvalue2 <- pvalue1
pvalue2 <- wilcox_test(input_figureA, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
legend_A <- get_legend(list_plots[[1]])
plot_A_1 <- ggplot_build(list_plots[[1]] + guides(fill = FALSE))
plot_A_1$data[[5]][,"annotation"] <- as.character(plot_A_1$data[[5]][,"annotation"])
plot_A_1$data[[5]][,"annotation"][1:3] <- rep(round(pvalue1,2),3)
plot_A_1$data[[5]][,"annotation"][4:6] <- rep(round(pvalue2,2),3)
plot_A_1$data[[5]][,"annotation"] <- as.factor(plot_A_1$data[[5]][,"annotation"])
plot_A_1 <- ggplot_gtable(plot_A_1)
plot_A_1 <- ggplotify::as.ggplot(plot_A_1)
pvalue3 <- wilcox_test(input_figureA, mutation_type_char = "Indel", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
pvalue4 <- wilcox_test(input_figureA, mutation_type_char = "Missense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
plot_A_2 <- ggplot_build(list_plots[[2]] + guides(fill = FALSE))
plot_A_2$data[[5]][,"annotation"] <- as.character(plot_A_2$data[[5]][,"annotation"])
plot_A_2$data[[5]][,"annotation"][1:3] <- rep(round(pvalue3,2),3)
plot_A_2$data[[5]][,"annotation"][4:6] <- rep(round(pvalue4,2),3)
plot_A_2$data[[5]][,"annotation"] <- as.factor(plot_A_2$data[[5]][,"annotation"])
plot_A_2 <- ggplot_gtable(plot_A_2)
plot_A_2 <- ggplotify::as.ggplot(plot_A_2)

# plot_A <- ggarrange(plot_A_1, legend.grob = legend_A,
plot_A <- ggarrange(plot_A_1,plot_A_2, labels = c("",""), legend.grob = legend_A,
                font.label = list(size = 20), #widths = c(0.45,0.28,0.28),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
# plot_A <- annotate_figure(plot_A, top = text_grob("Tumor Supressor Genes", 
#                             size = 16, face = "plain")) #+ ggplot_theme_bw()
       
# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, plot_A, width = 250, height = 200, units = "mm") 

################################
############ FIG 6B ############
################################

# Data
input_figureB <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/panel_B.RData")

# Calculate Wilcoxon test for each facet
facet_pvalues <- input_figureB %>%
  group_by(cancer_type_new) %>% # Group by facets
  rstatix::wilcox_test(
    quanTIseq_lsei_TIL10 ~ NMD_type,
    alternative = "less" # Always test High < Low
  ) %>%
  mutate(p.format = format(p, digits = 2, nsmall = 2)) %>% # Add formatted p-value column
  mutate(NMD_type = "High")

# Plot quanTIseq_lsei_TIL10 with percentiles on the x-axis
plot_B <- ggplot(input_figureB, aes(x = as.factor(NMD_type), y = quanTIseq_lsei_TIL10 * 100, fill = NMD_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
  geom_point(
    aes(color = NMD_type), 
    position = position_dodge(0.8), 
    size = 2, 
    alpha = 0.7
  ) +
  facet_wrap(.~cancer_type_new, scale = "free_y") +
  labs(
    title = "",
    color = "",
    fill = "ETG iNMDeff",
    x = "",
    y = "CD8 T %"
  ) +
  ggplot_theme() +
  guides(color = "none") + # Remove the color legend
  coord_cartesian(ylim = c(0, 10)) +
  theme(axis.text.x = element_blank(),
    legend.text = element_text(colour = "black", size = 10)) +
  stat_pvalue_manual(
    facet_pvalues, 
    label = "p.format", # Use formatted p-values
    y.position = 7
  )

##################################
############ FIG 6C ############
##################################

# Data
input_figureC_E <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/panel_C_E.RData")

surv_obj <- input_figureC_E[["SKCM_no_treatment_PFS"]]$surv_obj
km_fit <- input_figureC_E[["SKCM_no_treatment_PFS"]]$km_fit
df <- input_figureC_E[["SKCM_no_treatment_PFS"]]$df
# Correlation
cor.test(df$endogenous_NMD_Consensus,df$PFS_months, method = "pearson")
cor.test(df$ASE_PTC_NMD_triggering_0.2,df$PFS_months, method = "pearson")

plot_C <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS",
        title = "TCGA - SKCM \nw/o available treament",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_C_to_save <- ggsave_workaround(plot_C)

################################
############ FIG 6D ############
################################

surv_obj <- input_figureC_E[["SKCM_immunotherapy_PFS"]]$surv_obj
km_fit <- input_figureC_E[["SKCM_immunotherapy_PFS"]]$km_fit
df <- input_figureC_E[["SKCM_immunotherapy_PFS"]]$df
# Correlation
cor.test(df$endogenous_NMD_Consensus,df$PFS_months, method = "pearson")
cor.test(df$ASE_PTC_NMD_triggering_0.2,df$PFS_months, method = "pearson")

plot_D <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS",
        title = "TCGA - SKCM with\n at least immunotherapy",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_D_to_save <- ggsave_workaround(plot_D)

################################
############ FIG 6E ############
################################

surv_obj <- input_figureC_E[["SKCM_chemotherapy_PFS"]]$surv_obj
km_fit <- input_figureC_E[["SKCM_chemotherapy_PFS"]]$km_fit
df <- input_figureC_E[["SKCM_chemotherapy_PFS"]]$df
# Correlation
cor.test(df$endogenous_NMD_Consensus,df$PFS_months, method = "pearson")
cor.test(df$ASE_PTC_NMD_triggering_0.2,df$PFS_months, method = "pearson")

plot_E <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS",
        title = "TCGA - SKCM with\n chemotherapy, w/o immunotherapy",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_E_to_save <- ggsave_workaround(plot_E)

################################
############ FIG 6F ############
################################

# Data
input_figureF_H <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/panel_F_H.RData")

surv_obj <- input_figureF_H[["Liu_SKCM_PFS"]][["endogenous_NMD_Consensus"]]$surv_obj
km_fit <- input_figureF_H[["Liu_SKCM_PFS"]][["endogenous_NMD_Consensus"]]$km_fit
df <- input_figureF_H[["Liu_SKCM_PFS"]][["endogenous_NMD_Consensus"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS, method = "pearson")

plot_F <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "PFS",
        title = "Liu 2019 - SKCM\nwith immunotherapy (P20)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_F_to_save <- ggsave_workaround(plot_F)

################################
############ FIG 6G ############
################################

surv_obj <- input_figureF_H[["Carrol_EAC_PFS"]][["endogenous_NMD_Consensus"]]$surv_obj
km_fit <- input_figureF_H[["Carrol_EAC_PFS"]][["endogenous_NMD_Consensus"]]$km_fit
df <- input_figureF_H[["Carrol_EAC_PFS"]][["endogenous_NMD_Consensus"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS, method = "pearson")

plot_G <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (days)", ylab = "PFS",
        title = "Carrol 2023 - EAC\nwith immunotherapy (P25)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_G_to_save <- ggsave_workaround(plot_G)

################################
############ FIG 6H ############
################################

surv_obj <- input_figureF_H[["Motzer_RCC_PFS"]][["endogenous_NMD_Consensus"]]$surv_obj
km_fit <- input_figureF_H[["Motzer_RCC_PFS"]][["endogenous_NMD_Consensus"]]$km_fit
df <- input_figureF_H[["Motzer_RCC_PFS"]][["endogenous_NMD_Consensus"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS_P, method = "pearson")

plot_H <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS",
        title = "Motzer 2020 - RCC\nwith immunotherapy (P30)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_H_to_save <- ggsave_workaround(plot_H)

################################
############ FIG 6I ############
################################

# Data
input_figureI <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/panel_I.RData")

plot_I <- input_figureI %>%
            filter(!dataset %in% "Hartwig_all") %>%
            ggplot(aes(x = factor(dataset_n), y = AUC, fill = factor(model))) + 
                geom_bar( stat = 'identity', position = position_dodge(width=0.9) ) +
                # Add error bars for confidence intervals
                geom_errorbar(aes(ymin = AUC_CI_low, ymax = AUC_CI_high), 
                               position = position_dodge(width=0.9), width = 0.2) +  # Adjust width as needed
                facet_grid(NMD_method ~ .) + ylim(c(0,1))+
                geom_hline(yintercept=0.5, linetype="dashed", color = "black", size=1) +
                labs(title = "", x = "", y = "AUC", fill = "Model") +                 
                ggplot_theme_bw() + scale_fill_brewer(palette = "OrRd", direction = -1) +
                theme(legend.position = "top", 
                        # strip.background = element_rect(),
                        legend.text = element_text(colour = "black", size = 10),
                      axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.75, size = 12),
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
                geom_text(aes(label = ifelse(AUC_p_value < 0.05, "*", "")), 
                                        position = position_dodge(width=0.9), size = 13, color = "black", hjust = 0.5, vjust = 0.5)

# For the Manuscript

input_figureI %>% 
        filter(!dataset %in% "Hartwig_all") %>%
        group_by(NMD_method, model) %>% 
        summarize(
                min(AUC, na.rm = TRUE),
                max(AUC, na.rm = TRUE),
                mean(AUC, na.rm = TRUE),
                mean(accuracy, na.rm = TRUE),
                mean(specificity, na.rm = TRUE),
                mean(sensitivity, na.rm = TRUE),
                mean(precision, na.rm = TRUE)) %>% as.data.frame()
# Increase in AUC for each dataset separately
input_figureI %>% 
        filter(!dataset %in% "Hartwig_all") %>%
        group_by(NMD_method,dataset) %>%
        summarize(
                AUC_diff = AUC[model == "TMB,iNMDeff,TMB*iNMDeff"] -  AUC[model == "TMB"],
        ) %>% arrange(dataset)

###############################################
############ PANEL UP --> FIG 6A ##############
###############################################

# Create blank filler plots
blank_space <- ggplot() + theme_void()

# Center plot_A in the middle of the row
# fig_up <- plot_grid(
#   blank_space, 
#   plot_A, 
#   blank_space, 
#   ncol = 3, 
#   rel_widths = c(0.2, 0.55, 0.2), # Adjust these values to control spacing
#   labels = c("", "A", ""), 
#   label_size = 20
# )

fig_up <- plot_grid(
  plot_A, 
  plot_B, 
  ncol = 2, 
  rel_widths = c(0.5,0.5), # Adjust these values to control spacing
  labels = c("A", "B"), 
  label_size = 20
)

###############################################
############ PANEL MID --> FIG C-H ############
###############################################

fig_mid <- plot_grid(plot_C_to_save, plot_D_to_save, plot_E_to_save,
                      plot_F_to_save, plot_G_to_save, plot_H_to_save, 
                      nrow = 2, ncol = 3,  
                      labels = c("C","D","E","F","G","H"), 
                      label_size = 20)

##################################################
############ PANEL BOTTOM --> FIG H ##############
##################################################

# Center plot_I in the middle of the row
fig_bottom <- plot_grid(
  blank_space, 
  plot_I, 
  blank_space, 
  ncol = 3, 
  rel_widths = c(0.2, 0.6, 0.2), # Adjust these values to control spacing
  labels = c("", "I", ""), 
  label_size = 20
)

#############################################
############ FINAL FIG --> FIG 6 ############
#############################################

Fig6_final <- plot_grid(fig_up, fig_mid, fig_bottom, 
                        nrow = 3, 
                        rel_heights = c(0.2, 0.4, 0.3), 
                        align = "v") # Align all vertically

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig6_complete.png"
ggsave(final_figure_path, Fig6_final, width = 250, height = 330, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig6_complete.pdf"
ggsave(final_figure_path, Fig6_final, width = 250, height = 330, units = "mm", device = "pdf")
