library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(corrplot)
library(ggcorrplot)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(scales)
library(ggh4x)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

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
############ PLOT A ############
################################

# Data
input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig26/panel_A.RData")
input_figureA$gene_name <- factor(input_figureA$gene_name, levels = input_figureA$gene_name)
# input_figureC$gene_name <- paste0(input_figureC$gene_name," (n=",input_figureC$n_non,")")
input_figureA$n_non_total <- input_figureA$n_non_High + input_figureA$n_non_Low
plot_A <- input_figureA %>%
      filter(n_non_total >= 10) %>%
        ggplot(aes(x = gene_name, 
                      y = dNdScv_diff, fill = factor(dNdScv_group))) +
          geom_bar(stat = "identity", position=position_dodge(width=0.9)) + #coord_flip() +
          labs(title = "", x = "", y = "dNdS difference", fill = "dNdS mean") + 
          theme_bw() + ggplot_theme() + ylim(-6.5,10) +
          scale_fill_brewer(palette = "Set2", labels = c("pos_sel" = " >1 (Positive selection)", "neg_sel" = " <1 (No positive selection)"), direction = -1) +
          theme( legend.position = "bottom",
                  # plot.title = element_text(size = 12, hjust = 0.5),
                  legend.text = element_text(colour = "black", size = 10),
                  plot.margin = unit(c(1,0.5,1,0.5), "cm"),
                  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 9)
                ) +
          geom_text(aes(label = n_non_total), position = position_dodge(width=0.9), 
              size = 1.75, color = "black", hjust = 0.5, vjust = 0)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig26/panel_B.RData")

# Change
input_figureB$mutation_type <- as.character(input_figureB$mutation_type)
input_figureB$mutation_type <- ifelse(input_figureB$mutation_type == "wmis_cv", "Missense",input_figureB$mutation_type)
input_figureB$mutation_type <- ifelse(input_figureB$mutation_type == "wnon_cv", "Nonsense",input_figureB$mutation_type)
input_figureB$mutation_type <- ifelse(input_figureB$mutation_type == "wind_cv", "Indel",input_figureB$mutation_type)
table(input_figureB$mutation_type)

list_plots <- list() 
combinations <- combn(names(table(input_figureB$NMDeff)), 2, simplify = FALSE)
for (NMD_method_char in c("ASE")) {
  for (mutation_type_char in c("Nonsense","Missense")) {
    if ( mutation_type_char == "Nonsense" ) {
      if (NMD_method_char == "ASE") { 
            ylim <- c(0,6)
      } else if (NMD_method_char == "END") {
            ylim <- c(0,5)
      }
    } else if ( mutation_type_char == "Missense" ) {
      if (NMD_method_char == "ASE") { 
        ylim <- c(0,3.5)
      } else if (NMD_method_char == "END") {
        ylim <- c(0,3)
      }
    }
    p <- input_figureB %>%
        filter(gene_type == "TSG") %>%
        mutate(NMD_variant_type = factor(NMD_variant_type, levels = c("NMD-triggering","NMD-evading"))) %>%
        # filter(percentile %in% c(20,30,40,50)) %>%
        filter(percentile %in% c(50)) %>%
        filter(NMD_method %in% NMD_method_char) %>%
        filter(mutation_type %in% mutation_type_char) %>%
          ggplot(aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                  geom_violin() + scale_fill_brewer(palette = "Pastel1") + xlab("") + ylab("dNdS ratio") +
                  geom_jitter(aes(fill = factor(NMDeff)), 
                    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                    alpha = 0.1, size = 2) + xlab("") +
                  labs(fill = "ASE iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                  ggtitle(paste0(mutation_type_char," in TSG")) +
                  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                  facet_grid(. ~ NMD_variant_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                  geom_boxplot(width=0.3, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                  ggplot_theme_bw() +
                  theme(strip.text = element_text(size = 9)) +
                  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
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
wilcox_test <- function(input_figure, mutation_type_char, gene_type_char, NMD_variant_type_char,
                        NMD_method_char, percentile_char) {
  input_figure_filt <- input_figure %>% 
            filter(mutation_type == mutation_type_char & gene_type == gene_type_char &
                  percentile == percentile_char & NMD_method == NMD_method_char & 
                  NMD_variant_type == NMD_variant_type_char)
  # Subset data for 'High' and 'Low' NMDeff
  high_group <- subset(input_figure_filt, NMDeff == "High")
  colnames(high_group)[colnames(high_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_high"
  low_group <- subset(input_figure_filt, NMDeff == "Low")
  colnames(low_group)[colnames(low_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_low"
  # Merge subsets by 'gene_name'
  paired_data <- merge(high_group[,c("gene_name","dNdScv_ratio_high")], low_group[,c("gene_name","dNdScv_ratio_low")], by = c("gene_name"))       
  wilcox_res <- wilcox.test(paired_data$dNdScv_ratio_high, paired_data$dNdScv_ratio_low, paired = TRUE, 
                            method = "wilcox.test", alternative = "greater")
  #wilcox.test(dNdScv_ratio ~ NMDeff, data = input_figureB_filt, paired = FALSE, alternative = "greater")
  return(wilcox_res$p.value)
}

# Wilcoxon tests with pairing data
pvalue1 <- wilcox_test(input_figureB, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "ASE", percentile_char = 50)
pvalue2 <- wilcox_test(input_figureB, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "ASE", percentile_char = 50)
legend_B <- get_legend(list_plots[[1]])
plot_B_1 <- ggplot_build(list_plots[[1]] + guides(fill = FALSE))
plot_B_1$data[[5]][,"annotation"] <- as.character(plot_B_1$data[[5]][,"annotation"])
plot_B_1$data[[5]][,"annotation"][1:3] <- rep(round(pvalue1,2),3)
plot_B_1$data[[5]][,"annotation"][4:6] <- rep(round(pvalue2,2),3)
plot_B_1$data[[5]][,"annotation"] <- as.factor(plot_B_1$data[[5]][,"annotation"])
plot_B_1 <- ggplot_gtable(plot_B_1)
plot_B_1 <- ggplotify::as.ggplot(plot_B_1)
pvalue3 <- wilcox_test(input_figureB, mutation_type_char = "Missense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "ASE", percentile_char = 50)
pvalue4 <- wilcox_test(input_figureB, mutation_type_char = "Missense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "ASE", percentile_char = 50)
# plot_B_2 <- ggplot_Cuild(list_plots[[2]])
plot_B_2 <- ggplot_build(list_plots[[2]] + guides(fill = FALSE))
plot_B_2$data[[5]][,"annotation"] <- as.character(plot_B_2$data[[5]][,"annotation"])
plot_B_2$data[[5]][,"annotation"][1:3] <- rep(round(pvalue3,2),3)
plot_B_2$data[[5]][,"annotation"][4:6] <- rep(round(pvalue4,2),3)
plot_B_2$data[[5]][,"annotation"] <- as.factor(plot_B_2$data[[5]][,"annotation"])
plot_B_2 <- ggplot_gtable(plot_B_2)
plot_B_2 <- ggplotify::as.ggplot(plot_B_2)

plot_B <- ggarrange(plot_B_1,plot_B_2, labels = c("",""), legend.grob = legend_B,
                font.label = list(size = 20), #widths = c(0.45,0.28,0.28),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
# plot_B <- annotate_figure(plot_B, top = text_grob("Tumor Supressor Genes", 
#                             size = 16, face = "plain")) #+ ggplot_theme_bw()

################################
############ PLOT C ############
################################

# Data
input_figureC <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig26/panel_C.RData")
input_figureC$gene_name <- factor(input_figureC$gene_name, levels = input_figureC$gene_name)
input_figureC$n_non_total <- input_figureC$n_non_High + input_figureC$n_non_Low

plot_C <- input_figureC %>%
      filter(n_non_total >= 10) %>%
      ggplot(aes(x = gene_name, 
                    y = dNdScv_diff, fill = factor(dNdScv_group))) +
        geom_bar(stat = "identity", position=position_dodge()) + #coord_flip() +
        labs(title = "", x = "", y = "dNdS difference", fill = "dNdS mean") + 
        theme_bw() + ggplot_theme() +
        scale_fill_brewer(palette = "Set2", labels = c("pos_sel" = " > 1 (Positive Selection)", "neg_sel" = " < 1 (Negative Selection)"), direction = -1) +
        theme( legend.position = "bottom",
                legend.text = element_text(colour = "black", size = 9),
                plot.margin = unit(c(0, 1, 0, 1), "cm"),
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
              ) +
          geom_text(aes(label = n_non_total), position = position_dodge(width=0.9), 
              size = 2, color = "black", hjust = 0.5, vjust = 0)

################################
############ PLOT D ############
################################

# Data
input_figureD <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig26/panel_D.RData")
# Change
input_figureD$mutation_type <- as.character(input_figureD$mutation_type)
input_figureD$mutation_type <- ifelse(input_figureD$mutation_type == "wmis_cv", "Missense",input_figureD$mutation_type)
input_figureD$mutation_type <- ifelse(input_figureD$mutation_type == "wnon_cv", "Nonsense",input_figureD$mutation_type)
input_figureD$mutation_type <- ifelse(input_figureD$mutation_type == "wind_cv", "Indel",input_figureD$mutation_type)

table(input_figureD$mutation_type)
list_plots <- list() 
for (NMD_method_char in c("END","ASE")) {
  for (mutation_type_char in c("Nonsense","Missense")) {
    if ( mutation_type_char == "Nonsense" ) {
      if (NMD_method_char == "ASE") { 
            ylim <- c(0,6)
      } else if (NMD_method_char == "END") {
            ylim <- c(0,3)
      }
    } else if ( mutation_type_char == "Missense" ) {
      if (NMD_method_char == "ASE") { 
        ylim <- c(0,4)
      } else if (NMD_method_char == "END") {
        ylim <- c(0,5)
      }
    }
    p <- input_figureD %>%
        filter(gene_type == "oncogene") %>%
        # filter(percentile %in% c(20,30,40,50)) %>%
        mutate(NMD_variant_type = factor(NMD_variant_type, levels = c("NMD-triggering","NMD-evading"))) %>%
        filter(percentile %in% c(50)) %>%
        filter(NMD_method %in% NMD_method_char) %>%
        filter(mutation_type %in% mutation_type_char) %>%
          ggplot(aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                  geom_violin() + scale_fill_brewer(palette = "Pastel1") + xlab("") + ylab("dNdS ratio") +
                  geom_jitter(aes(fill = factor(NMDeff)), 
                    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                    alpha = 0.1, size = 2) + xlab("") +
                  labs(fill = "ETG iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                  ggtitle(paste0(mutation_type_char," in OGs")) +
                  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                  facet_grid(. ~ NMD_variant_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                  # ylim(ylim)+
                  geom_boxplot(width=0.3, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                  ggplot_theme_bw() +
                  theme(strip.text = element_text(size = 9))+
                  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
                stat_compare_means(aes(group = NMDeff),
                                  method.args = list(alternative = "greater"),
                                  comparisons = combinations, size = 5, 
                                  label.y = c(1.25), na.rm = TRUE, 
                                  label = "p.format", method = "wilcox.test", hide.ns = FALSE)
      list_plots[[length(list_plots) + 1]] <- p
  }
}

# Wilcoxon tests with pairing data
pvalue1 <- wilcox_test(input_figureD, mutation_type_char = "Nonsense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
pvalue2 <- wilcox_test(input_figureD, mutation_type_char = "Nonsense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
legend_D <- get_legend(list_plots[[1]])
plot_D_1 <- ggplot_build(list_plots[[1]] + guides(fill = FALSE))
plot_D_1$data[[5]][,"annotation"] <- as.character(plot_D_1$data[[5]][,"annotation"])
plot_D_1$data[[5]][,"annotation"][1:3] <- rep(round(pvalue1,2),3)
plot_D_1$data[[5]][,"annotation"][4:6] <- rep(round(pvalue2,2),3)
plot_D_1$data[[5]][,"annotation"] <- as.factor(plot_D_1$data[[5]][,"annotation"])
plot_D_1 <- ggplot_gtable(plot_D_1)
plot_D_1 <- ggplotify::as.ggplot(plot_D_1)
pvalue3 <- wilcox_test(input_figureD, mutation_type_char = "Missense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
pvalue4 <- wilcox_test(input_figureD, mutation_type_char = "Missense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
plot_D_2 <- ggplot_build(list_plots[[2]] + guides(fill = FALSE))
plot_D_2$data[[5]][,"annotation"] <- as.character(plot_D_2$data[[5]][,"annotation"])
plot_D_2$data[[5]][,"annotation"][1:3] <- rep(round(pvalue3,2),3)
plot_D_2$data[[5]][,"annotation"][4:6] <- rep(round(pvalue4,2),3)
plot_D_2$data[[5]][,"annotation"] <- as.factor(plot_D_2$data[[5]][,"annotation"])
plot_D_2 <- ggplot_gtable(plot_D_2)
plot_D_2 <- ggplotify::as.ggplot(plot_D_2)

plot_D <- ggarrange(plot_D_1,plot_D_2, labels = c("D",""), legend.grob = legend_D,
                font.label = list(size = 20), #widths = c(0.45,0.28,0.28),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
# plot_C <- annotate_figure(plot_C, top = text_grob("Oncogenes",  
#                             size = 16, face = "plain"), ) #+ ggplot_theme_bw()

################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig26/panel_E.RData")
# Change
input_figureE$mutation_type <- as.character(input_figureE$mutation_type)
input_figureE$mutation_type <- ifelse(input_figureE$mutation_type == "wmis_cv", "Missense",input_figureE$mutation_type)
input_figureE$mutation_type <- ifelse(input_figureE$mutation_type == "wnon_cv", "Nonsense",input_figureE$mutation_type)
input_figureE$mutation_type <- ifelse(input_figureE$mutation_type == "wind_cv", "Indel",input_figureE$mutation_type)
table(input_figureE$mutation_type)

list_plots <- list() 
for (NMD_method_char in c("ASE","END")) {
  for (mutation_type_char in c("Nonsense","Missense")) {
    if ( mutation_type_char == "Nonsense" ) {
      if (NMD_method_char == "ASE") { 
            ylim <- c(0,6)
      } else if (NMD_method_char == "END") {
            ylim <- c(0,3)
      }
    } else if ( mutation_type_char == "Missense" ) {
      if (NMD_method_char == "ASE") { 
        ylim <- c(0,4)
      } else if (NMD_method_char == "END") {
        ylim <- c(0,5)
      }
    }
    p <- input_figureE %>%
        filter(gene_type == "oncogene") %>%
        mutate(NMD_variant_type = factor(NMD_variant_type, levels = c("NMD-triggering","NMD-evading"))) %>%
        # filter(percentile %in% c(20,30,40,50)) %>%
        filter(percentile %in% c(50)) %>%
        filter(NMD_method %in% NMD_method_char) %>%
        filter(mutation_type %in% mutation_type_char) %>%
          ggplot(aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                  geom_violin() + scale_fill_brewer(palette = "Pastel1") + xlab("") + ylab("dNdS ratio") +
                  geom_jitter(aes(fill = factor(NMDeff)), 
                    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                    alpha = 0.1, size = 2) + xlab("") +
                  labs(fill = "ASE iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                  ggtitle(paste0(mutation_type_char," in OGs")) +
                  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                  facet_grid(. ~ NMD_variant_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                  # ylim(ylim)+
                  geom_boxplot(width=0.3, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                  ggplot_theme_bw() +
                  theme(strip.text = element_text(size = 9)) +
                  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
                stat_compare_means(aes(group = NMDeff),
                                  method.args = list(alternative = "greater"),
                                  comparisons = combinations, size = 5, 
                                  label.y = c(1.25), na.rm = TRUE, 
                                  label = "p.format", method = "wilcox.test", hide.ns = FALSE)
      list_plots[[length(list_plots) + 1]] <- p
  }
}

# Wilcoxon tests with pairing data
pvalue1 <- wilcox_test(input_figureE, mutation_type_char = "Nonsense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "ASE", percentile_char = 50)
pvalue2 <- wilcox_test(input_figureE, mutation_type_char = "Nonsense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "ASE", percentile_char = 50)
legend_E <- get_legend(list_plots[[1]])
plot_E_1 <- ggplot_build(list_plots[[1]] + guides(fill = FALSE))
plot_E_1$data[[5]][,"annotation"] <- as.character(plot_E_1$data[[5]][,"annotation"])
plot_E_1$data[[5]][,"annotation"][1:3] <- rep(round(pvalue1,2),3)
plot_E_1$data[[5]][,"annotation"][4:6] <- rep(round(pvalue2,2),3)
plot_E_1$data[[5]][,"annotation"] <- as.factor(plot_E_1$data[[5]][,"annotation"])
plot_E_1 <- ggplot_gtable(plot_E_1)
plot_E_1 <- ggplotify::as.ggplot(plot_E_1)
pvalue3 <- wilcox_test(input_figureE, mutation_type_char = "Missense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "ASE", percentile_char = 50)
pvalue4 <- wilcox_test(input_figureE, mutation_type_char = "Missense", gene_type_char = "oncogene",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "ASE", percentile_char = 50)
plot_E_2 <- ggplot_build(list_plots[[2]] + guides(fill = FALSE))
plot_E_2$data[[5]][,"annotation"] <- as.character(plot_E_2$data[[5]][,"annotation"])
plot_E_2$data[[5]][,"annotation"][1:3] <- rep(round(pvalue3,2),3)
plot_E_2$data[[5]][,"annotation"][4:6] <- rep(round(pvalue4,2),3)
plot_E_2$data[[5]][,"annotation"] <- as.factor(plot_E_2$data[[5]][,"annotation"])
plot_E_2 <- ggplot_gtable(plot_E_2)
plot_E_2 <- ggplotify::as.ggplot(plot_E_2)

plot_E <- ggarrange(plot_E_1,plot_E_2, labels = c("E",""), legend.grob = legend_E,
                font.label = list(size = 20), #widths = c(0.45,0.28,0.28),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
# plot_E <- annotate_figure(plot_C, top = text_grob("Oncogenes",  
#                             size = 16, face = "plain"), ) #+ ggplot_theme_bw()

######################################################
############ PANEL UP --> PLOTS A-B ##################
######################################################

fig_up <- plot_grid(plot_B, plot_A, labels = c("A","B"), label_size = 20,
        ncol = 2, nrow = 1, rel_widths = c(0.5,0.5))

######################################################
############ PANEL MID --> PLOTS C-D ##################
######################################################

fig_mid <- plot_grid(plot_C, plot_D, labels = c("C","D"), label_size = 20,
        ncol = 2, nrow = 1, rel_widths = c(0.5,0.5))
        
#################################################
############ PANEL BOTTOM --> PLOT E ##########
#################################################

fig_bottom <- plot_grid(plot_E, ncol = 2, nrow = 1, rel_widths = c(0.5,0.5))

#############################################
############ FINAL FIG --> FIG 5 ############
#############################################

SuppFig_final <- plot_grid(fig_up, fig_mid, fig_bottom, nrow = 3, rel_heights = c(0.5,0.5)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig26_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 250, height = 220, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig26_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 250, height = 220, units = "mm")
