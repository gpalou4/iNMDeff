library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(corrplot)
library(ggcorrplot)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA_1 <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/panel_A.RData")
input_figureA_1$type <- ifelse(input_figureA_1$type == "NMD_efficiency","Observed",input_figureA_1$type)

# Plot
plot_A_1 <- ggplot(data = input_figureA_1, aes(x = tissue, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","Observed")))) +
  geom_boxplot() + coord_flip(ylim = c(-2.5,2.5)) + #coord_cartesian(ylim = c(-5,5)) +
  facet_wrap(. ~ NMD_method) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  labs(y = "iNMDeff", x = "", title = "", fill = "iNMDeff") +
  ggplot_theme_bw() + scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "top",
        legend.title = element_text(size = 12),
        strip.background = element_blank(),
        legend.text = element_text(size = 11))

plot_A <- annotate_figure(plot_A_1, top = text_grob("TCGA",
                            size = 20, face = "plain"))

################################
############ PLOT B ############
################################

input_figureB_1 <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/panel_B.RData")
input_figureB_1$type <- ifelse(input_figureB_1$type == "NMD_efficiency","Observed",input_figureB_1$type)

# Plot
plot_B_1 <- ggplot(data = input_figureB_1, aes(x = tissue, y = NMD_efficiency, fill = factor(type, levels = c("Randomization","Observed")))) +
  geom_boxplot() + coord_flip(ylim = c(-2.5,2.5)) + #coord_cartesian(ylim = c(-5,5)) +
  facet_wrap(. ~ NMD_method) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
  labs(y = "iNMDeff", x = "", title = "", fill = "iNMDeff") +
  ggplot_theme_bw() + scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "top",
        legend.text = element_text(size = 11),
        strip.background = element_blank(),
        legend.title = element_text(size = 12))

plot_B <- annotate_figure(plot_B_1, top = text_grob("GTex",
                            size = 20, face = "plain")) 

########################################################
############ PANEL UP --> PLOTS A-B ####################
########################################################

# fig_up <- plot_grid(plot_A, plot_B, labels = c("A","B"), nrow = 1, ncol = 2, label_size = 20)

fig_up <- ggarrange(plot_A,plot_B,
                labels = c("A","B","",""),
                font.label = list(size = 20),
                ncol=2, common.legend = TRUE)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/panel_C.RData")

input_figureC$database <- ifelse(input_figureC$database == "GTEx","GTex","TCGA")

# input_figureC$NMD_geneset

# Plot
plot_C_top <- input_figureC %>% filter(NMD_method == "ASE") %>% 
        ggplot(aes(x = NMD_geneset, y = median_SD_of_medians_diff_mean, fill = factor(controls))) +
        geom_bar(stat = "identity", linewidth = 0.6, color = "black") + coord_flip() + #coord_cartesian(ylim = c(-5,5)) +
        # Add error bars for confidence intervals
        geom_errorbar(aes(ymin = perc5_SD_of_medians_diff_mean, ymax = perc95_SD_of_medians_diff_mean), 
                        width = 0.2) +  # Adjust width as needed
        guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
        facet_grid(NMD_method ~ database, scale = "free") +
        ylab("") + ggtitle("") + xlab("") +
        ggplot_theme() +
         scale_x_discrete(labels = c("PTC NMD-evading\n0.2" = "NMD-evading PTCs",
                                        "Synonymous 0.2" = "Synonymous",
                                        "PTC NMD-\ntriggering 0.2" = "NMD-triggering PTCs")) +
        scale_fill_brewer(palette = "Pastel1") + scale_color_viridis(discrete=TRUE, option="inferno", direction = -1) +
        theme(plot.title = element_text(hjust = 0.5, size = 35),
                plot.margin = unit(c(-1.25, 1, -0.25, 1), "cm"),
                axis.text.y = element_text(size = 11),
                legend.text = element_text(size = 10),
                legend.position='none') +
        geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
                                position = position_dodge(), size = 9, color = "black", hjust = -0.2, vjust = 0.5)

plot_C_bottom <- input_figureC %>% filter(NMD_method == "ETG") %>% 
        ggplot(aes(x = NMD_geneset, y = median_SD_of_medians_diff_mean, fill = factor(controls))) +
        geom_bar(stat = "identity", linewidth = 0.6, color = "black") + coord_flip() +
        facet_grid(NMD_method ~ database, scale = "free") +
        geom_errorbar(aes(ymin = perc5_SD_of_medians_diff_mean, ymax = perc95_SD_of_medians_diff_mean), 
                        width = 0.2) +  # Adjust width as needed
        guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
        ylab("Inter-tissue iNMDeff variability deviation") + ggtitle("") + xlab("") +
        ggplot_theme() + 
        scale_x_discrete(labels = c("RandomGenes\nwithout NMD\nfeatures" = "RandomGenes\nw/o NMD features",
                                "RandomGenes\nwith NMD\nfeatures" = "RandomGenes\nwith NMD features")) +
        scale_fill_brewer(palette = "Pastel1") + scale_color_viridis(discrete=TRUE, option="inferno", direction = -1) +
        theme(plot.title = element_text(hjust = 0.5, size = 35),
                axis.text.y = element_text(size = 9, lineheight = 0.7),
                plot.margin = unit(c(-1.25, 1, 0.25, 1), "cm"),
                legend.text = element_text(size = 11),
                legend.title = element_text(size = 12),
                legend.position='bottom') +
        geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
                                position = position_dodge(), size = 9, color = "black", hjust = -0.2, vjust = 0.5)

plot_C_tmp <- plot_grid(plot_C_top,plot_C_bottom, nrow = 2,  align = "hv",
                        rel_heights = c(0.3, 0.7), axis = "right")

empty_plot <- plot_grid(NULL, NULL, labels = "")

plot_C <- plot_grid(empty_plot, plot_C_tmp, empty_plot, nrow = 3, labels = c("C","",""), label_size = 20, 
                        rel_heights = c(0.04,1,0.04), vjust = 1, label_y = 1.2)

##########################################################
############ PANEL bottom --> PLOTS C ####################
##########################################################

fig_bottom <- plot_grid(empty_plot,plot_C,empty_plot, nrow = 1, ncol = 3, label_size = 20, 
                        rel_widths = c(0.2,0.6,0.2))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.65,0.35)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig7_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 275, height = 350, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig7_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 275, height = 350, units = "mm")
