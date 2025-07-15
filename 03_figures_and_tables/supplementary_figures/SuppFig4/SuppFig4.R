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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_A.RData")
# Stuff
formula <- y ~ x
method1 <- paste0("ASE_PTC_NMD_triggering_0.2")
method2 <- paste0("endogenous_NMD_Consensus")
# Check sample size
m1_ss <- sum(!is.na(input_figureA[,method1]))
m2_ss <- sum(!is.na(input_figureA[,method2]))
if (m1_ss >= m2_ss) {
    sample_size <- m2_ss
} else {
    sample_size <- m1_ss
}

# Plot
plot_A <- ggplot(data = input_figureA, mapping = aes(x = eval(parse(text = method1)), y = eval(parse(text = method2)) )) +
    geom_point(alpha = 0.15, color = "black") +
    geom_smooth(method = "lm", formula = formula, se = TRUE, size = 1, color = "#69b3a2") +
    labs(x = "ASE iNMDeff", y = "ETG iNMDeff", title = paste0("TCGA pan-cancer, n = ",sample_size)) +
    coord_cartesian(xlim = c(-3,3), ylim = c(-3,3)) +
    ggplot_theme() + 
    scale_colour_discrete(na.translate = FALSE) +
    stat_cor(p.accuracy = NULL, r.accuracy = 0.01, size = 5, label.x = -2.75, label.y = 2)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_B.RData")
# Stuff
formula <- y ~ x
method1 <- paste0("ASE_PTC_NMD_triggering_0.2")
method2 <- paste0("endogenous_NMD_Consensus")
# Check sample size
m1_ss <- sum(!is.na(input_figureB[,method1]))
m2_ss <- sum(!is.na(input_figureB[,method2]))
if (m1_ss >= m2_ss) {
    sample_size <- m2_ss
} else {
    sample_size <- m1_ss
}

# Plot
plot_B <- ggplot(data = input_figureB, mapping = aes(x = eval(parse(text = method1)), y = eval(parse(text = method2)) )) +
    geom_point(alpha = 0.15, color = "black") +
    geom_smooth(method = "lm", formula = formula, se = TRUE, size = 1, color = "#69b3a2") +
    labs(x = "ASE iNMDeff", y = "ETG iNMDeff", title = paste0("GTex pan-tissue, n = ",sample_size)) +
    coord_cartesian(xlim = c(-4,4), ylim = c(-4,4)) +
    ggplot_theme() + 
    scale_colour_discrete(na.translate = FALSE) +
    stat_cor(p.accuracy = NULL, r.accuracy = 0.01, size = 5, label.x = -3.75, label.y = 3)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_C.RData")

plot_C <- ggplot(data = input_figureC, aes(x = tissue, y = R)) +
    geom_bar(stat = "identity", linewidth = 0.7, fill = "skyblue", width = 0.6) + coord_flip(ylim = c(-0.5,0.5)) +
    facet_grid(. ~ database) +
    labs(y = "Pearson correlation (R)", x = "", title = "") +
    geom_errorbar( aes(ymin = R_conf_int_low, ymax = R_conf_int_high), width = 0.2) +
    ggplot_theme_bw() + scale_fill_brewer(palette = "Pastel1") +
    theme(legend.position='top',
            strip.background = element_blank(),
            axis.text.y = element_text(colour = "black", size = 9)) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
            position = position_dodge(width = 1), size = 9, 
            hjust = -0.1, color = "black", vjust = 0.75)

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_D.RData")

input_figureD$database <- ifelse(input_figureD$database == "GTEx","GTex",input_figureD$database )

plot_D <- ggplot(data = input_figureD, aes(x = tissue, y = R)) +
    geom_bar(stat = "identity", linewidth = 0.7, fill = "skyblue", width = 0.6) + coord_flip(ylim = c(-0.5,0.5)) +
    facet_grid(. ~ database) +
    labs(y = "Pearson correlation (R)", x = "", title = "") +
    geom_errorbar( aes(ymin = R_conf_int_low, ymax = R_conf_int_high), width = 0.2) +
    ggplot_theme_bw() + scale_fill_brewer(palette = "Pastel1") +
    theme(legend.position='top',
        strip.background = element_blank(),
        axis.text.y = element_text(colour = "black", size = 9)) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
            position = position_dodge(width = 1), size = 9, 
            hjust = -0.1, color = "black", vjust = 0.75)


################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_E.RData")

plot_E <- ggplot(data = input_figureE, aes(x = bins, y = NMDeff, fill = factor(NMD_method))) +
        geom_boxplot(color="black", alpha=0.5) + ylim(c(-3,3)) +
        labs(x = "Percentiles (%), ordered by ASE method", y = "iNMDeff", fill = "NMD method") +
        ggtitle("TCGA") +
        scale_x_discrete(labels = paste0(1:10,"0")) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        theme_bw() + ggplot_theme() +
        theme(axis.title.x = element_text(size=11)) +
        scale_fill_brewer(palette = "Accent", labels = c("ETG","ASE"), direction = -1) +
        guides(fill = guide_legend(override.aes = list(size = 8)))

################################
############ PLOT F ############
################################

input_figureF <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/panel_F.RData")

plot_F <- ggplot(data = input_figureF, aes(x = bins, y = NMDeff, fill = factor(NMD_method))) +
        geom_boxplot(color="black", alpha=0.5) + ylim(c(-3,3)) +
        labs(x = "Percentiles (%), ordered by ASE method", y = "iNMDeff", fill = "NMD method") +
        ggtitle("GTex") +
        scale_x_discrete(labels = paste0(1:10,"0")) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        theme_bw() + ggplot_theme() +
        theme(axis.title.x = element_text(size=11)) +
        scale_fill_brewer(palette = "Accent", labels = c("ETG","ASE"), direction = -1) +
        guides(fill = guide_legend(override.aes = list(size = 8)))

########################################################
############ PANEL UP --> PLOTS A and E ################
########################################################

fig_up <- plot_grid(plot_A, plot_E, nrow = 1, ncol = 2, labels = c("A", "E"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

#########################################################
############ PANEL MID --> PLOTS B and F ################
#########################################################

fig_mid <- plot_grid(plot_B,plot_F, nrow = 1, ncol = 2, labels = c("B", "F"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

########################################################
############ PANEL BOTTOM --> PLOTS C-D ################
########################################################

fig_bottom <- plot_grid(plot_C, plot_D, nrow = 1, ncol = 2, labels = c("C","D"), label_size = 20, vjust = 1)

##################################################
############ FINAL FIG --> SUPP FIG 2 ############
##################################################

SuppFig_final <- plot_grid(fig_up, fig_mid, fig_bottom, nrow = 3, rel_heights = c(0.175,0.175,0.45))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig4_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 325, units = "mm")
