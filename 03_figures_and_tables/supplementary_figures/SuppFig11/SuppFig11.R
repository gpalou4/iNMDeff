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
library(purrr)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

####################################
############ PLOT A ################
####################################

# Data
input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/panel_A.RData")
# Colors
my_colors <- setNames(c("#BEAED4","#aed4d3"), c("GBM","LGG"))

# Plots A1-5
plot_A_1 <- ggplot(data = input_figureA, aes(x = cancer_type_strat, 
                    y = ETG, fill = factor(cancer_type))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-1.5,1.5)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ETG", x = "", y = "iNMDeff", fill = "") +
        theme_bw() +  ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
              plot.margin = unit(c(0,0,0.6,0), "cm"),
                plot.title = element_text(size = 12, hjust = 0.5))

plot_A_2 <- ggplot(data = input_figureA, aes(x = cancer_type_strat, 
                    y = log2(RBFOX3), fill = factor(cancer_type))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip() +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        #facet_grid(. ~ NMD_method , scale = "free") +
        scale_fill_manual(
          values = c("LGG" = "#89a7ea", "GBM" = "#89a7ea")  # Assign colors to genes
        ) +      
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "RBFOX3 (neuron)", x = "", y = expression(log[10](TPM)), fill = "") +
        theme_bw() + ggplot_theme() +
        theme(legend.position = "top",
                plot.margin = unit(c(0,0.1,0.5,0), "cm"),
                plot.title = element_text(size = 12, hjust = 0.5),
                axis.text.y = element_blank())

plot_A_3 <- ggplot(data = input_figureA, aes(x = cancer_type_strat, 
                    y = log2(AQP4), fill = factor(cancer_type))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(3,12)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        #facet_grid(. ~ NMD_method , scale = "free") +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "AQP4 (glia)", x = "", y = expression(log[10](TPM)), fill = "") +
        theme_bw() + ggplot_theme() +
        scale_fill_manual(
          values = c("LGG" = "#e4b56e", "GBM" = "#e4b56e")  # Assign colors to genes
        ) +
            # scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
              plot.margin = unit(c(0,0.1,0.5,0), "cm"),
              plot.title = element_text(size = 12, hjust = 0.5),
              axis.text.y = element_blank())

plot_A_1_2_3 <- ggarrange(plot_A_1,plot_A_2,plot_A_3, labels = c("A","",""), font.label = list(size = 20),
                widths = c(0.45,0.28,0.28),
                ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

# Reshape data to long format
input_long <- input_figureA %>%
  pivot_longer(
    cols = c(AQP4, RBFOX3),  # Columns to reshape
    names_to = "cell_type",  # New column for cell type
    values_to = "value"      # New column for values
  )

# Create the plot
plot_A_4 <- ggplot(data = input_long, aes(x = scale(ETG), y = log10(value), color = cell_type)) +
  geom_smooth(method = lm, se = FALSE, alpha = 1, size = 1) +
  geom_point(size = 1, alpha = 0.65) +
  scale_color_manual(
    values = c("AQP4" = "#e4b56e", "RBFOX3" = "#89a7ea")  # Assign colors to genes
  ) +
  labs(
    title = "",
    color = "Gene",
    x = "ETG iNMDeff",
    y = expression(log[10](TPM))
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  ggplot_theme() +
  ggpubr::stat_cor(aes(group = cell_type), method = "pearson", label.y = c(-2,-1.5), label.x = -3, size = 4) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(1.25, 0.5, 1.25, 0), "cm")
  )

plot_A <- plot_grid(plot_A_1_2_3, plot_A_4, nrow = 1, ncol = 2,
                        rel_widths = c(0.72,0.26), vjust = 1)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
ggsave(final_figure_path, plot_A, width = 300, height = 125, units = "mm") 

##################################################
############ PANEL TOP --> PLOT A ################
##################################################

fig_up <- plot_A

####################################
############ PLOT B ################
####################################

# Data
input_figureB <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/panel_B.RData")
# Colors
my_colors <- setNames(rep(c("#BEAED4"),,14), unique(input_figureB$acronyms))

# Plots B1-5
plot_B_1 <- ggplot(data = input_figureB, aes(x = acronyms, 
                    y = ETG, fill = factor(acronyms))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-3,1)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ETG", x = "", y = "iNMDeff", fill = "") +
        theme_bw() +  ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.title = element_text(size = 12, hjust = 0.5))

my_colors <- setNames(rep(c("#89a7ea"),,14), unique(input_figureB$acronyms))

plot_B_2 <- ggplot(data = input_figureB, aes(x = acronyms, 
                    y = log2(RBFOX3), fill = factor(acronyms))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(0,8)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        #facet_grid(. ~ NMD_method , scale = "free") +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "RBFOX3 (neuron)", x = "", y = expression(log[10](TPM)), fill = "") +
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.margin = unit(c(0.25,0.1,0.15,0), "cm"),
                plot.title = element_text(size = 12, hjust = 0.5),
                        axis.text.y = element_blank())

my_colors <- setNames(rep(c("#e4b56e"),,14), unique(input_figureB$acronyms))

plot_B_3 <- ggplot(data = input_figureB, aes(x = acronyms, 
                    y = log2(AQP4), fill = factor(acronyms))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(0,7)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        #facet_grid(. ~ NMD_method , scale = "free") +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "AQP4 (glia)", x = "", y = expression(log[10](TPM)), fill = "") +
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.margin = unit(c(0.25,0.1,0.15,0), "cm"),
                plot.title = element_text(size = 12, hjust = 0.5),
                        axis.text.y = element_blank())

plot_B_1_2_3 <- ggarrange(plot_B_1,plot_B_2,plot_B_3,labels = c("B","","",""), font.label = list(size = 20),
                widths = c(0.45,0.28,0.28), vjust = 1,
                ncol=3, nrow=1, common.legend = FALSE, legend="none")

# Reshape data to long format
input_long <- input_figureB %>%
  pivot_longer(
    cols = c(AQP4, RBFOX3),  # Columns to reshape
    names_to = "cell_type",  # New column for cell type
    values_to = "value"      # New column for values
  )

# Create the plot
plot_B_4 <- ggplot(data = input_long, aes(x = scale(ETG), y = log10(value), color = cell_type)) +
  geom_smooth(method = lm, se = FALSE, alpha = 1, size = 1) +
  geom_point(size = 1, alpha = 0.35) +
  scale_color_manual(
    values = c("AQP4" = "#e4b56e", "RBFOX3" = "#89a7ea")  # Assign colors to genes
  ) +
  labs(
    title = "",
    color = "Gene",
    x = "ETG iNMDeff",
    y = expression(log[10](TPM))
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  ggplot_theme() +
  ggpubr::stat_cor(aes(group = cell_type), method = "pearson", label.y = c(-2,-1.5), label.x = -3, size = 4) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(1.25, 0.5, 1.25, 0), "cm")
  )

plot_B <- plot_grid(plot_B_1_2_3, plot_B_4, nrow = 1, ncol = 2,
                        rel_widths = c(0.72,0.26), vjust = 1)



final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
ggsave(final_figure_path, plot_B, width = 300, height = 125, units = "mm") 


#################################################
############ PANEL MID --> PLOT B ###############
#################################################

fig_middle <- plot_B

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/panel_C.RData")
input_figureC$dataset <- gsub("_"," ", input_figureC$dataset)

plot_C <- ggplot(data = input_figureC, aes(x = factor(type_var), y = correlation, fill = factor(randomization))) +
    geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(-0.1,1)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    # geom_jitter(aes(fill = factor(randomization)), 
    #             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
    #             alpha = 0.01, size = 1) +
    facet_wrap(. ~ dataset) +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    labs(y = "Spearman correlation", x = "", title = "ASE PTC-NMDeff", fill = "") +
    ggplot_theme_bw() +
    theme(legend.position = "right",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 0.5, vjust = 0.5)) +
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(nrow = 1)) #+
    # stat_compare_means(aes(group=type_var), size = 3, paired = TRUE,
    #                     label.y = c(0.75),
    #                     label.x = 1.5,
    #                     label = "p.format", method = "wilcox.test", hide.ns = TRUE)

mann_whitney_test <- function(data) {
  test_result <- wilcox.test(correlation ~ randomization, data = data)
  p_value = test_result$p.value
  # Calculate difference in variance median
  medians <- tapply(data$correlation, data$randomization, median)
  if (length(medians) == 2) { # Ensure there are exactly two groups
    median_diff <- diff(medians)
  } else {
    median_diff <- NA # Assign NA if there aren't two groups
  }
  return(list(p_value = p_value, median_diff = median_diff))
}

results <- input_figureC %>%
  group_by(dataset, type_var) %>%
  nest() %>%
  mutate(test_results = map(data, mann_whitney_test),
         p_value = map_dbl(test_results, 'p_value'),
         median_diff = map_dbl(test_results, 'median_diff')) %>%
  select(-data, -test_results)
results

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig11/panel_D.RData")
input_figureD$dataset <- gsub("_"," ", input_figureD$dataset)

plot_D <- ggplot(data = input_figureD, aes(x = factor(type_var), y = variance, fill = factor(randomization))) +
    geom_violin(position = position_dodge(width = 0.9)) + coord_cartesian(ylim = c(0,5)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    geom_jitter(aes(fill = factor(randomization)), 
                position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                alpha = 0.05, size = 1) +
    facet_wrap(.~dataset) +
    geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
    labs(y = "Variance", x = "", title = "ASE PTC-NMDeff", fill = "") +
    ggplot_theme_bw() +
    theme(legend.position = "right",
      strip.background = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 0.5, vjust = 0.5)) +
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(nrow = 1)) #+
    # stat_compare_means(aes(group=type_var), size = 3,
    #                     label.y = c(3), #paired = TRUE,
    #                     label.x = 1.5,
    #                     label = "p.format", method = "wilcox.test", hide.ns = TRUE)

# For the manuscript
# Function to perform Mann-Whitney U test
mann_whitney_test <- function(data) {
  test_result <- wilcox.test(variance ~ randomization, data = data)
  p_value = as.numeric(test_result$p.value)
  # Calculate difference in variance median
  medians <- tapply(data$variance, data$randomization, median)
  if (length(medians) == 2) { # Ensure there are exactly two groups
    median_diff <- diff(medians)
  } else {
    median_diff <- NA # Assign NA if there aren't two groups
  }
  return(list(p_value = p_value, median_diff = median_diff))
}

results <- input_figureD %>%
  group_by(dataset, type_var) %>%
  nest() %>%
  mutate(test_results = map(data, mann_whitney_test),
         p_value = map_dbl(test_results, 'p_value'),
         median_diff = map_dbl(test_results, 'median_diff')) %>%
  select(-data, -test_results)

####################################################
############ PANEL UP --> PLOTS C-D ################
####################################################

fig_bottom <- ggarrange(plot_D,plot_C, labels = c("C","D"), font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

########################################
############ FINAL SUPP FIG ############
########################################

##########################################
############ FINAL FIG ###################
##########################################

SuppFig_final <- plot_grid(fig_up, fig_middle, fig_bottom, nrow = 3, rel_heights = c(0.35,0.35,0.29))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig11_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 230, height = 270, units = "mm")
