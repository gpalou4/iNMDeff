# conda activate R_figures
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(stringr)
library(dplyr)
library(png)
library(ggpubr)
library(ggbreak)
library(scales)
library(tibble)
library(corrplot)
library(ggcorrplot)
library(MASS)
library(ggh4x)
library(survminer)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

# Function to save ggsurv plots with "ggsave"
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)} 

################################
############ PLOT A ############
################################

input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig28/panel_A.RData")

# Define specific y-axis ranges for each cell type
cell_type_ranges <- list(
  "CD8T_cytotoxic" = c(0, 0.025),
  "CD8T_t_cell" = c(0, 20),
  "CD8T_positive_t_cell" = c(0, 0.02),
  "CD8T_memory" = c(0, 0.01),
  "CD8T_cytokine" = c(0, 20)
)

# Create a named vector for renaming cell types
cell_type_labels <- c(
  "CD8T_cytotoxic" = "CD8+T Cytotoxic",
  "CD8T_t_cell" = "CD8+T",
  "CD8T_positive_t_cell" = "CD8+T activated",
  "CD8T_memory" = "CD8+T memory",
  "CD8T_cytokine" = "CD8+T cytokine"
)

# Initialize an empty list to store the plots
plots <- list()
is_first_plot <- TRUE

# Loop over each immune cell type to create individual plots
for (cell_type in names(cell_type_ranges)) {
  # Subset the data for the current immune cell type
  subset_data <- input_figureA %>% 
    filter(immune_cell_type == cell_type)
  
  # Get the y-axis range for the current cell type
  y_range <- cell_type_ranges[[cell_type]]

  legend_position <- if (is_first_plot) "top" else "none"

  # Create the plot
  plot <- ggplot(subset_data, aes(x = factor(NMD_type), y = value * 100, fill = NMD_type)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.5) +
    geom_point(
      aes(color = NMD_type), 
      position = position_dodge(0.8), 
      size = 2, 
      alpha = 0.7
    ) +
    facet_grid(immune_cell_type~primaryTumorLocation, scale = "free_y",
            labeller = labeller(immune_cell_type = cell_type_labels)) +
    labs(
      title = "",
      x = "",
      y = "cell type (%)",
      fill = "ETG iNMDeff"
    ) +
    ggplot_theme() +
    coord_cartesian(ylim = y_range) +
    stat_compare_means(
      # aes(group = interaction(TIB_type, NMD_type)),
      label.y = y_range[2]-(y_range[2]*0.1),
      label.x = 1.5,
      size = 3.5,
      method = "wilcox.test",
      # method.args = list(alternative = "less"),
      label = "p.signif"
      # label = function(x) sprintf("%.3g", x$p)
    ) + guides(color = "none") + 
    theme(strip.text = element_text(size = 7),
            legend.position = legend_position,
            axis.text.x = element_blank(),
            axis.title.y =  element_text(size = 10)
            )
  # After the first plot, set is_first_plot to FALSE
  is_first_plot <- FALSE
  
  # Save the plot in the list
  plots[[cell_type]] <- plot
}

################################################
############ PANEL UP --> PLOTS A ############
################################################

fig_up <- plot_grid(plots[[1]], plots[[2]],plots[[3]],plots[[4]],plots[[5]], nrow = 5,
                       labels = c(""), label_size = 20, align = "v", rel_heights = c(0.23,0.2,0.2,0.2,0.2))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig28_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 130, height = 260, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig28_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 130, height = 260, units = "mm")
