library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(corrplot)
library(ggcorrplot)
library(tidyr)
library(dplyr)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

# Data
input_figureA <- read.table(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/panel_A.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
input_figureA <- input_figureA %>%
                        filter(NMD_gene_excluded != "all")

# Calculate mean and CI for each group (Type and CL)
summary_df <- input_figureA %>%
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

plot_A <- ggplot(input_figureA, aes(x = CL, y = nb_coeff, fill = Type)) +
    geom_bar(
    data = input_figureA %>%
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
        label.y = max(input_figureA$nb_coeff) + 0.5  # Adjust for p-value placement
    ) +
    labs(
        title = "",
        fill = "",
        x = "Cell lines",
        y = "ETG cNMDeff"
    ) +
    ggplot_theme() +
    scale_fill_manual(values = c("WT" = "#56B4E9", "UPF1 KD" = "#E69F00")) +
    theme(
        legend.title = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "top"
    )

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/panel_B_C.RData")

input_figureB <- input_figureB[,c("score","TCGA_brain_enriched","GTEx_brain_enriched","final_consensus","cell_type")]

# Pivot the first two columns into a long format
input_figureB <- input_figureB %>%
  pivot_longer(
    cols = c(TCGA_brain_enriched, GTEx_brain_enriched),  # Columns to pivot
    names_to = "Cohort",                                # New column for names
    values_to = "Brain_enrichment"                      # New column for values
  ) %>%
  mutate(Cohort = ifelse(Cohort == "TCGA_brain_enriched", "TCGA", "GTex")) %>%
  mutate(Brain_enrichment = ifelse(Brain_enrichment == "Brain-enriched", "Brain\nenriched", "Non-brain\nenriched"))

# Plot dREG score density for TCGA (using TCGA_brain_enriched for fill)
plot_B <- ggplot(input_figureB, 
               aes(y = score, x = Brain_enrichment, fill = factor(Brain_enrichment),
                   color = final_consensus)) +
  # Boxplot with transparency
  geom_boxplot(alpha = 0.7) +
  # Jittered points
  geom_jitter(aes(color = final_consensus), 
              position = position_jitterdodge(), 
              size = 2, alpha = 0.5) +
  # Modern color palette for fill and color
  scale_fill_brewer(palette = "Set2") +  # Updated softer palette
  scale_color_brewer(palette = "Dark2") +  # Updated softer palette for consensus
  # Facet customization
#   facet_wrap(Cohort ~ cell_type, 
#              scales = "free") +  # Move panel titles below
#   facet_wrap(Cohort ~ cell_type, 
#              scales = "free", 
#              strip.position = "bottom") +  # Move panel titles below
  theme(strip.background = element_blank(),  # Remove panel background
        strip.text = element_text(size = 10, face = "bold")) +  # Keep titles
  # Add horizontal line and p-values
  geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3, label.y = 2 ) +
  # Labels and axis adjustments
  labs(title = "", fill = "", color = "Transcript",
       x = "", y = expression(log[10]("PRO-Seq score"))) +
  ggplot_theme() +
  guides(fill = "none") +
  theme(
    # legend.title = element_blank(), 
    plot.title = element_text(hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5), # Rotate X-axis labels
    strip.background = element_blank(),  # Remove strip background
    strip.placement = "outside"          # Keep facet titles as outer labels
  ) +
  # Reorder facets
  facet_grid(Cohort ~ factor(cell_type, levels = c("SH_SY5Y_1", "SH_SY5Y_2", "U2OS", "A549")))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
ggsave(final_figure_path, plot_B, width = 100, height = 150, units = "mm") 

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

##################################################
############ PANEL BOTTOM --> PLOTS B ############
##################################################

fig_bottom <- plot_grid(plot_B, nrow = 1, ncol = 1, labels = c("B"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

##################################################
############ FINAL FIG --> SUPP FIG 2 ############
##################################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, rel_heights = c(0.5,0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig3_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 175, height = 225, units = "mm")
