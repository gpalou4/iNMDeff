
source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ FIG 3A ############
################################

# Data
figureA_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures/Schematic_Randomization_test_inter_tissues/Schematic_Randomization_test_inter_tissues.png"
# Read the PNG image using the png package
figureA_png <- readPNG(figureA_path)
# Convert the PNG image to a raster object
raster_grob <- rasterGrob(figureA_png, interpolate=TRUE)
# Create a ggplot object with the raster image
plot_A <- ggplot() + 
  annotation_custom(raster_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void()  # Remove axes and labels

################################
############ FIG 3B ############
################################

# Data
input_figureB <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3B.RData")
input_figureB$database <- ifelse(input_figureB$database == "GTEx","GTex","TCGA")

# Plot
plot_B_top <- input_figureB %>% filter(NMD_method == "ASE") %>% 
        ggplot(aes(x = NMD_geneset, y = median_SD_of_medians_diff_mean, fill = factor(controls))) +
        geom_bar(stat = "identity", linewidth = 0.6, color = "black") + coord_flip() + #coord_cartesian(ylim = c(-5,5)) +
        # Add error bars for confidence intervals
        geom_errorbar(aes(ymin = perc5_SD_of_medians_diff_mean, ymax = perc95_SD_of_medians_diff_mean), 
                        width = 0.2) +  # Adjust width as needed
        guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
        facet_grid(NMD_method ~ database, scale = "free") +
        ylab("") + ggtitle("") + xlab("") +
        theme_bw() + ggplot_theme() +
        # scale_y_break(c(0.38, 0.55)) +
         scale_x_discrete(labels = c("PTC NMD-evading\n0.2" = "NMD-evading PTCs",
                                        "Synonymous 0.2" = "Synonymous",
                                        "PTC NMD-\ntriggering 0.2" = "NMD-triggering PTCs")) +
        #scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2), limits = c(-0.1, 0.2)) +
        scale_fill_brewer(palette = "Pastel1") + scale_color_viridis(discrete=TRUE, option="inferno", direction = -1) +
        theme(plot.title = element_text(hjust = 0.5, size = 35),
                # plot.margin = unit(c(-0.75, 0, -0.5, 0.8), "cm"),
                plot.margin = unit(c(-1.25, 0, -0.25, 0), "cm"),
                axis.text.y = element_text(size = 9),
                legend.position='top') +
        geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
                                position = position_dodge(), size = 9, color = "black", hjust = -0.2, vjust = 0.5)
        # geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
        #                         position = position_dodge(width = 1), size = 10, hjust = -0.1, color = "black", vjust = -0.01)

plot_B_bottom <- input_figureB %>% filter(NMD_method == "ETG") %>% 
        ggplot(aes(x = NMD_geneset, y = median_SD_of_medians_diff_mean, fill = factor(controls))) +
        geom_bar(stat = "identity", linewidth = 0.6, color = "black") + coord_flip() + #coord_cartesian(ylim = c(-5,5)) +
        #ggbreak::scale_y_break(breaks = c(0.35, 0.6), space = 0.5) +
        facet_grid(NMD_method ~ database, scale = "free") +
        geom_errorbar(aes(ymin = perc5_SD_of_medians_diff_mean, ymax = perc95_SD_of_medians_diff_mean), 
                        width = 0.2) +  # Adjust width as needed
        guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
        ylab("Inter-tissue iNMDeff variability deviation") + ggtitle("") + xlab("") +
        theme_bw() + ggplot_theme() + 
        scale_x_discrete(labels = c("RandomGenes\nwithout NMD\nfeatures" = "RandomGenes\nw/o NMD features",
                                "RandomGenes\nwith NMD\nfeatures" = "RandomGenes\nw/ NMD features")) +
        scale_fill_brewer(palette = "Pastel1") + scale_color_viridis(discrete=TRUE, option="inferno", direction = -1) +
        theme(plot.title = element_text(hjust = 0.5, size = 35),
                axis.text.y = element_text(size = 9, lineheight = 0.7),
                plot.margin = unit(c(-1.25, 0, 0.25, 0), "cm"),
                legend.position='none') +
        geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
                                position = position_dodge(), size = 9, color = "black", hjust = -0.2, vjust = 0.5)
        # geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
        #                         position = position_dodge(width = 1), size = 10, hjust = -0.1, color = "black", vjust = -0.01)

plot_B_tmp <- plot_grid(plot_B_top,plot_B_bottom,nrow = 2,  align = "v",
                        rel_heights = c(1, 1), axis = "top")

empty_plot <- plot_grid(NULL, NULL, labels = "")

plot_B <- plot_grid(empty_plot, plot_B_tmp, empty_plot, nrow = 3, labels = c("","",""), label_size = 20, 
                        rel_heights = c(0.04,1,0.04), vjust = 1)

###############################################
############ PANEL UP --> FIG 3A-B ############
###############################################

fig_up <- plot_grid( 
        plot_A,
        plot_B,
        nrow = 1, labels = c("A", "B"), label_size = 20, 
        rel_widths = c(0.55,0.45), rel_heights = c(0.5,0.5))

################################
############ FIG 3C ############
################################

my_colors <- c(
  "Other" = "#F0027F",
  "Pan-GI" = "#FDC086",
  "Pan-kidney" = "#7FC97F",
  "Pan-nervous" = "#BEAED4",
  "Pan-reproductive" = "#FFFF99",
  "Pan-squamous" = "#386CB0"
)

# Data
input_figureC <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3C.RData")

# Stuff
table(input_figureC$pan_organ_system)
input_figureC$pan_organ_system <- gsub("Brain", "Pan-nervous",input_figureC$pan_organ_system)
input_figureC$pan_organ_system <- gsub("Pan-Kidney", "Pan-kidney",input_figureC$pan_organ_system)
input_figureC$pan_organ_system <- gsub("Pan-gyn", "Pan-reproductive",input_figureC$pan_organ_system)
input_figureC$pan_organ_system <- gsub("Pan-Squamous", "Pan-squamous",input_figureC$pan_organ_system)

# Extract colors from the Accent palette
# magma_colors <- viridis(n = length(unique(input_figureC$pan_organ_system)), direction = -1, option = "viridis")
# accent_colors <- brewer.pal(n = length(unique(input_figureC$pan_organ_system)), name = "Accent")[length(unique(input_figureC$pan_organ_system)):1]  # Reverse the colors with [length(...):1]
# # Create a named vector of colors
# my_colors <- setNames(accent_colors, unique(input_figureC$pan_organ_system))
# Correct tissue p-values by FDR
tissues_p_value <- input_figureC[,c("NMD_method","tissues","p_value_above_NMDeff_median","p_value_below_NMDeff_median")]
tissues_p_value <- unique(tissues_p_value)
tissues_p_value <- tissues_p_value %>%
        rowwise() %>%
        mutate(p_value = if_else(p_value_above_NMDeff_median < p_value_below_NMDeff_median, 
                        p_value_above_NMDeff_median, p_value_below_NMDeff_median)) %>%
        group_by(NMD_method) %>%
        mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr"))
# Manual check for manuscript
# df <- tissues_p_value[tissues_p_value$NMD_method == "ASE",]
# table(df$p_value_FDR_adjust < 0.05)

# Add adjusted p-values to the dataframe
input_figureC <- merge(input_figureC,tissues_p_value[,c("NMD_method","tissues","p_value","p_value_FDR_adjust")], by = c("NMD_method","tissues"), all.x = TRUE)
input_figureC <- input_figureC
input_figureC <- input_figureC %>%
  group_by(tissues, NMD_method) %>%
  dplyr::mutate(
    median_score = median(values, na.rm = TRUE),  # Calculate median for each group
    closest = which.min(abs(values - median_score)),  # Find closest row index in each group
    p_value_to_show = ifelse(row_number() == closest, p_value_FDR_adjust, NA)  # Set p_values to NA except for the closest row
  ) %>%
  ungroup() #%>%
#   select(-median_score, -closest)  # Optionally remove the helper columns

# Manual TND and p-value for the Manuscript
df <- input_figureC %>%
        filter(tissues == "KIRC") %>%
        filter(NMD_method == "ASE") %>% as.data.frame()
median(df$values,na.rm=TRUE)
number <- df$p_value_to_show[!is.na(df$p_value_to_show)]
format(number, scientific = TRUE, digits = 3)

plot_C <- input_figureC %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.3,0.3)) +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    labs(y = "Tissue iNMDeff deviation", x = "", fill = "Pan-Organ System", title = "TCGA") +
    ggplot_theme() +
    theme(legend.position = "top") +
    scale_fill_manual(values = my_colors) +
   geom_text(aes(label = ifelse(p_value_to_show < 0.05, "*", "")), 
                position = position_dodge(), size = 6, color = "black", hjust = 0.5, vjust = 0.4)

################################
############ FIG 3D ############
################################

# Data
input_figureD <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3D.RData")
# input_figureD <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3D_subset.RData")
input_figureD$pan_organ_system <- gsub("Brain", "Pan-nervous",input_figureD$pan_organ_system)
input_figureD$pan_organ_system <- gsub("Pan-Kidney", "Pan-kidney",input_figureD$pan_organ_system)
input_figureD$pan_organ_system <- gsub("Pan-gyn", "Pan-reproductive",input_figureD$pan_organ_system)
input_figureD$pan_organ_system <- gsub("Pan-Squamous", "Pan-squamous",input_figureD$pan_organ_system)
# Create a named vector of colors
# accent_colors <- brewer.pal(n = length(unique(input_figureD$pan_organ_system)), name = "Accent")[length(unique(input_figureD$pan_organ_system)):1]  # Reverse the colors with [length(...):1]
# my_colors <- setNames(accent_colors, unique(input_figureD$pan_organ_system))
my_colors <- c(
#   "Other" = "#F0027F",
  "Pan-GI" = "#FDC086",
  "Pan-kidney" = "#7FC97F",
  "Pan-nervous" = "#BEAED4",
  "Pan-reproductive" = "#FFFF99",
  "Pan-squamous" = "#386CB0"
)
# Correct tissue p-values by FDR
tissues_p_value <- input_figureD[,c("NMD_method","tissues","p_value_above_NMDeff_median","p_value_below_NMDeff_median")]
tissues_p_value <- unique(tissues_p_value)
tissues_p_value <- tissues_p_value %>%
        rowwise() %>%
        mutate(p_value = if_else(p_value_above_NMDeff_median < p_value_below_NMDeff_median, 
                        p_value_above_NMDeff_median, p_value_below_NMDeff_median)) %>%
        group_by(NMD_method) %>%
        mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr"))
# Manual check for manuscript (here there are only 33 tissues, look at raw data...)
# df <- tissues_p_value[tissues_p_value$NMD_method == "ETG",]
# table(df$p_value_FDR_adjust < 0.05)
# Add adjusted p-values to the dataframe
input_figureD <- merge(input_figureD,tissues_p_value[,c("NMD_method","tissues","p_value","p_value_FDR_adjust")], by = c("NMD_method","tissues"), all.x = TRUE)
input_figureD <- input_figureD
input_figureD <- input_figureD %>%
  group_by(tissues, NMD_method) %>%
  mutate(
    median_score = median(values, na.rm = TRUE),  # Calculate median for each group
    closest = which.min(abs(values - median_score)),  # Find closest row index in each group
    p_value_to_show = ifelse(row_number() == closest, p_value_FDR_adjust, NA)  # Set p_values to NA except for the closest row
  ) %>%
  ungroup()# %>%
#   select(-median_score, -closest)  # Optionally remove the helper columns


# df1 %>%
#         filter(tissues == "BRNSPC") %>% data.frame()
# df1_pval <- df1 %>%
#         select(NMD_method, tissues,p_value_FDR_adjust) %>%
#         data.frame()
# df2_pval <- df2 %>%
#         select(NMD_method, tissues,p_value_FDR_adjust) %>%
#         data.frame()
# a <- merge(df1_pval,df2_pval, by = c("tissues","NMD_method"),all.x=TRUE)
# a <- a %>%
#         filter(NMD_method=="ETG") %>%
#   mutate(Category = ifelse(grepl("^BRN", tissues), "Brain", "Non-Brain"))
# aggregate(p_value_FDR_adjust.y ~ Category, data = a, median)

# Assuming your data frame is named `data`
# brain_pvalue_counts <- a %>%
#   group_by(NMD_method) %>%
#   filter(Category == "Brain") %>%  # Filter only "Brain" tissues
#   summarise(
#     p_value_FDR_adjust_x_below_0.05 = sum(p_value_FDR_adjust.x < 0.05, na.rm = TRUE),
#     p_value_FDR_adjust_y_below_0.05 = sum(p_value_FDR_adjust.y < 0.05, na.rm = TRUE)
#   )

# # Display the result
# print(brain_pvalue_counts)

# Manual p-value for the Manuscript
df <- input_figureD %>%
        filter(tissues == "LCL") %>%
        filter(NMD_method == "ASE") %>% 
        as.data.frame()
median(df$values,na.rm=TRUE)
number <- df$p_value_to_show[!is.na(df$p_value_to_show)]
format(number, scientific = TRUE, digits = 3)

plot_D <- input_figureD %>%
        ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.6,0.3)) +
        facet_wrap(. ~ NMD_method) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
        labs(y = "Tissue iNMDeff deviation", x = "", fill = "Pan-organ system", title = "GTex") +
        theme_bw() + ggplot_theme() +
        theme(legend.position = "top") +
        #scale_color_brewer(palette = "Accent", direction = -1) +
        scale_fill_manual(values = my_colors) +
        geom_text(aes(label = ifelse(p_value_to_show < 0.05, "*", "")), 
                # position = position_dodge(), size = 9, color = "black", hjust = -0.7, vjust = 0.775)
                position = position_dodge(), size = 6, color = "black", hjust = 0.5, vjust = 0.4)

################################################
############ PANEL MID --> FIG 3C-D ############
################################################

fig_middle <- ggarrange(plot_C, plot_D, labels = c("C","D"), font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, plot_D, width = 250, height = 200, units = "mm") 

################################
############ FIG 3E ############
################################

# Data
input_figureE <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3E.RData")
# Colors
my_colors <- setNames(c("#BEAED4","#aed4d3"), c("GBM","LGG"))

# Plots E 1-3
plot_E_1 <- ggplot(data = input_figureE, aes(x = cancer_type_strat, 
                    y = ETG, fill = factor(cancer_type))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-1.5,1.5)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ETG", x = "", y = "iNMDeff", fill = "") +
        theme_bw() +  ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.title = element_text(size = 12, hjust = 0.5),)
plot_E_2 <- ggplot(data = input_figureE, aes(x = cancer_type_strat, 
                    y = ASE, fill = factor(cancer_type))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-1.5,1.5)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        #facet_grid(. ~ NMD_method , scale = "free") +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ASE", x = "", y = "iNMDeff", fill = "") +
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.title = element_text(size = 12, hjust = 0.5),
                        axis.text.y = element_blank())
plot_E_3 <- ggplot(data = input_figureE, aes(x = cancer_type_strat, 
                    y = neuron, fill = factor(cancer_type))) +
        geom_boxplot() + coord_flip(ylim = c(0,0.5)) + #coord_cartesian(ylim = c(-3,3)) +
        scale_y_continuous(labels=c(0,10,20,30,40,50)) +
        labs(title = "Neuron cell type", x = "", y = "Proportion (%)", fill = "") + 
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "top",
                plot.title = element_text(size = 12, hjust = 0.5),
                axis.text.y = element_blank())

plot_E_1_2_3 <- ggarrange(plot_E_1,plot_E_2,plot_E_3, labels = c("F","",""), font.label = list(size = 20),
                widths = c(0.45,0.28,0.28),
                ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
# Plot E 4
# controls <- c("epithelial.cell","oligodendrocyte","immature.astrocyte","astrocyte","hematopoietic.cell","glial.cell","leukocyte"
# "fibroblast" "macroglial.cell" "precursor.cell" "somatic.cell" "macrophage" "endothelial.cell")
cell_type <- "neuron"
## Split by cancer type
# plot_E_4 <- ggplot(data = input_figureE, aes(x = scale(ETG), y =  eval(parse(text=cell_type))*100, 
#                         color = factor(cancer_type) )) +
#                         # )) +
#         # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white", alpha = 0.35)+
#         geom_smooth(color = "black", method = lm, se = FALSE, alpha = 0.75, size = 0.5) +
#         # geom_smooth(aes(group = factor(cancer_type)),fill = "black", method = lm, se = FALSE) +
#         # geom_point(size = 1, alpha = 0.65, color = "#BEAED4") + 
#         geom_point(size = 1, alpha = 0.65) + 
#         labs(title = "", color = "Cancer type", x = "ETG iNMDeff", 
#                 y = paste0(cell_type," cell type (%)")) +
#         scale_color_manual(values = my_colors) +
#         coord_cartesian(xlim = c(-3,3)) +
#         scale_y_continuous(breaks = c(0,25,50,75,100)) +
#         ggplot_theme() + 
#         ggpubr::stat_cor(aes(color = "black"), color = "black",
#                 method = "pearson", label.y.npc = "top", label.x = -2) +
#         theme(legend.position = "top",
#                 legend.title = element_text(hjust = 0.5),
#                 legend.text = element_text(size = 7),
#                 legend.spacing.y = unit(0.01, "cm"),
#                 plot.margin = unit(c(1, 0, 1, 0), "cm")) +
#         guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
#                                         nrow = 1,
#                                         title.vjust = 2,
#                                         title.position = "top",
#                                         keywidth = unit(0.02, "cm"),
#                                         keyheight = unit(0.02, "cm")) ) 
## All merged
plot_E_4 <- ggplot(data = input_figureE, aes(x = scale(ETG), y =  eval(parse(text=cell_type))*100, 
                        )) +
        geom_smooth(color = "black", method = lm, se = FALSE, alpha = 0.75, size = 0.5) +
        geom_point(size = 1, alpha = 0.65, color = "#BEAED4") + 
        labs(title = "", color = "Cancer type", x = "ETG iNMDeff", 
                y = paste0(cell_type," cell type (%)")) +
        coord_cartesian(xlim = c(-3,3)) +
        scale_y_continuous(breaks = c(0,25,50,75,100)) +
        ggplot_theme() + 
        ggpubr::stat_cor(aes(color = "black"), color = "black",
                method = "pearson", label.y.npc = "top", label.x = -2) +
        theme(legend.position = "top",
                plot.margin = unit(c(1, 0, 1, 0), "cm"))

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test2.png"
# ggsave(final_figure_path, plot_E_4, width = 125, height = 100, units = "mm") 

plot_E <- plot_grid(plot_E_1_2_3, plot_E_4, nrow = 1, ncol = 2,
                        rel_widths = c(0.72,0.26), vjust = 1)

#################################################
############ PANEL BOTTOM 1 --> FIG 3E ############
#################################################

fig_bottom_2 <- plot_E

################################
############ FIG 3F ############
################################

# Data
input_figureF <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3F.RData")
# Colors
my_colors <- setNames(rep(c("#BEAED4"),,14), unique(input_figureF$acronyms))
# Remove NERVET (it has not data for UCD deconvolution)
input_figureF <- input_figureF %>% filter(acronyms != "NERVET")

# Plots E 1-3
plot_F_1 <- ggplot(data = input_figureF, aes(x = acronyms, 
                    y = ETG, fill = factor(acronyms))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-3,1)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ETG", x = "", y = "iNMDeff", fill = "") +
        theme_bw() +  ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "none",
                plot.title = element_text(size = 12, hjust = 0.5),)
plot_F_2 <- ggplot(data = input_figureF, aes(x = acronyms, 
                    y = ASE,fill = factor(acronyms))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_flip(ylim = c(-2,1)) +
        geom_boxplot(width =  0.3, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(title = "ASE", x = "", y = "iNMDeff", fill = "") +
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "none",
                plot.title = element_text(size = 12, hjust = 0.5),
                        axis.text.y = element_blank())
plot_F_3 <- ggplot(data = input_figureF, aes(x = acronyms, 
                    y = neuron, fill = factor(acronyms))) +
        geom_boxplot() + coord_flip(ylim = c(0,1)) + #coord_cartesian(ylim = c(-3,3)) +
        scale_y_continuous(labels=c(0,25,50,75,100)) +
        labs(title = "Neuron cell type", x = "", y = "Proportion (%)", fill = "") + 
        theme_bw() + ggplot_theme() +
        scale_fill_manual(values = my_colors) +
        theme(legend.position = "none",
                plot.title = element_text(size = 12, hjust = 0.5),
                axis.text.y = element_blank())

plot_F_1_2_3 <- ggarrange(plot_F_1,plot_F_2,plot_F_3, labels = c("E","",""), font.label = list(size = 20),
                widths = c(0.42,0.28,0.28), ncol=3, nrow=1, common.legend = TRUE, legend="none")

# Plot F 4
cell_type <- "neuron"
## Splitting by brain subregion
plot_F_4 <- ggplot(data = input_figureF, aes(x = scale(ETG), y =  eval(parse(text=cell_type))*100, 
                        color = factor(acronyms) )) +
                        # )) +
        # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white", alpha = 0.35)+
        geom_smooth(color = "black", method = lm, se = FALSE, alpha = 0.75, size = 0.5) +
        # geom_smooth(aes(group = factor(cancer_type)),fill = "black", method = lm, se = FALSE) +
        # geom_point(size = 1, alpha = 0.65, color = "#BEAED4") + 
        geom_point(size = 1, alpha = 0.65) + 
        labs(title = "", color = "Brain subregion", x = "ETG iNMDeff", 
                y = paste0(cell_type," cell type (%)")) +
        # scale_color_manual(values = my_colors) +
        coord_cartesian(xlim = c(-3,3)) +
        scale_y_continuous(breaks = c(0,25,50,75,100)) +
        ggplot_theme() + 
        ggpubr::stat_cor(aes(color = "black"), color = "black",
                method = "pearson", label.y.npc = "top", label.x = -2) +
        theme(legend.position = "top",
                legend.title = element_text(hjust = 0.5),
                legend.text = element_text(size = 7),
                legend.spacing.y = unit(0.01, "cm"),
                plot.margin = unit(c(0.5, 0, 0.5, 0), "cm")) +
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                        nrow = 5,
                                        title.vjust = 2,
                                        title.position = "top",
                                        keywidth = unit(0.02, "cm"),
                                        keyheight = unit(0.02, "cm")) ) 
## All merged
plot_F_4 <- ggplot(data = input_figureF, aes(x = scale(ETG), y =  eval(parse(text=cell_type))*100, 
                        )) +
        geom_smooth(color = "black", method = lm, se = FALSE, alpha = 0.75, size = 0.5) +
        geom_point(size = 1, alpha = 0.65, color = "#BEAED4") + 
        labs(title = "", color = "Brain subregion", x = "ETG iNMDeff", 
                y = paste0(cell_type," cell type (%)")) +
        coord_cartesian(xlim = c(-3,3)) +
        scale_y_continuous(breaks = c(0,25,50,75,100)) +
        ggplot_theme() + 
        ggpubr::stat_cor(aes(color = "black"), color = "black",
                method = "pearson", label.y.npc = "top", label.x = -2) +
        theme(legend.position = "top",
                plot.margin = unit(c(1, 0, 1, 0), "cm"))

plot_F <- plot_grid(plot_F_1_2_3, plot_F_4, nrow = 1, ncol = 2,
                        rel_widths = c(0.75,0.25), vjust = 1)

###################################################
############ PANEL BOTTOM 2 --> FIG 3E ############
###################################################

fig_bottom_1 <- plot_F

#############################################
############ FINAL FIG --> FIG 3 ############
#############################################

# Fig3_final <- plot_grid(fig_up, fig_middle, fig_bottom, nrow = 3, rel_heights = c(0.4,0.6,0.4))
Fig3_final <- plot_grid(fig_up, fig_middle, fig_bottom_1, fig_bottom_2, 
                 nrow = 4, rel_heights = c(0.4,0.6,0.4,0.4)) + ggplot_theme()
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig3_complete.png"
ggsave(final_figure_path, Fig3_final, width = 250, height = 400, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig3_complete.pdf"
ggsave(final_figure_path, Fig3_final, width = 250, height = 400, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig3_complete.tiff"
ggsave(final_figure_path, Fig3_final, width = 250, height = 400, units = "mm", device = "tiff")
