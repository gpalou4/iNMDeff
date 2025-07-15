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
library(ggforce)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

# Data
input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig20/panel_A.RData")

# Colors
input_figureA$label_color <- "#8A8A8A"
input_figureA[input_figureA$gene_symbol %in% "SMG5","label_color"] <- "#2015D0"
tmp <- subset(input_figureA, gene_symbol %in% top_hits)
label_color <- tmp$label_color
input_figureA$controls_vs_rest <- NA
input_figureA[which(input_figureA$selected_bands == "other"),"controls_vs_rest"] <- "other"
input_figureA[which(input_figureA$selected_bands != "other"),"controls_vs_rest"] <- "1q21-23.1"
control_bands <- c("q24.1","q24.2","q24.3","q25.1","q25.2","q25.3","q31.1","q31.2","q31.3","q32.1","q32.2","q32.3","q41","q42.11","q42.12","q42.13","q42.2",
                "q42.3","q43","q44")
selected_bands_char <-  c("q21.1","q21.2","q21.3","q22","q23.1","q23.2","q23.3")
selected_bands_2_char <- c("q21.3")
input_figureA <- input_figureA %>%
        mutate(selected_bands_2 = if_else(band %in% control_bands, "1q-end",band)) %>%
        mutate(selected_bands_2 = if_else(selected_bands_2 %in% selected_bands_char, "1q21.1-23.3",selected_bands_2)) #%>%
table(input_figureA$selected_bands_2)
# Extract colors from the Brewer palette
brewer_colors <- scales::brewer_pal(palette = "Set2")(length(unique(input_figureA$selected_bands_2)))
# Modify the color for the "other" group
names(brewer_colors) <- unique(input_figureA$selected_bands_2)
brewer_colors[names(brewer_colors) == "other"] <- "#808080"
# Manual change
brewer_colors[names(brewer_colors) == "1q-end"] <- "#181A1D"
brewer_colors[names(brewer_colors) == "1q21.1-23.3"] <- "#DCA620"


# Assign categories based on gene membership
input_figureC <- input_figureC %>%
  mutate(
    category = case_when(
      gene %in% candidate_genes_nonNMD ~ "Candidates",
      gene %in% candidate_genes_NMD ~ "Candidates NMD",
      gene %in% NMD_factors ~ "NMD-related",
      TRUE ~ "Chr 1q - other"
    ),
    label_colors = case_when(
      category == "Candidates" ~ "orange",
      category == "Candidates NMD" ~ "red",
      category == "NMD-related" ~ "blue",
      category == "Chr 1q - other" ~ "black",
      TRUE ~ "gray"  # Assign gray for any unknown genes
    )
  )

plot_A1 <- ggplot(data = input_figureA, aes(x = CNA_RNA_corr, y = ASE_iNMDeff_RNA_corr, 
                        color = factor(selected_bands_2)), fill = factor(selected_bands_2)) +
        annotate("rect", xmin = 0.25, xmax = Inf, ymin = -Inf, ymax = -0.15, fill= "#FCF1BC") + 
        annotate("text", x = 0.45, y = -0.3, label = "Candidate genes") + 
        stat_density_2d(data = input_figureA %>% filter(selected_bands_2 == "1q-end"),
                        aes(fill = stat(level)), geom = "polygon", alpha = 0.5) +
        geom_smooth(aes(group = factor(selected_bands_2)),fill = "black", method = lm, se = FALSE) +
        geom_point(data = input_figureA %>% filter(selected_bands_2 != "1q-end"), size = 2, alpha = 0.35) + 
        geom_point(data = input_figureA %>% filter(selected_bands_2 == "1q-end"), shape = 3, size = 2, alpha = 0.5) + 
        geom_vline(xintercept= 0.25, size = 0.5, color = "black", linetype = "dashed") +
        geom_hline(yintercept= -0.15, size = 0.5, color = "black", linetype = "dashed") +
        labs(title = "Chr 1q genes", color = "Chr1 bands", x = "Gene exp vs CNA correlation", y = "Gene exp vs ASE iNMDeff correlation") +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = c(0,0.1,0.2,0.25,0.3,0.4,0.5)) +
        scale_y_continuous(breaks = c(0.1,0,-0.1,-0.15,-0.2,-0.3)) +
        # geom_text_repel(data = subset(input_figureA, gene_symbol %in% top_hits & !gene_symbol %in% c("SMG5")),
        #                 max.overlaps = nrow(input_figureA), #color = label_color,
        #                 aes(label = gene_symbol), size = 3,      arrow = arrow(length = unit(0.015, 'npc')),
        #                 point.padding = 0.05, nudge_x = .04, nudge_y = -.03) +
        geom_text_repel(data = subset(input_figureA, gene_symbol %in% c("GON4L","VPS72","PRPF3","PMF1","MSTO1")),
                        max.overlaps = nrow(input_figureA), #color = label_color,
                        aes(label = gene_symbol), size = 3,      arrow = arrow(length = unit(0.015, 'npc')),
                        point.padding = 0.05, nudge_x = .04, nudge_y = -.03) +
        geom_text_repel(data = subset(input_figureA, gene_symbol %in% c("SMG5","INTS3","SF3B4","RBM8A")),
                        max.overlaps = nrow(input_figureA), color = "#ff0000",
                        aes(label = gene_symbol), size = 3, arrow = arrow(length = unit(0.01, 'npc')),
                        point.padding = 0.05, nudge_x = .02, nudge_y = -.02) +
        scale_color_manual(values = brewer_colors) + # Use the modified colors
        ggplot_theme()+ 
        ggpubr::stat_cor(method = "pearson", label.y.npc = "bottom") +
        theme(legend.position = "bottom")

input_figureA2 <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig20/panel_A_2.RData")

input_figureA2_stack <- stack(input_figureA2[,c("SMG5_RNAseq","PMF1_RNAseq")])
input_figureA2_stack$ASE_PTC_NMD_triggering_0.2 <- rep(input_figureA2$ASE_PTC_NMD_triggering_0.2,2)
colnames(input_figureA2_stack) <- c("RNAseq","gene","ASE")
input_figureA2_stack$gene <- gsub("_RNAseq","",input_figureA2_stack$gene)

plot_A2 <- ggplot(data = input_figureA2_stack, aes(x = RNAseq, y = ASE)) +
                facet_wrap(.~gene) +
                geom_point(size = 1, alpha = 0.2, color = "#3b528b",  shape = 21) +
                stat_density_2d(aes(fill = ..level..), geom = "polygon") + guides(fill = FALSE) +
                geom_smooth(color = "#3b528b", method = lm, se = TRUE) +
                coord_cartesian(ylim = c(-5,5), xlim = c(0,150)) +
                ggplot_theme_bw()+ 
                scale_x_continuous(breaks = c(50,150)) +
                labs(title = "", x = "Gene expression (TPM)", y = "ASE iNMDeff") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom")

# Plot A3
input_figureA2$gene_1 <- "PMF1"
plot_A3 <- ggplot(data = input_figureA2, aes(x = PMF1_RNAseq, y = PMF1_CNA)) +
                facet_wrap(.~gene_1) +
                coord_cartesian(ylim = c(-1,4), xlim = c(0,200)) +
                ggplot_theme_bw()+ 
                geom_bin2d(bins = 50, alpha = 0.7) + guides(fill = FALSE) +   
                geom_smooth(method = lm, se = TRUE, color = "black") +
                scale_fill_gradient(low = "#B03838", high = "#FF0000") +
                labs(title = "", x = "", y = "CNA") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom",
                        plot.margin = unit(c(0, 0, 0, 0), "cm"),
                        axis.title.y = element_text(margin = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt")))
# Plot A4
input_figureA2$gene_2 <- "SMG5"
plot_A4 <- ggplot(data = input_figureA2, aes(x = SMG5_RNAseq, y = SMG5_CNA)) +
                facet_wrap(.~gene_2) +
                coord_cartesian(ylim = c(-1,4), xlim = c(0,200)) +
                ggplot_theme_bw()+ 
                geom_bin2d(bins = 50,  alpha = 0.7) + guides(fill = FALSE) +   
                geom_smooth(method = lm, se = TRUE, color = "black") +
                scale_fill_gradient(low = "#B03838", high = "#FF0000") +
                labs(title = "", x = "", y = "") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom",
                        axis.title.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

plot_A2_3_4 <- plot_grid(plot_A2, plot_A3, plot_A4, nrow = 1, ncol = 3, rel_widths = c(0.5,0.27,0.23), align = "hv")
plot_A_tmp <- plot_grid(plot_A1, plot_A2_3_4, nrow = 2, ncol = 1, rel_heights = c(0.65,0.35))
plot_A <- ggdraw() + 
  draw_plot(plot_A_tmp) +
  draw_label("Gene expression (TPM)", x = 0.76, y = 0.015, hjust = 0.5, size = 12) +
        theme(plot.margin = unit(c(0, 4, 0, 4), "cm"))

################################
############ PLOT B ############
################################

# Data
input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig20/panel_B.RData")

# NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
candidate_genes <- rownames(input_figureB[input_figureB$label_colors == "orange",])
candidate_genes_nonNMD <- candidate_genes[!candidate_genes %in% NMD_genes$gene_symbol]
candidate_genes_NMD <- c("SMG5","RBM8A",candidate_genes[candidate_genes %in% NMD_genes$gene_symbol])
control_genes <- rownames(input_figureB[input_figureB$label_colors == "blue",])
# Remove some genes
NMD_related_genes_to_remove <- NMD_genes[NMD_genes$NMD_type %in% c("NMD_related","NMD_ER"),"gene_symbol"]
NMD_factors <- NMD_genes[!NMD_genes$NMD_type %in% c("NMD_related","NMD_ER"),"gene_symbol"]
NMD_factors <- NMD_factors[!NMD_factors %in% "SMG5"]
# Filter selected genes
genes_to_keep <- c(NMD_factors,candidate_genes_nonNMD,candidate_genes_NMD,control_genes)
input_figureB <- input_figureB[rownames(input_figureB) %in% genes_to_keep, colnames(input_figureB) %in% c("label_colors",genes_to_keep)]
# Colors
input_figureB$gene <- rownames(input_figureB)
input_figureB$category <- NA

# Assign categories based on gene membership
input_figureB <- input_figureB %>%
  mutate(
    category = case_when(
      gene %in% candidate_genes_nonNMD ~ "Candidates",
      gene %in% candidate_genes_NMD ~ "Candidates NMD",
      gene %in% NMD_factors ~ "NMD-related",
      gene %in% control_genes ~ "Controls",
      TRUE ~ "Unknown"
    ),
    label_colors = case_when(
      category == "Candidates" ~ "orange",
      category == "Candidates NMD" ~ "red",
      category == "NMD-related" ~ "blue",
      category == "Controls" ~ "black",
      TRUE ~ "gray"  # Assign gray for any unknown genes
    )
  )
table(input_figureB$label_colors,input_figureB$category)

df <- unique(merge(unique(input_figureB[,c("gene","category")]),NMD_genes[,c("gene_symbol","NMD_type")], 
    by.x = "gene", by.y = "gene_symbol", all.x = TRUE))
table(df$category, df$NMD_type)
  
plot_B <- ggcorrplot(as.matrix(input_figureB[,!colnames(input_figureB) %in% c("label_colors","gene","category","color")]), method = "square", 
                        hc.order = TRUE, type = "full", lab = FALSE,
                        title = "CRISPR KO co-dependency scores ",
                        ggtheme = theme_bw(),
                        tl.cex = 8,
                        tl.srt= 55,
                        tl.col = input_figureB$label_colors,
                        show.diag = TRUE)

# Extract the data from the ggplot object
# plot_data <- levels(plot_B$data$Var1)

# Match the order of label_colors to the order of labels in the plot
ordered_colors <- input_figureB$label_colors[match(levels(plot_B$data$Var1), rownames(input_figureB))]

plot_B$data$label_colors <- ordered_colors 
plot_B$data$label_colors <- ifelse(plot_B$data$label_colors == "blue","NMD-related",plot_B$data$label_colors)
plot_B$data$label_colors <- ifelse(plot_B$data$label_colors == "black","Controls",plot_B$data$label_colors)
plot_B$data$label_colors <- ifelse(plot_B$data$label_colors == "orange","Candidates",plot_B$data$label_colors)
plot_B$data$label_colors <- ifelse(plot_B$data$label_colors == "red","Candidates NMD",plot_B$data$label_colors)

vector_colors <- c("#0000FF","#FFA500","#000000","#FF0000")
names(vector_colors) <- c("NMD-related","Candidates","Controls","Candidates NMD")

plot_B <- plot_B +
        geom_point(aes(color = factor(label_colors)), alpha = 0) +
        scale_color_manual(limits = names(vector_colors), values = vector_colors)+
        scale_fill_gradient2(limits = c(0,0.1), low = "#2A0BEE", mid = "#EDFEFE", high = "#EE0B0B") +
        labs(fill = "Scores", x = "", y = "", color = "Genes") + ggplot_theme() +
        guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
        theme(plot.title = element_text(hjust = 0.5),
                plot.background = element_rect("white"),
                # plot.margin = unit(c(0.1, 0.95, -0.1, 1.35), "cm"),
                # plot.margin = unit(c(-0.5, 0.5, -0.5, 0.5), "cm"),
                legend.position = "right",
                axis.text.x = element_text(color = ordered_colors, angle = 90, hjust = 1, vjust = 0.5, size = 4),
                axis.text.y = element_text(color = ordered_colors, size = 4))

################################
############ PLOT C ############
################################

# Data
input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig20/panel_C.RData")
input_figureC <- input_figureC[input_figureC$gene != "SMG7",]
input_figureC$gene_type <- ifelse(input_figureC$gene_type == "chr1q","other","NMD related")

# Assign categories based on gene membership
input_figureC <- input_figureC %>%
  mutate(
    category = case_when(
      gene %in% candidate_genes_nonNMD ~ "Candidates",
      gene %in% candidate_genes_NMD ~ "Candidates NMD",
#       gene %in% NMD_factors ~ "NMD-related",
      TRUE ~ "Other"
    ),
    label_colors = case_when(
      category == "Candidates" ~ "orange",
      category == "Candidates NMD" ~ "red",
#       category == "NMD-related" ~ "blue",
      category == "Other" ~ "#7d7d7d",
      TRUE ~ "gray"  # Assign gray for any unknown genes
    )
  )
table(input_figureC$label_colors,input_figureC$category)
vector_colors <- c("#FFA500","#7d7d7d","#FF0000")
names(vector_colors) <- c("Candidates","Other","Candidates NMD")

plot_C <- ggplot(data = input_figureC, aes(x = as.numeric(as.character(genome_location_chr)), 
                        y = -log10(p_value_FDR_adjust), fill = as.numeric(gene_type),
                        color = factor(category))) +
        geom_point(size = 3) +
        geom_text_repel(data = subset(input_figureC, p_value_FDR_adjust < 0.1 | gene %in% c("SF3B1","CWC22","SF3B4","FARSB")),
                aes(label = gene), size = 3, arrow = arrow(length = unit(0.005, 'npc'), type = "closed"),
                        point.padding = 0.005, nudge_x = .01, nudge_y = .003, max.overlaps = 200) +
        scale_color_manual(limits = names(vector_colors), values = vector_colors)+
        ggplot_theme()+
        guides(color = guide_legend(nrow = 2, override.aes = list(size = 4))) +  # Two rows in legend
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                legend.text = element_text(colour = "black", size = 10),
                legend.position = "top") +
        labs(title = "", x = "Chr 1q21.1-23.1", y = expression("-log"[10] * "(FDR)"), colour = "Genes") +
        geom_hline(yintercept= -log10(0.05), size = 0.75, color = "black", linetype = "dotted") +
        geom_hline(yintercept= -log10(0.1), size = 0.75, color = "red", linetype = "dotted") #+
        # geom_hline(yintercept= -log10(0.25), size = 0.75, color = "red", linetype = "dotted")

###################################################
############ PANEL UP --> FIG A & B ###############
###################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"))

###################################################
############ PANEL BOTTOM --> FIG C ###############
###################################################

# plot_B_E <- plot_grid(plot_B, plot_E, nrow = 2, ncol = 1, labels = c("B","D"))

fig_bottom <- plot_grid(plot_B, plot_C, nrow = 1, ncol = 2, labels = c("B","C"))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.55,0.45)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig20_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 225, height = 280, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig20_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 225, height = 280, units = "mm")
