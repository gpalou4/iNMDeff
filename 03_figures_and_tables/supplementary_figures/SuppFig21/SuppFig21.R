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

    strip.text = element_text(colour = "black", size = 12)
  )
}

################################
############ PLOT A ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/panel_A.RData")
# input_figureB <- input_figureB[-which(input_figureB$gene_name %in% c("SNORA44","RN7SL371P","RNA5SP162")),]

# Plot
plot_A <- ggplot(data = input_figureB, aes(x = genome_location, y = values, group = ind, color = ind)) +
        geom_point(size = 1) + ggtitle(paste0("")) +
        geom_line() +  
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.5,2)) +
        geom_hline(yintercept=c(0,0.5)) +
        labs(color = "CNA-PC52 groups", x = "Chr2 Genome location", y = "GISTIC CNA average state") +
        ggplot_theme() +
        scale_color_brewer(palette="Set2") +
        theme(legend.position = "top") +
        scale_x_discrete(guide=guide_axis(n.dodge = 3)) +
        guides(color = guide_legend(override.aes = list(size = 6), nrow = 1))

################################
############ PLOT B ############
################################

# Plot D1
input_figureB_1 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/panel_B.RData")
combinations_B_1 <- combn(names(table(input_figureB_1$bins)), 2, simplify = FALSE)

plot_B_1 <- input_figureB_1 %>% 
        filter(NMD_method == "Endogenous") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "", x = "CNA-PC52 groups", y = "ETG iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                        plot.margin = unit(c(1, 0, 1, 0), "cm")) +
                stat_compare_means(comparisons = combinations_B_1, size = 4,
                                label.y = c(-0.25,0.25,0.75),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

# Plot E2
input_figureB_2 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/panel_B.RData")
combinations_B_2 <- combn(names(table(input_figureB_2$bins)), 2, simplify = FALSE)

plot_B_2 <- input_figureB_2 %>% 
        filter(NMD_method == "ASE") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "", x = "CNA-PC52 groups", y = "ASE iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                        plot.margin = unit(c(1, 0, 1, 0), "cm")) +
                stat_compare_means(comparisons = combinations_B_2, size = 4,
                                label.y = c(0.5,1,1.35),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

plot_B <- plot_grid(plot_B_1, plot_B_2, nrow = 1, ncol = 2)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/panel_C.RData")

# Plot
plot_C <- ggplot(input_figureC, aes(y = round(samples_percentage,2), x = cancer_type, fill = CNA_PC_bins)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = paste0(""), x = "", y = "% of individuals", fill = paste0("Group")) + 
    scale_fill_brewer(labels = c("Low", "Mid", "High"), palette = "Set2", direction = -1) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5, size = 9),
            plot.margin = unit(c(1.5, 0.25, 1.5, 0.25), "cm"),
            legend.position = "top")

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/panel_D_1.RData")
# Merge ASE+ETG
# Pivot the specified columns
input_figureD <- input_figureD %>%
  pivot_longer(
    cols = c(NMDeff_endogenous_1, NMDeff_ASE_1),  # Columns to pivot
    names_to = "iNMDeff_type",                          # New column to store the names
    values_to = "iNMDeff"                           # New column to store the values
  ) %>% data.frame()

input_figureD$iNMDeff_type <- ifelse(input_figureD$iNMDeff_type == "NMDeff_endogenous_1","ETG","ASE")

head(input_figureD)
dim(input_figureD)

# Data
input_figureD <- input_figureD[order(input_figureD$genome_location),]
CNV_color <- brewer.pal(n = 2, name = "Paired")
CNV_color <- CNV_color[-3]
# Outlier
# input_figureA <- input_figureA[-which(input_figureA$gene_symbols == "SNORA44"),]

# Create a fake color only to output the legend
input_figureD$CNA_type <- "CNA_freq_1"
input_figureD$CNA_type[1:20] <- "CNA_freq_2"

plot_D1 <- ggplot(data = input_figureD, aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1, fill = factor(CNA_type) )) +
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1), stat = "identity", color = CNV_color[1]) +
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_2), stat = "identity", color = CNV_color[2]) +
                labs(title = "", x = "", y = "CNA freq") +
                geom_hline(yintercept= 0, size = 0.5, color = "black") +
                ggplot_theme() +
                theme(axis.text.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        plot.margin = unit(c(0, 0.5, -0.5, 0.25), "cm"),
                        legend.position = "top") +
                scale_fill_manual(values = c("CNA_freq_1" = CNV_color[1], "CNA_freq_2" = CNV_color[2]),
                                labels = c("CNA_freq_1" = "> 0", "CNA_freq_2" = ">= 1.5"),
                                  name = "GISTIC score")

# Plot A-2 (bottom)
# Colors and labels
coAmp_color <- brewer.pal(n = 8, name = "OrRd")
candidate_genes_NMD <- c("CWC22","SF3B1","NOP58","FARSB")
NMD_related_genes <- c("")
genes_label <- c("CWC22","SF3B1","NOP58","FARSB","TMEM237","MFF","SSB","KCTD18","TYW5","OSGEPL1","CIR1","PMS1","BCS1L","ORMDL1","CTDSP1","RNF25","PRKRA","METTL21A")
genes_label <- c("WDR75","PRKRA","SPC25","BARD1","CFLAR","CTDSP1","WDR12","RPL37A","WDR75")
genes_label <- c(candidate_genes_NMD,genes_label,NMD_related_genes)
input_figureD$label_color <- "#000000"
input_figureD[input_figureD$gene_symbols %in% candidate_genes_NMD,"label_color"] <- "#ff0000"
input_figureD[input_figureD$gene_symbols %in% NMD_related_genes,"label_color"] <- "#2015D0"
input_figureD[input_figureD$gene_symbols %in% genes_label & !input_figureD$gene_symbols %in% 
                        c(candidate_genes_NMD,NMD_related_genes),"label_color"] <- "#DCA620"
tmp <- subset(input_figureD, (gene_symbols %in% genes_label) & (iNMDeff_type %in% c("ASE")))
label_color <- tmp$label_color
input_figureD$genome_location <- factor(input_figureD$genome_location, levels = unique(input_figureD$genome_location))

plot_D2 <- ggplot(data = input_figureD, aes(x = as.numeric(as.character(genome_location)), y = iNMDeff, fill = iNMDeff_type,
                        colour = factor(label_color))) +
        geom_point(size = 1.25, shape = 21,  stroke = 0.4, alpha = 0.8) +
        labs(title = "", x = "Chr 2 genome location", y = "median iNMDeff", fill = "") +
        geom_hline(yintercept= 0, size = 0.5, color = "black") +
        # 2q31.1-36.3 bands
        geom_vline(xintercept = c(168834132,176664676), size = 0.5, color = "black", linetype = "dashed") + # q31.1
        geom_vline(xintercept = c(177138699,179101692), size = 0.5, color = "black", linetype = "dashed") + # q31.2
        geom_vline(xintercept = c(179441982,182048822), size = 0.5, color = "black", linetype = "dashed") + # q31.3
        geom_vline(xintercept = c(182140036,188297492), size = 0.5, color = "black", linetype = "dashed") + # q32.1
        geom_vline(xintercept = c(188734155,190964358), size = 0.5, color = "black", linetype = "dashed") + # q32.2
        geom_vline(xintercept = c(191029576,196289409), size = 0.5, color = "black", linetype = "dashed") + # q32.3
        geom_vline(xintercept = c(196639554,202376936), size = 0.5, color = "black", linetype = "dashed") + # q33.1
        geom_vline(xintercept = c(202635188,203936748), size = 0.5, color = "black", linetype = "dashed") + # q33.2
        geom_vline(xintercept = c(204545793,208165343), size = 0.5, color = "black", linetype = "dashed") + # q33.3
        geom_vline(xintercept = c(208236227,214411065), size = 0.5, color = "black", linetype = "dashed") + # q34
        geom_vline(xintercept = c(214725646,219906502), size = 0.5, color = "black", linetype = "dashed") + # q35
        geom_vline(xintercept = c(221418027,223975112), size = 0.5, color = "black", linetype = "dashed") + # q36.1
        geom_vline(xintercept = c(224378698,225010461), size = 0.5, color = "black", linetype = "dashed") + # q36.2
        geom_vline(xintercept = c(225399710,230058229), size = 0.5, color = "black", linetype = "dashed") + # q36.3
        # input_figureB[which(input_figureB$band == "q36.3"),c("genome_location")]
        scale_x_continuous(
                breaks = c(#1.50e+08, 
                        mean(c(168834132,176664676)),
                        mean(c(177138699,179101692)),
                        mean(c(179441982,182048822)),
                        mean(c(182140036,188297492)),
                        mean(c(188734155,190964358)),
                        mean(c(191029576,196289409)),
                        mean(c(196639554,202376936)),
                        mean(c(202635188,203936748)),
                        mean(c(204545793,208165343)),
                        mean(c(208236227,214411065)),
                        mean(c(214725646,219906502)),
                        mean(c(221418027,223975112)),
                        mean(c(224378698,225010461)),
                        mean(c(225399710,230058229))
                        #160555568
                        ),
                # labels = c("q31.1", "q31.2", "q31.3", "q32.1", "q32.2", "q32.3","q33.1","q33.2","q33.3","q34","q35","q36.1","q36.2","q36.3")
                labels = c("q31.1", "", "", "q32.1", "", "","q33.1","","","q34","","q36.1","","")
        ) +
        # scale_fill_gradient(name = paste0("Frequency of co-amp\n with SF3B1"), 
        #         low = coAmp_color[3], high = coAmp_color[8],
        #         labels = scales::percent_format(scale = 100),
        #         guide = guide_colorbar(barheight = 1, barwidth = 6)) +
        scale_colour_manual(values = c("#2015D0" = "blue", "#DCA620" = "#DCA620", "#ff0000" = "#ff0000", "#000000" = "black"),
                        labels = c("#DCA620" = "Candidates", "#2015D0" = "NMD-related", "#ff0000" = "Candidates NMD", "#000000" = "Other"),
                                name = "") +
        geom_text_repel(data = subset(input_figureD, (gene_symbols %in% genes_label) & (iNMDeff_type %in% c("ASE")) ),
                aes(label = gene_symbols), color = label_color, size = 2, arrow = arrow(length = unit(0.005, 'npc'), type = "closed"),
                        point.padding = 0.005, nudge_x = .01, nudge_y = .03, max.overlaps = 200) +
        scale_fill_brewer(palette = "Accent", labels = c("ETG","ASE"), direction = -1) +
        scale_y_continuous(breaks = c(0,-0.05,-0.1,-0.15, -0.2,-0.25, -0.3)) +
        # scale_x_break(breaks = c(121167646, 143955364))+
        ggplot_theme_bw() +
        theme(plot.margin = unit(c(-0.5, 0.5, 0, 0), "cm"),
                legend.position = "bottom") +
        guides(colour = guide_legend(override.aes = list(size = 4), nrow = 2),
               fill = guide_legend(override.aes = list(size = 4), nrow = 2) ) +
        facet_zoom(xlim = c(168834132-100, 230058229+100),
                        zoom.size = 0.5,
                        ylim = c(-0.17, -0.3),
                        show.area = TRUE,
                        horizontal = FALSE)

# Remove gene symbol labels and vertical dashes from the original panel
plot_D2 <- ggplot_build(plot_D2)
for (i in c(2:17)) {
        plot_D2$data[[i]][plot_D2$data[[i]]$PANEL %in% c(1,2), 'alpha'] <- 0
}

plot_D2 <- ggplot_gtable(plot_D2)
plot_D <- plot_grid(plot_D1, plot_D2, nrow = 2, ncol = 1,  
                        rel_heights = c(0.3,0.7), labels = c("",""), label_size = 20)

################################################
############ PANEL UP --> PLOTS A-B ############
################################################

fig_up <- plot_grid( 
        plot_A,
        plot_B,
        nrow = 1, labels = c("A", "B"), label_size = 20, 
        rel_widths = c(0.5,0.5), rel_heights = c(0.5,0.5))

##########################################################
############ PANEL BOTTOM --> PLOTS C-D ##################
##########################################################

fig_bottom <- plot_grid(plot_D, plot_C, nrow = 1, ncol = 2, label_size = 20, 
                labels = c("C","D"))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.4,0.6)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig21_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 240, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig21_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 240, units = "mm")
