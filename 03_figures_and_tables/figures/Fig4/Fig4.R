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
############ FIG 4A ############
################################

# Data
input_figureA <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4A_1.RData")

input_figureA$replicated_hits <- ifelse(input_figureA$replicated_hits %in% c("11","27","71"),"Not replicated",input_figureA$replicated_hits)
table(input_figureA$replicated_hits)

# Colors
# Extract colors from the Brewer palette
brewer_colors <- scales::brewer_pal(palette = "Dark2")(length(unique(input_figureA$replicated_hits)))
# Modify the color for the "other" group
names(brewer_colors) <- unique(input_figureA$replicated_hits)
brewer_colors[names(brewer_colors) == "Not replicated"] <- "#787B80"
# Manual change
# brewer_colors[names(brewer_colors) == "71"] <- "#BDCB14"
# CNA-PCs with different directionality
rem_CNAs <- c("pancancer_27","pancancer_71","pancancer_11")

NMD_method <- "ASE"
plot_A <- ggplot(input_figureA, aes(x = eval(parse(text=paste0(NMD_method,"_coefficient"))), 
              y = -log10(eval(parse(text=paste0(NMD_method,"_p_value")))), color = replicated_hits )) +
                  xlim(-0.5,0.5) +
                  geom_point(data = input_figureA %>% filter(replicated_hits == "Not replicated"),
                      size = 2, alpha = 0.3, color = "grey") + 
                #   geom_point(data = input_figureA %>% filter(replicated_hits != "Not replicated" & ! (replicated_hits_labels %in% rem_CNAs)),
                #       size = 2, alpha = 0.75) +
                  geom_point(data = input_figureA %>% filter(replicated_hits != "Not replicated"),
                      size = 2, alpha = 1) +
                  geom_label_repel(data = input_figureA %>% filter(replicated_hits != "Not replicated" & ! (replicated_hits_labels %in% rem_CNAs)),
                                aes(label=replicated_hits_labels),hjust=0.5, vjust=0.5, size = 3, 
                      max.overlaps = nrow(input_figureA), alpha = 1, key_glyph = draw_key_blank) +
                #   geom_text_repel(data = input_figureA %>% filter(replicated_hits != "Not replicated" & (replicated_hits_labels %in% rem_CNAs)),
                #                 aes(label=replicated_hits_labels),hjust=0.5, vjust=0.5, size = 3, 
                #       max.overlaps = nrow(input_figureA), alpha = 1, key_glyph = draw_key_blank) +
                  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
                  geom_vline(xintercept=0, linetype="dashed", color = "red") +
                  labs(x = "Association effect size", y =  expression("-log"[10] * "(p-value)"), color = "Replicated CNA-PCs",
                   title = paste0("CNA-PCs associations to ",NMD_method," iNMDeff")) +
                  scale_color_manual(values = brewer_colors) + # Use the modified colors
                  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
                  ggplot_theme() +
                  theme(legend.position="top",
                        legend.spacing.y = unit(0.05, "cm"),
                        legend.text = element_text(size = 9),
                        legend.key.height = unit(0.01, "cm")) +
                  guides(color = guide_legend(override.aes = list(size=4)))

################################
############ FIG 4B ############
################################

# Data
input_figureB <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4B.RData")

# Outlier
input_figureB <- input_figureB[-which(input_figureB$gene_name %in% c("SNORA44","RN7SL371P","RNA5SP162")),]
# input_figureB[which(input_figureB$chr_arm == "p" & input_figureB$values > 0.2),]

plot_B <- ggplot(data = input_figureB, aes(x = genome_location, y = values, group = ind, color = ind)) +
        geom_point(size = 1) + ggtitle(paste0("")) +
        geom_line() +  
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
        geom_hline(yintercept=c(0,1)) +
        labs(color = "CNA-PC3 groups", x = "Chr1 genome location", y = "GISTIC CNA average state") +
        ggplot_theme() +
        scale_color_brewer(palette="Set2") +
        theme(legend.position = "top",
                legend.text = element_text(size = 9),
                plot.margin = unit(c(0, 1, 0, 1), "cm")) +
        scale_x_discrete(guide=guide_axis(n.dodge = 3)) +
        guides(color = guide_legend(override.aes = list(size = 6), nrow = 1))

##################################
############ FIG 4C ##############
##################################

# Data
input_figureC <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4C.RData")
input_figureC <- input_figureC[order(input_figureC$genome_location),]
CNV_color <- brewer.pal(n = 2, name = "Paired")
CNV_color <- CNV_color[-3]
# Outlier
input_figureC <- input_figureC[-which(input_figureC$gene_symbols == "SNORA44"),]

# Plot C-1 (up)

# Create a fake color only to output the legend
input_figureC$CNA_type <- "CNA_freq_1"
input_figureC$CNA_type[1:20] <- "CNA_freq_2"

plot_C1 <- ggplot(data = input_figureC, aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1, fill = factor(CNA_type) )) +
                # geom_bar(stat="identity", width = 1.2) +
                # geom_point()+
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1), stat = "identity", color = CNV_color[1]) +
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_2), stat = "identity", color = CNV_color[2]) +
                labs(title = "", x = "", y = "CNA freq") +
                geom_hline(yintercept= 0, size = 0.5, color = "black") +
                ggplot_theme() +
                theme(axis.text.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        legend.text = element_text(size = 9),
                        plot.margin = unit(c(0, 0.5, -0.5, 0.25), "cm"),
                        legend.position = "top") +
                scale_fill_manual(values = c("CNA_freq_1" = CNV_color[1], "CNA_freq_2" = CNV_color[2]),
                                labels = c("CNA_freq_1" = "> 0", "CNA_freq_2" = ">= 1.5"),
                                  name = "GISTIC score")

# Plot C-2 (bottom)

# Colors and labels
coAmp_color <- brewer.pal(n = 8, name = "OrRd")
candidate_genes_NMD <- c("SMG5","RBM8A","SF3B4","INTS3")
NMD_related_genes <- c("SMG7")
genes_label <- c("SMG5","SMG7","RBM8A","SF3B4","INTS3","PMF1","BOLA1","CLK2","PRPF3","DENND4B","TPR","PI4KB","USF1","SETDB1","CRTC2","GORAB","ZNF687",
                "SCNM1","POGZ","YY1AP1","LYSMD1","DUSP12","USP21","GON4L","MSTO1","CHTOP","CHTOP","PPOX","VPS45","VPS72")
genes_label <- c("PMF1","GON4L","CRTC2","POGZ","PRPF3","VPS72")
genes_label <- c(candidate_genes_NMD,genes_label,NMD_related_genes)
input_figureC$label_color <- "#000000"
input_figureC[input_figureC$gene_symbols %in% candidate_genes_NMD,"label_color"] <- "#ff0000"
input_figureC[input_figureC$gene_symbols %in% NMD_related_genes,"label_color"] <- "#2015D0"
input_figureC[input_figureC$gene_symbols %in% genes_label & !input_figureC$gene_symbols %in% 
                        c(candidate_genes_NMD,NMD_related_genes),"label_color"] <- "#DCA620"
tmp <- subset(input_figureC, gene_symbols %in% genes_label)
label_color_char <- tmp[!tmp$gene_name %in% c("PMF1","GON4L"),"label_color"]
input_figureC$genome_location <- factor(input_figureC$genome_location, levels = unique(input_figureC$genome_location))
table(input_figureC$label_color)

plot_C2 <- ggplot(data = input_figureC, aes(x = as.numeric(as.character(genome_location)), y = eval(parse(text=paste0("NMDeff_endogenous_1"))), fill = as.numeric(coCNV_freq_1),
                        colour = factor(label_color))) +
        geom_point(size = 1.25, shape = 21,  stroke = 0.4, alpha = 0.7) +
        geom_text_repel(data = subset(input_figureC,  gene_symbols %in% genes_label[!genes_label %in% c("PMF1","GON4L")] ),
                aes(label = gene_symbols), color = label_color_char, size = 2.5, arrow = arrow(length = unit(0.005, 'npc'), type = "closed"),
                        point.padding = 0.005, nudge_x = .01, nudge_y = .003, max.overlaps = 200) +
        geom_text_repel(data = subset(input_figureC, gene_symbols %in% c("PMF1","GON4L")), fontface = "bold",
                aes(label = gene_symbols), color = c("#DCA620","#DCA620"), size = 2.5, arrow = arrow(length = unit(0.05, 'npc'), type = "closed"),
                        nudge_x = .09, nudge_y = -.006, max.overlaps = 200) +
        labs(title = "", x = "Chr 1 genome location", y = "median ETG iNMDeff") +
        geom_hline(yintercept= 0, size = 0.5, color = "black") +
        # 1q21.1-23.1 bands
        geom_vline(xintercept = c(143955364,147242641), size = 0.5, color = "black", linetype = "dashed") + # q21.1
        geom_vline(xintercept = c(147258885,150574551), size = 0.5, color = "black", linetype = "dashed") + # q21.2
        geom_vline(xintercept = c(150600539,155078872), size = 0.5, color = "black", linetype = "dashed") + # q21.3
        geom_vline(xintercept = c(155127460,156579727), size = 0.5, color = "black", linetype = "dashed") + # q22
        geom_vline(xintercept = c(156594487,158999968), size = 0.5, color = "black", linetype = "dashed") + # q23.1
        scale_x_continuous(
                breaks = c(1.50e+08, 
                        mean(c(150600539,155078872)),
                        mean(c(155127460,156579727)),
                        mean(c(156594487,158999968)),
                        160555568),
                labels = c("q21.2", "q21.3", "q22", "q23.1", "q23.2")
        ) +
        scale_fill_gradient(name = paste0("Co-amp\n with SMG5"), 
                low = coAmp_color[3], high = coAmp_color[8],
                labels = scales::percent_format(scale = 100),
                guide = guide_colorbar(barheight = 1, barwidth = 6)) +
        scale_colour_manual(values = c("#2015D0" = "blue", "#DCA620" = "#DCA620", "#ff0000" = "#ff0000", "#000000" = "black"),
                        labels = c("#DCA620" = "Candidates", "#2015D0" = "NMD-related", "#ff0000" = "Candidates NMD", "#000000" = "Other"),
                                name = "") +
        scale_y_continuous(breaks = c(0,-0.02,-0.05,-0.07, -0.1)) +
        scale_x_break(breaks = c(121167646, 143955364))+
        ggplot_theme_bw() +
        theme(plot.margin = unit(c(-0.5, 0.5, 0, 0), "cm"),
                legend.position = "bottom"
                #axis.title.y = element_text(margin = margin(t = 0, r = -5, b = 0, l = 0))
                ) +
        guides(colour = guide_legend(override.aes = list(size = 6), nrow = 2)) +
        facet_zoom(xlim = c(149887890-100, 161749758+100),
                        zoom.size = 0.5,
                        ylim = c(-0.08, -0.1),
                        show.area = TRUE,
                        horizontal = FALSE)

# Remove gene symbol labels and vertical dashes from the original panel
plot_C2 <- ggplot_build(plot_C2)
# plot_C2$data[[1]][plot_C2$data[[1]]$PANEL %in% c(1,2) & plot_C2$data[[1]]$x %in% c(160555568), 'alpha'] <- 0
for (i in c(2,3,4,5,6,7,8,9)) {
        plot_C2$data[[i]][plot_C2$data[[i]]$PANEL %in% c(1,2), 'alpha'] <- 0
}

plot_C2 <- ggplot_gtable(plot_C2)
plot_C <- plot_grid(plot_C1, plot_C2, nrow = 2, ncol = 1,  
                        rel_heights = c(0.3,0.7), labels = c("",""), label_size = 20)

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, plot_C, width = 250, height = 200, units = "mm") 

##################################
############ FIG 4D ##############
##################################

# Data
input_figureD <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4D.RData")
# Colors
input_figureD$label_color <- "#8A8A8A"
input_figureD[input_figureD$gene_symbol %in% "SMG5","label_color"] <- "#2015D0"
tmp <- subset(input_figureD, gene_symbol %in% top_hits)
label_color <- tmp$label_color
input_figureD$controls_vs_rest <- NA
input_figureD[which(input_figureD$selected_bands == "other"),"controls_vs_rest"] <- "other"
input_figureD[which(input_figureD$selected_bands != "other"),"controls_vs_rest"] <- "1q21-23.1"
control_bands <- c("q24.1","q24.2","q24.3","q25.1","q25.2","q25.3","q31.1","q31.2","q31.3","q32.1","q32.2","q32.3","q41","q42.11","q42.12","q42.13","q42.2",
                "q42.3","q43","q44")
selected_bands_char <-  c("q21.1","q21.2","q21.3","q22","q23.1","q23.2","q23.3")
selected_bands_2_char <- c("q21.3")
input_figureD <- input_figureD %>%
        mutate(selected_bands_2 = if_else(band %in% control_bands, "1q-end",band)) %>%
        mutate(selected_bands_2 = if_else(selected_bands_2 %in% selected_bands_char, "1q21.1-23.3",selected_bands_2)) #%>%
        # mutate(selected_bands_2 = if_else(selected_bands_2 %in% selected_bands_2_char, "q21.3",selected_bands_2)) %>%
        # mutate(selected_bands_2 = if_else(selected_bands_2 %in% c("1q21.1-23.1","1q-end"), selected_bands_2,"NA"))
table(input_figureD$selected_bands_2)
# Extract colors from the Brewer palette
brewer_colors <- scales::brewer_pal(palette = "Set2")(length(unique(input_figureD$selected_bands_2)))
# Modify the color for the "other" group
names(brewer_colors) <- unique(input_figureD$selected_bands_2)
brewer_colors[names(brewer_colors) == "other"] <- "#808080"
# Manual change
brewer_colors[names(brewer_colors) == "1q-end"] <- "#181A1D"
brewer_colors[names(brewer_colors) == "1q21.1-23.3"] <- "#DCA620"

plot_D <- ggplot(data = input_figureD, aes(x = CNA_RNA_corr, y = ETG_iNMDeff_RNA_corr, 
                        color = factor(selected_bands_2)), fill = factor(selected_bands_2)) +
        annotate("rect", xmin = 0.25, xmax = Inf, ymin = -Inf, ymax = -0.15, fill= "#FCF1BC") + 
        annotate("text", x = 0.47, y = -0.33, label = "Candidate genes") + 
        stat_density_2d(data = input_figureD %>% filter(selected_bands_2 == "1q-end"),
                        aes(fill = stat(level)), geom = "polygon", alpha = 0.5) +
        geom_smooth(aes(group = factor(selected_bands_2)),fill = "black", method = lm, se = FALSE) +
        geom_point(data = input_figureD %>% filter(selected_bands_2 != "1q-end"), size = 2, alpha = 0.35) + 
        geom_point(data = input_figureD %>% filter(selected_bands_2 == "1q-end"), shape = 3, size = 2, alpha = 0.5) + 
        geom_vline(xintercept= 0.25, size = 0.5, color = "black", linetype = "dashed") +
        geom_hline(yintercept= -0.15, size = 0.5, color = "black", linetype = "dashed") +
        labs(title = "Chr 1q genes", color = "Chr1 bands", x = "Gene exp vs CNA correlation", y = "Gene exp vs ETG iNMDeff correlation") +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = c(0,0.1,0.2,0.25,0.3,0.4,0.5)) +
        scale_y_continuous(breaks = c(0.1,0,-0.1,-0.15,-0.2,-0.3)) +
        # geom_text_repel(data = subset(input_figureA, gene_symbol %in% top_hits & !gene_symbol %in% c("SMG5")),
        #                 max.overlaps = nrow(input_figureA), #color = label_color,
        #                 aes(label = gene_symbol), size = 3,      arrow = arrow(length = unit(0.015, 'npc')),
        #                 point.padding = 0.05, nudge_x = .04, nudge_y = -.03) +
        geom_text_repel(data = subset(input_figureD, gene_symbol %in% c("GON4L","VPS72","PRPF3","PMF1","POGZ","CHTOP")),
                        max.overlaps = nrow(input_figureD), #color = label_color,
                        aes(label = gene_symbol), size = 3,      arrow = arrow(length = unit(0.015, 'npc')),
                        point.padding = 0.05, nudge_x = .04, nudge_y = -.03) +
        geom_text_repel(data = subset(input_figureD, gene_symbol %in% c("SMG5","INTS3","SF3B4","RBM8A")),
                        max.overlaps = nrow(input_figureD), color = "#ff0000",
                        aes(label = gene_symbol), size = 3, arrow = arrow(length = unit(0.01, 'npc')),
                        point.padding = 0.05, nudge_x = .01, nudge_y = -.01) +
        #scale_color_brewer(palette = "Set2", direction = -1, values = c(other = "grey")) +
        scale_color_manual(values = brewer_colors) + # Use the modified colors
        ggplot_theme()+ 
        #   stat_poly_eq(aes(label = paste(stat(rr.label), stat(p.value.label), sep = "~~~")), 
        #        formula = y ~ x, parse = TRUE, label.x = "right") +
        ggpubr::stat_cor(method = "pearson", label.y.npc = "bottom") +
        theme(legend.position = "bottom",
                legend.text = element_text(size = 9))

# Plot D2
input_figureD2 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4D_2.RData")

input_figureD2_stack <- stack(input_figureD2[,c("SMG5_RNAseq","PMF1_RNAseq")])
input_figureD2_stack$endogenous_NMD_Consensus <- rep(input_figureD2$endogenous_NMD_Consensus,2)
colnames(input_figureD2_stack) <- c("RNAseq","gene","ETG")
input_figureD2_stack$gene <- gsub("_RNAseq","",input_figureD2_stack$gene)

plot_D2 <- ggplot(data = input_figureD2_stack, aes(x = RNAseq, y = ETG)) +
                facet_wrap(.~gene) +
                geom_point(size = 1, alpha = 0.2, color = "#3b528b",  shape = 21) +
                stat_density_2d(aes(fill = ..level..), geom = "polygon") + guides(fill = FALSE) +
                # geom_density_2d()
                # geom_bin2d() + 
                geom_smooth(color = "#3b528b", method = lm, se = TRUE) +
                #coord_cartesian(ylim = c(-5,5), xlim = c(2,8)) +
                coord_cartesian(ylim = c(-5,5), xlim = c(0,150)) +
                ggplot_theme_bw()+ 
                scale_x_continuous(breaks = c(50,150)) +
                labs(title = "", x = "Gene expression (TPM)", y = "ETG iNMDeff") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom",
                        axis.text.x = element_text(size = 8),
                        strip.background = element_blank())

# Plot D3
input_figureD2$gene_1 <- "PMF1"
plot_D3 <- ggplot(data = input_figureD2, aes(x = PMF1_RNAseq, y = PMF1_CNA)) +
                facet_wrap(.~gene_1) +
                # geom_point(size = 2, alpha = 0.15, color = "#B03838",  shape = 21) +
                # stat_density_2d(aes(fill = ..level..), geom = "polygon", na.rm = TRUE, bins = 30)+
                coord_cartesian(ylim = c(-1,4), xlim = c(0,200)) +
                ggplot_theme_bw()+ 
                geom_bin2d(bins = 50, alpha = 0.7) + guides(fill = FALSE) +   
                geom_smooth(method = lm, se = TRUE, color = "black") +
                scale_fill_gradient(low = "#B03838", high = "#FF0000") +
                # scale_fill_continuous(type = "#B03838") +
                # scale_x_continuous(breaks = c(50,150)) +
                labs(title = "", x = "", y = "CNA") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom",
                        plot.margin = unit(c(0, 0, 0, 0), "cm"),
                        strip.background = element_blank(),
                        axis.text.x = element_text(size = 8),
                        axis.title.y = element_text(margin = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt")))
# Plot D4
input_figureD2$gene_2 <- "SMG5"
plot_D4 <- ggplot(data = input_figureD2, aes(x = SMG5_RNAseq, y = SMG5_CNA)) +
                facet_wrap(.~gene_2) +
                # geom_point(size = 2, alpha = 0.15, color = "#B03838",  shape = 21) +
                coord_cartesian(ylim = c(-1,4), xlim = c(0,200)) +
                ggplot_theme_bw()+ 
                geom_bin2d(bins = 50,  alpha = 0.7) + guides(fill = FALSE) +   
                geom_smooth(method = lm, se = TRUE, color = "black") +
                scale_fill_gradient(low = "#B03838", high = "#FF0000") +
                # scale_x_continuous(breaks = c(50,150)) +
                labs(title = "", x = "", y = "") +
                ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 2, label.y = 2.5, label.sep = "\n", size = 3)+
                theme(legend.position = "bottom",
                        axis.title.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.text.x = element_text(size = 8),
                        strip.background = element_blank(),
                        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

##################################
############ FIG 4E ##############
##################################

# Plot E-1
input_figureE_1 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4E.RData")
combinations_E_1 <- combn(names(table(input_figureE_1$bins)), 2, simplify = FALSE)

plot_E_1 <- input_figureE_1 %>% 
        filter(NMD_method == "Endogenous") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "TCGA", x = "CNA-PC3 groups", y = "ETG iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                        plot.margin = unit(c(1, 0, 1, 0), "cm")) +
                stat_compare_means(comparisons = combinations_E_1, size = 4,
                                label.y = c(-0.25,0.25,0.75),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

# Plot E-2
input_figureE_2 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4E_2.RData")
combinations_E_2 <- combn(names(table(input_figureE_2$Chr_1q_CNV)), 2, simplify = FALSE)

plot_E_2 <- ggplot(input_figureE_2, aes(x = Chr_1q_CNV, y = -(NMDeff), fill = Chr_1q_CNV)) +
                geom_violin() + scale_fill_brewer(palette = "Dark2", direction = 1) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                scale_x_discrete(labels = c("Gain","Neutral","Deletion")) +
                labs( title = "Validation in cell lines",x = "CNA state - chr 1q", y = "TCGA-model ETG cNMDeff") +
                ggplot_theme() + coord_cartesian(ylim = c(-2,2)) +
                theme(legend.position = "none",
                        plot.title = element_text(hjust = 0.90),
                        plot.margin = unit(c(1, 0, 1, 0), "cm")) +
                stat_compare_means(comparisons = combinations_E_2, size = 4,
                                label.y = c(0.25,0.75,1.25),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE) +
                annotate("text",
                        #fontface = "bold",
                        x = 1:length(table(input_figureE_2$Chr_1q_CNV)),
                        y = aggregate( -(NMDeff) ~ Chr_1q_CNV, input_figureE_2, median)[ , 2],
                        label = table(input_figureE_2$Chr_1q_CNV),
                        col = "black",
                        vjust = - 0.5,
                        size = 4)

plot_E <- plot_grid(plot_E_1, plot_E_2, nrow = 1, ncol = 2)

##################################
############ FIG 4F ##############
##################################

input_figureF <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/Fig4F.RData")

# NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
candidate_genes <- rownames(input_figureF[input_figureF$label_colors == "orange",])
candidate_genes_nonNMD <- candidate_genes[!candidate_genes %in% NMD_genes$gene_symbol]
candidate_genes_NMD <- c("SMG5","RBM8A",candidate_genes[candidate_genes %in% NMD_genes$gene_symbol])
control_genes <- rownames(input_figureF[input_figureF$label_colors == "blue",])
# Remove some genes
NMD_related_genes_to_remove <- NMD_genes[NMD_genes$NMD_type %in% c("NMD_related","NMD_ER"),"gene_symbol"]
NMD_factors <- NMD_genes[!NMD_genes$NMD_type %in% c("NMD_related","NMD_ER"),"gene_symbol"]
NMD_factors <- NMD_factors[!NMD_factors %in% "SMG5"]

# Filter selected genes
genes_to_keep <- c(NMD_factors,candidate_genes_NMD,c("CHTOP","POGZ","TRP","PRPF3","GON4L","VPS72","PMF1"))
input_figureF <- input_figureF[rownames(input_figureF) %in% genes_to_keep, 
                colnames(input_figureF) %in% c("label_colors",genes_to_keep)]

# Colors
input_figureF$gene <- rownames(input_figureF)
input_figureF$category <- NA

# Assign categories based on gene membership
input_figureF <- input_figureF %>%
  mutate(
    category = case_when(
      gene %in% candidate_genes_nonNMD ~ "Candidates",
      gene %in% candidate_genes_NMD ~ "Candidates NMD",
      gene %in% NMD_factors ~ "NMD-related",
#       gene %in% control_genes ~ "Controls",
      TRUE ~ "Unknown"
    ),
    label_colors = case_when(
      category == "Candidates" ~ "orange",
      category == "Candidates NMD" ~ "red",
      category == "NMD-related" ~ "blue",
#       category == "Controls" ~ "black",
      TRUE ~ "gray"  # Assign gray for any unknown genes
    )
  )
table(input_figureF$label_colors,input_figureF$category)

plot_F <- ggcorrplot(as.matrix(input_figureF[,!colnames(input_figureF) %in% c("label_colors","gene","category","color")]), method = "square", 
                        hc.order = TRUE, type = "full", lab = FALSE,
                        title = "CRISPR KO co-dependency scores",
                        ggtheme = theme_bw(),
                        tl.cex = 8,
                        show.legend = TRUE,
                        tl.srt= 55,
                        tl.col = input_figureF$label_colors,
                        show.diag = TRUE)

# Extract the data from the ggplot object
plot_data <- levels(plot_F$data$Var1)

ordered_colors <- input_figureF$label_colors[match(levels(plot_F$data$Var1), rownames(input_figureF))]

plot_F$data$label_colors <- ordered_colors 
plot_F$data$label_colors <- ifelse(plot_F$data$label_colors == "blue","NMD-related",plot_F$data$label_colors)
# plot_F$data$label_colors <- ifelse(plot_F$data$label_colors == "black","Controls",plot_F$data$label_colors)
plot_F$data$label_colors <- ifelse(plot_F$data$label_colors == "orange","Candidates",plot_F$data$label_colors)
plot_F$data$label_colors <- ifelse(plot_F$data$label_colors == "red","Candidates NMD",plot_F$data$label_colors)

vector_colors <- c("#0000FF","#FFA500","#FF0000")
names(vector_colors) <- c("NMD-related","Candidates","Candidates NMD")

plot_F <- plot_F +
        geom_point(aes(color = factor(label_colors)), alpha = 0) +
        scale_color_manual(limits = names(vector_colors), values = vector_colors)+
        scale_fill_gradient2(limits = c(0,0.1), low = "#2A0BEE", mid = "#EDFEFE", high = "#EE0B0B") +
        labs(fill = "Scores", x = "", y = "", color = "Genes") + ggplot_theme() +
        guides(
                color = guide_legend(override.aes = list(size = 4, alpha = 1)),
                fill = guide_colourbar(
                        barwidth = 0.75,  # Adjust the width of the bar
                        barheight = 5  # Adjust the height of the bar
                )) +
        theme(plot.title = element_text(hjust = 0.5),
                plot.background = element_rect("white"),
                # plot.margin = unit(c(0.1, 0.95, -0.1, 1.35), "cm"),
                # plot.margin = unit(c(-0.5, 0.5, -0.5, 0.5), "cm"),
                legend.position = "right",
                axis.text.x = element_text(color = ordered_colors, angle = 90, hjust = 1, vjust = 0.5, size = 7),
                axis.text.y = element_text(color = ordered_colors, size = 7))

# plot_F_final <- plot_grid(NULL,plot_F, nrow = 2, ncol = 1,  rel_heights = c(0,1))

# Manual for the manuscript
# input_figureF["PMF1","SMG5", drop = FALSE]
# input_figureF["GON4L","SMG5", drop = FALSE]
# input_figureF["PMF1","SMG6", drop = FALSE]
# input_figureF["GON4L","SMG6", drop = FALSE]
# input_figureF["PMF1","SMG7", drop = FALSE]
# input_figureF["GON4L","SMG7", drop = FALSE]

# input_figureF["VPS72","SMG5", drop = FALSE]
# input_figureF["PRPF3","SMG5", drop = FALSE]

###############################################
############ PANEL UP --> FIG 4A-B ############
###############################################

fig_up <- plot_grid( 
        plot_A,
        plot_B,
        nrow = 1, labels = c("A", "B"), label_size = 20, 
        rel_widths = c(0.5,0.5), rel_heights = c(0.5,0.5))

#################################################
############ PANEL LEFT --> FIG 4C-D ############
#################################################

# plot_E_final <- plot_grid(NULL,plot_Es, nrow = 2, ncol = 1,  rel_heights = c(0,1))

plot_D2_3_4 <- plot_grid(plot_D2, plot_D3, plot_D4, nrow = 1, ncol = 3, rel_widths = c(0.5,0.27,0.23), align = "hv")
plot_D_tmp <- plot_grid(plot_D, plot_D2_3_4, nrow = 2, ncol = 1, rel_heights = c(0.7,0.3))
plot_D_tmp <- ggdraw() + 
  draw_plot(plot_D_tmp) +
  draw_label("Gene expression (TPM)", x = 0.76, y = 0.03, hjust = 0.5, size = 12)

fig_left <- plot_grid(plot_E, plot_D_tmp, nrow = 2, ncol = 1,  
                        rel_heights = c(0.35,0.65), labels = c("C","E"), label_size = 20)

##################################################
############ PANEL RIGHT --> FIG 4C-D ############
##################################################

# fig_bottom <- plot_grid(plot_E, plot_F, nrow = 1, ncol = 2, # align = "hv",  
#                         rel_widths = c(0.5,0.5), labels = c("E","F"), label_size = 20)

fig_right <- plot_grid(plot_C, plot_F, nrow = 2, ncol = 1, # align = "hv",  
                        rel_heights = c(0.65,0.35), labels = c("D","F"), label_size = 20)

#############################################
############ FINAL FIG --> FIG 4 ############
#############################################

# Fig4_final <- plot_grid(fig_up, fig_middle, fig_bottom, nrow = 3, rel_heights = c(0.3,0.5,0.3))
Fig4_final_tmp <- plot_grid(fig_left, fig_right, nrow = 1, ncol = 2)
Fig4_final <- plot_grid(fig_up, Fig4_final_tmp, nrow = 2, rel_heights = c(0.3,0.7))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig4_complete.png"
ggsave(final_figure_path, Fig4_final, width = 250, height = 330, units = "mm", device = "png")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig4_complete.pdf"
ggsave(final_figure_path, Fig4_final, width = 250, height = 330, units = "mm", device = "pdf")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig4_complete.tiff"
ggsave(final_figure_path, Fig4_final, width = 250, height = 330, units = "mm", device = "tiff")

