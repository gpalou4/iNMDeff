# conda activate R_figures

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ FIG 2A ############
################################

# Data
input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2A.RData")
# Change stuff
NMD_genesets <- c("RandomGenes\nwith NMD\nfeatures","RandomGenes\nwithout NMD\nfeatures","NMD Consensus", "NMD All")
input_figureA <- input_figureA %>%
        filter(ind %in% NMD_genesets)
input_figureA$ind <- gsub("RandomGenes\nwith NMD\nfeatures","RandomGenes\nw/ NMD feat.",input_figureA$ind)
input_figureA$ind <- gsub("RandomGenes\nwithout NMD\nfeatures","RandomGenes\nw/o NMD feat.",input_figureA$ind)
# combinations <- list(c("NMD All","NMD Consensus"), c("RandomGenes\nwith NMD\nfeatures","NMD Consensus"),c("RandomGenes\nwithout NMD\nfeatures","NMD Consensus"))
# combinations <- list(c("NMD All","NMD Consensus"), c("RandomGenes\nw/ NMD feat.","NMD Consensus"),c("RandomGenes\nw/o NMD feat.","NMD Consensus"))
combinations <- list(c("NMD All","NMD Consensus"),c("RandomGenes\nw/o NMD feat.","NMD Consensus"))
tmp_df <- input_figureA[!is.na(input_figureA$values),] # NAs
ylim <- c(-1,4)
annotate_hjust <- 0.5
pval_ylabel <- c(2.35,2.85,3.35)
boxplot_width <- 0.3

input_figureA[input_figureA$ind == "RandomGenes\nw/ NMD feat.","NMD_geneset"] <- "NMD"

# Plot
plot_A <- ggplot(data = input_figureA, aes(x = factor(ind), y = values, fill = factor(NMD_geneset, levels = c("RandomGenes","NMD")))) +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_cartesian(ylim = ylim) +
        geom_boxplot(width = boxplot_width, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(y = "\nETG iNMDeff", x = "", title = "TCGA", fill = "Gene set") +
        ggplot_theme() + 
         scale_color_manual(values = c("Random" = "blue", "NMD" = "red"),
                       labels = c("Random", "NMD Genes")) +
        theme(axis.text.x = element_text(angle = 20, hjust = 0.65, vjust = 0.75, size = 8, lineheight = 0.65),
              legend.margin = margin(t = -20, r = -0, b = 0, l = -0, unit = "pt"),
              legend.title = element_blank(),
              legend.text = element_text(size = 7),
        #       axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
              plot.margin = unit(c(0.6, 0, 0.6, -0.3), "cm"),
              legend.position = "bottom") +
        guides(fill = guide_legend(keyheight = unit(0.5, "cm"), keywidth = unit(0.5, "cm"))) +
        scale_fill_brewer(palette = "Dark2", labels = c("Negative control","NMD gene sets")) +
        annotate("text",
                x = 1:length(table(tmp_df$ind)),
                y = aggregate( values ~ ind, tmp_df, median)[ , 2],
                label = table(tmp_df$ind),
                col = "black",
                hjust = annotate_hjust,
                vjust = 3,
                size = 3) +
        stat_compare_means(comparisons = combinations, size = 5, vjust = 0.5,
                      label.y = pval_ylabel,
                      label = "p.signif", method = "wilcox.test", hide.ns = TRUE)

aggregate(values ~ ind + NMD_geneset, data = input_figureA, median)

################################
############ FIG 2B ############
################################

# Data
input_figureB <- read.table(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2B.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#Change stuff
# combinations <- list(c("PTC NMD-triggering","Synonymous"),c("RandomGenes\nw/o NMD feat.","NMD Consensus"))
ylim <- c(-1,2)
pval_ylabel <- c(0.5,1,1.5)
boxplot_width <- 0.12
annotate_hjust <- 0.5
tmp_df <- input_figureB[!is.na(input_figureB$values),] # NAs
input_figureB$NMD_geneset <- ifelse(input_figureB$NMD_geneset == "PTC NMD-triggering","NMD-triggering PTCs",input_figureB$NMD_geneset)
input_figureB$NMD_geneset <- ifelse(input_figureB$NMD_geneset == "PTC NMD-evading","NMD-evading PTCs",input_figureB$NMD_geneset)
input_figureB$ind <- as.character(input_figureB$ind)
input_figureB$ind <- ifelse(input_figureB$ind == "PTC NMD-triggering","NMD-triggering PTCs",input_figureB$ind)
input_figureB$ind <- factor(ifelse(input_figureB$ind == "PTC NMD-evading","NMD-evading PTCs",input_figureB$ind))
input_figureB$NMD_geneset <- factor(input_figureB$NMD_geneset, levels = c("Synonymous","NMD-triggering PTCs","NMD-evading PTCs"))
combinations <- combn(names(table(input_figureB$ind)), 2, simplify = FALSE)
combinations <- combinations[-2]

# Plot
plot_B <- ggplot(data = input_figureB, aes(x = ind, y = values, fill = NMD_geneset)) +
        geom_violin(draw_quantiles = TRUE, scale = "width", #bw = 0.1,
                        lwd = 0.25, na.rm = TRUE, width = 0.7) +
        geom_boxplot(width = boxplot_width, 
                        color="black", outlier.shape = NA, alpha=0.1) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(y = "ASE iNMDeff", x = "", title = "TCGA") +
        ggplot_theme() + coord_cartesian(ylim = ylim) +
        theme(axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.75, size = 8),
                legend.position = "none",
                plot.margin = unit(c(0.6, 0, 0.6, 0), "cm"),) +
        scale_fill_brewer(palette = "Dark2") +
        annotate("text",
                x = 1:length(table(tmp_df$ind)),
                y = aggregate( values ~ ind, tmp_df, median)[ , 2],
                label = table(tmp_df$ind),
                col = "black",
                hjust = annotate_hjust,
                vjust = 1.75,
                size = 3) +
        stat_compare_means(comparisons = combinations, size = 5, vjust = 0.5,
                        label.y = pval_ylabel,
                        label = "p.signif", method = "wilcox.test", hide.ns = TRUE)

#################################
############ PANEL C ############
#################################


input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2C.RData")
input_figureE$ind <- gsub("NMD Consensus", "NMD \n Consensus",input_figureE$ind)
#input_figureE$ind <- factor(input_figureE$ind, levels = c("NMD \n Consensus","RP9P", "GAS5", "SMG5"))
list_plots <- list() 
selected_genes <- c("RP9P","GAS5")

for (variable in c("NMD \n Consensus",selected_genes)) {

        plot_title <- ""
        face_text <- "italic"
        X_title_size <- 0
        size_title <- 11
        input_figureE_filt <- input_figureE %>%
                                filter(ind == variable)
        if (variable == "SMG5") {
                ylim <- c(4,10.5)
                X_title_size <- 12
        } else if (variable == "SMG1") {
                ylim <- c(3,6.5)
        } else if (variable == "RP9P") {
                ylim <- c(1,6)
        } else if (variable == "GAS5") {
                ylim <- c(6,19)
        } else if (variable == "SMG6") {
                ylim <- c(1,4)
        } else if (variable == "SMG7") {
                ylim <- c(4,7)
        } else if (variable == "SMG8") {
                ylim <- c(1.5,4)
        } else if (variable == "SMG9") {
                ylim <- c(2,5)
        } else if (variable == "UPF1") {
                ylim <- c(5,9)
        } else if (variable == "UPF2") {
                ylim <- c(3,6)
        } else if (variable == "UPF3B") {
                ylim <- c(3,5)
        } else if (variable == "UPF3A") {
                ylim <- c(4,8)
        } else {
                ylim <- c(-4,4)
                plot_title <- "GTex"
                face_text <- "plain"
                size_title <- 7
        }

        p <- ggplot(data = input_figureE_filt, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
                geom_smooth(data = subset(input_figureE_filt, ind %in% selected_genes), 
                                aes(y=values, group = factor(NMDeff_type, levels = c("Low","High"))), 
                                method = "gam", se = TRUE, linewidth = 1.5) +
                # ggrastr::rasterize(geom_point(alpha=0.1, size = 0.05, shape = 1)) +
                geom_point(alpha=0.1, size = 0.05, shape = 1) +
                coord_cartesian(ylim = ylim)+
                facet_grid(ind ~ ., scales = "free_y") +
                labs(x = "Individuals sorted by iNMDeff", color = "ETG iNMDeff", 
                                y = "", title = plot_title) +
                ggplot_theme() +        
                theme(plot.title = element_text(hjust = 0.5),
                        panel.spacing = grid::unit(0.3, "lines"),
                        axis.text.x = element_blank(),
                        axis.title.x = element_text(size = 10),
                        strip.text = element_text(colour = "black", size = size_title, face = face_text),
                        plot.margin = unit(c(-0.55, 0, 0, 0), "cm"),
                        strip.background = element_rect()) +
                scale_color_manual(  
                        values = c("Low" = "#2D3263", "High" = "#F07626"), 
                        labels = c("Low", "High")
                        )   

        if (!variable %in% "GAS5") {
                p <- p + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = "none") +
                                guides(x = "none")     
        }

      list_plots[[length(list_plots) + 1]] <- p
}

plot_C <- ggarrange(list_plots[[1]],list_plots[[2]],list_plots[[3]],
                labels = c("","","",""), align = "v",
                font.label = list(size = 20), heights = c(0.17,0.22,0.42),
                ncol=1, nrow=length(list_plots), common.legend = FALSE)

# Add custom titles using grid

plot_C <- ggdraw(plot_C) +
  draw_label("iNMDeff", x = 0.04, y = 0.88, hjust = 0.5, vjust = 0.05, angle = 90, size = 11) +
  draw_label("Gene expression (TPM)", x = 0.04, y = 0.5, hjust = 0.5, vjust = 0.05, angle = 90, size = 11)

################################
############ FIG 2D ############
################################

# Data
input_figureD <- read.table(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2D.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
input_figureD <- input_figureD %>%
                        filter(NMD_gene_excluded == "all")

# Calculate mean and CI for each group (Type and CL)
summary_df <- input_figureD %>%
  group_by(Type, CL) %>%
  summarize(
    nb_coeff = nb_coeff,
    mean_nb_coeff = mean(nb_coeff),
    sd_nb_coeff = sd(nb_coeff),
    n = n(),
    se_nb_coeff = sd_nb_coeff / sqrt(n),   # Standard error
#     CI_lower = mean_nb_coeff - qt(0.975, df = n - 1) * se_nb_coeff, # 95% CI lower
#     CI_upper = mean_nb_coeff + qt(0.975, df = n - 1) * se_nb_coeff  # 95% CI upper
    CI_lower = mean(CI_2.5),
    CI_upper = mean(CI_97.5)
  )

plot_D <- ggplot(input_figureD, aes(x = CL, y = nb_coeff, fill = Type)) +
    geom_bar(
    data = input_figureD %>%
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
        label.y = max(input_figureD$nb_coeff) + 0.5  # Adjust for p-value placement
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
        legend.position = "top"
    )

################################
############ FIG 2E ############
################################

# Data
input_figureE <- read.table(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2E.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
input_figureE$allele <- factor(input_figureE$allele, levels = c("refCount_perc","altCount_perc"))

plot_E <- input_figureE %>%
    ggplot(aes(y = ASE, x = variant, fill = allele, color = TCGA_sample_type)) + 
    geom_bar(position="stack", stat="identity", size = 0.5) + 
    theme_classic() + ggplot_theme() +
    labs(title = "", x = "PTC variants", y = "Variant frequencies (%)", 
        fill = paste0("PTC"), color = "ASE\niNMDeff") + 
    scale_fill_brewer(labels = c("WT", "MUT"), palette = "Paired", direction = -1) +
    scale_color_manual(  
            values = c("Low" = "#2D3263", "High" = "#F07626"), 
            labels = c("High", "Low")
            ) +
        scale_x_continuous(
                breaks = c(2.5,9.5),
                labels = c("TCGA-F4-6856","TCGA-A6-6781")
        ) +
    theme(plot.title = element_text(size = 13),
                plot.margin = unit(c(0.5, 0, 0.25, 0), "cm"),
                legend.title = element_text(colour = "black", size = 11),
                # legend.box = "horizontal",
                legend.position = "right",
                legend.box = "vertical",
                legend.margin = margin(0,10,0,0),
                axis.text.x = element_text(size = 7)
        )


################################
############ PLOT F ############
################################

input_figureF <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2F.RData")

my_colors <- c(
  "Other" = "#F0027F",
  "Pan-GI" = "#FDC086",
  "Pan-kidney" = "#7FC97F",
  "Pan-nervous" = "#BEAED4",
  "Pan-reproductive" = "#FFFF99",
  "Pan-squamous" = "#386CB0"
)

# Organ system
input_figureF$pan_organ_system <- "Other"
input_figureF[input_figureF$tissue %in% c("KDNCTX","KDNMDL"),"pan_organ_system"] <- "Pan-kidney"
input_figureF[input_figureF$tissue %in% c("CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","STMACH","SNTTRM"),"pan_organ_system"] <- "Pan-GI"
input_figureF[input_figureF$tissue %in% c("BREAST","UTERUS","TESTIS","PRSTTE","VAGINA","OVARY","UCS","CVXECT","CVSEND","FLLPNT"),"pan_organ_system"] <- "Pan-reproductive"
input_figureF[input_figureF$tissue %in% c("NERVET","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG"),"pan_organ_system"] <- "Pan-nervous"
input_figureF[input_figureF$tissue %in% c("LUNG","BLDDER"),"pan_organ_system"] <- "Pan-squamous"
table(input_figureF$pan_organ_system)

formula <- "y ~ x"
plot_F <- ggplot(data = input_figureF, mapping = aes(x = GTEx_ETG_iNMDeff_method_ranking, y = GTEx_ASE_iNMDeff_method_ranking)) +
    geom_point(aes(color = pan_organ_system), alpha = 1) +
    geom_text_repel(aes(label = tissue), 
                    color = "black",  # Set color outside of aes()
                    size = 2, nudge_y = 0.05, max.overlaps = nrow(input_figureF)) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
    labs(x = "ETG iNMDeff GTex tissues ranking", y = "ASE iNMDeff GTex tissues ranking", 
            title = "", color = "") +
    ggplot_theme() +
    scale_color_manual(values = my_colors) + 
    theme(legend.position='top',
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9)
    #   plot.margin = unit(c(0, 2, 0, 2), "cm")
      ) + 
    guides(color = guide_legend(override.aes = list(size = 4))) +  # Increase point size in the legend
    stat_cor(p.accuracy = 0.0000000001, r.accuracy = 0.01, size = 4, label.y = 62)

###############################################
############ PANEL UP --> FIG 2A-B ############
###############################################

# fig_up <- plot_grid(plot_A, plot_B, nrow = 1, labels = c("A", "B"), label_size = 20, rel_widths = c(0.6,0.4), vjust = 1)
fig_up <- plot_grid(plot_A, plot_B , nrow = 1, ncol = 2, labels = c("A","B"), 
                        rel_widths = c(0.55,0.45), label_size = 20, vjust = 1)

################################################
############ PANEL MID --> FIG 2C-D ############
################################################

fig_mid <- plot_grid(plot_C, plot_D, nrow = 1, labels = c("C","D"), 
                                label_size = 20, rel_widths = c(0.45,0.55),
                                vjust = 0)

###################################################
############ PANEL BOTTOM --> FIG 2E-F ############
###################################################

fig_bottom <- plot_grid(plot_E, plot_F, nrow = 1, labels = c("E","F"), 
                                label_size = 20, rel_widths = c(0.45,0.55),
                                vjust = 0, label_y = 0.9)

#############################################
############ FINAL FIG --> FIG 2 ############
#############################################

Fig2_final <- plot_grid(fig_up, fig_mid, fig_bottom, nrow = 3, rel_heights = c(0.5,0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig2_complete.png"
ggsave(final_figure_path, Fig2_final, width = 170, height = 240, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig2_complete.pdf"
ggsave(final_figure_path, Fig2_final, width = 170, height = 240, units = "mm")

