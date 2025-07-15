library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/panel_A.RData")
# Change stuff
input_figureA$ind <- gsub("RandomGenes\nwith NMD\nfeatures","RandomGenes\nw/ NMD feat.",input_figureA$ind)
input_figureA$ind <- gsub("RandomGenes\nwithout NMD\nfeatures","RandomGenes\nw/o NMD feat.",input_figureA$ind)
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
        #geom_count()+#(alpha=0.1, size = 1) + 
        # geom_jitter(aes(fill = factor(factor(NMD_geneset, levels = c("RandomGenes","NMD")))), 
        #   position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
        #   alpha = 0.1, size = 1) + xlab("Percentiles") +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_cartesian(ylim = ylim) +
        geom_boxplot(width = boxplot_width, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(y = "ETG iNMDeff", x = "", title = "TCGA", fill = "Gene set") +
        ggplot_theme() + 
         scale_color_manual(values = c("Random" = "blue", "NMD" = "red"),
                       labels = c("Random", "NMD Genes")) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10, lineheight = 0.65),
              legend.margin = margin(t = -20, r = -0, b = 0, l = -0, unit = "pt"),
              plot.margin = unit(c(1, 0, 1, 0), "cm"),
              legend.text = element_text(size = 9),
              legend.title = element_blank(),
              legend.position = "bottom") +
        guides(fill = guide_legend(keyheight = unit(0.7, "cm"), keywidth = unit(0.7, "cm"))) +
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

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/panel_B.RData")
# Change stuff
input_figureB$ind <- gsub("RandomGenes\nwith NMD\nfeatures","RandomGenes\nw/ NMD feat.",input_figureB$ind)
input_figureB$ind <- gsub("RandomGenes\nwithout NMD\nfeatures","RandomGenes\nw/o NMD feat.",input_figureB$ind)
# combinations <- list(c("NMD All","NMD Consensus"), c("RandomGenes\nw/ NMD feat.","NMD Consensus"),c("RandomGenes\nw/o NMD feat.","NMD Consensus"))
combinations <- list(c("NMD All","NMD Consensus"),c("RandomGenes\nw/o NMD feat.","NMD Consensus"))
tmp_df <- input_figureB[!is.na(input_figureB$values),] # NAs
ylim <- c(-1,4)
annotate_hjust <- 0.5
pval_ylabel <- c(2.35,2.85,3.35)
boxplot_width <- 0.3
input_figureB[input_figureB$ind == "RandomGenes\nw/ NMD feat.","NMD_geneset"] <- "NMD"

# Plot
plot_B <- ggplot(data = input_figureB, aes(x = factor(ind), y = values, fill = factor(NMD_geneset, levels = c("RandomGenes","NMD")))) +
        #geom_count()+#(alpha=0.1, size = 1) + 
        # geom_jitter(aes(fill = factor(factor(NMD_geneset, levels = c("RandomGenes","NMD")))), 
        #   position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
        #   alpha = 0.1, size = 1) + xlab("Percentiles") +
        geom_violin(draw_quantiles = TRUE, na.rm = TRUE,lwd = 0.25) + coord_cartesian(ylim = ylim) +
        geom_boxplot(width = boxplot_width, color="black", alpha=0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(y = "ETG iNMDeff", x = "", title = "GTex", fill = "Gene set") +
        ggplot_theme() + 
         scale_color_manual(values = c("Random" = "blue", "NMD" = "red"),
                       labels = c("Random", "NMD Genes")) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10, lineheight = 0.65),
              legend.margin = margin(t = -20, r = -0, b = 0, l = -0, unit = "pt"),
              plot.margin = unit(c(1, 0, 1, 0), "cm"),
              legend.text = element_text(size = 9),
              legend.title = element_blank(),
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

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/panel_C.RData")
# Change stuff
ylim <- c(-1,2)
pval_ylabel <- c(0.5,1,1.5)
boxplot_width <- 0.12
annotate_hjust <- 0.5
tmp_df <- input_figureC[!is.na(input_figureC$values),] # NAs
input_figureC$NMD_geneset <- ifelse(input_figureC$NMD_geneset == "PTC NMD-triggering","NMD-triggering PTCs",input_figureC$NMD_geneset)
input_figureC$NMD_geneset <- ifelse(input_figureC$NMD_geneset == "PTC NMD-evading","NMD-evading PTCs",input_figureC$NMD_geneset)
input_figureC$ind <- as.character(input_figureC$ind)
input_figureC$ind <- ifelse(input_figureC$ind == "PTC NMD-triggering","NMD-triggering PTCs",input_figureC$ind)
input_figureC$ind <- factor(ifelse(input_figureC$ind == "PTC NMD-evading","NMD-evading PTCs",input_figureC$ind))
input_figureC$NMD_geneset <- factor(input_figureC$NMD_geneset, levels = c("Synonymous","NMD-triggering PTCs","NMD-evading PTCs"))
combinations <- combn(names(table(input_figureC$ind)), 2, simplify = FALSE)
combinations <- combinations[-2]
# Plot
plot_C <- ggplot(data = input_figureC, aes(x = ind, y = values, fill = NMD_geneset)) +
        geom_violin(draw_quantiles = TRUE, scale = "width", #bw = 0.1,
                        lwd = 0.25, na.rm = TRUE, width = 0.7) +
        geom_boxplot(width = boxplot_width, 
                        color="black", outlier.shape = NA, alpha=0.1) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        labs(y = "ASE iNMDeff", x = "", title = "GTex") +
        ggplot_theme() + coord_cartesian(ylim = ylim) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
                plot.margin = unit(c(1, 0, 1, 0), "cm"),
                legend.position = "none") +
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

####################################################
############ PANEL UP --> SUPP FIG 1A-C ############
####################################################

fig_up <- plot_grid(plot_A, plot_B, plot_C, nrow = 1, ncol = 3, labels = c("A", "B", "C"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.335,0.335,0.23))

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/panel_D.RData")
input_figureD$ind <- gsub("NMD Consensus", "NMD \n Consensus",input_figureD$ind)
#input_figureD$ind <- factor(input_figureD$ind, levels = c("NMD \n Consensus","RP9P", "GAS5", "SMG5"))
# selected_genes <- c("SMG1","SMG6","SMG7","UPF1","UPF2","UPF3A","UPF3B","SMG8","SMG9","CASC3","RBM8A","EIF4A3","MAGOH")
selected_genes <- c("RP9P","GAS5","SMG1","SMG5","SMG6","SMG7","SMG8","SMG9","UPF1","UPF2","UPF3B","UPF3A")
list_plots <- list() 

for (variable in c("NMD \n Consensus",selected_genes)) {

        plot_title <- ""
        face_text <- "italic"
        X_title_size <- 0
        input_figureD_filt <- input_figureD %>%
                                filter(ind == variable)
        size_title <- 11
        if (variable == "SMG5") {
                ylim <- c(4,10)
                # X_title_size <- 12
        } else if (variable == "SMG1") {
                ylim <- c(1,5)
        } else if (variable == "RP9P") {
                ylim <- c(1,4)
        } else if (variable == "GAS5") {
                ylim <- c(8,30)
        } else if (variable == "SMG6") {
                ylim <- c(2.5,7.5)
        } else if (variable == "SMG7") {
                ylim <- c(2.5,6.5)
        } else if (variable == "SMG8") {
                ylim <- c(1,4)
        } else if (variable == "SMG9") {
                ylim <- c(4,8)
        } else if (variable == "UPF1") {
                ylim <- c(4,8)
        } else if (variable == "UPF2") {
                ylim <- c(1.5,4.5)
        } else if (variable == "UPF3B") {
                ylim <- c(2,4)
        } else if (variable == "UPF3A") {
                ylim <- c(3,8)
        } else {
                ylim <- c(-4,4)
                plot_title <- "TCGA"
                face_text <- "plain"
                size_title <- 7
        }

        p <- ggplot(data = input_figureD_filt, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
                geom_smooth(data = subset(input_figureD_filt, ind %in% selected_genes), 
                                aes(y=values, group = factor(NMDeff_type, levels = c("Low","High"))), 
                                method = "gam", se = TRUE, linewidth = 1.5) +
                geom_point(alpha=0.1, size = 0.1) +
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

        if (!variable %in% "UPF3A") {
                p <- p + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = "none") +
                                guides(x = "none")     
        }

      list_plots[[length(list_plots) + 1]] <- p
}

plot_D <- ggarrange(list_plots[[1]],list_plots[[2]], list_plots[[3]], list_plots[[4]],
                list_plots[[5]],list_plots[[6]], list_plots[[7]],list_plots[[8]],list_plots[[9]],
                list_plots[[10]],list_plots[[11]],list_plots[[12]],list_plots[[13]],
                labels = c("","","",""), align = "v",
                font.label = list(size = 20), heights = c(0.17,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.42),
                ncol=1, nrow=length(list_plots), common.legend = FALSE)

# Add custom titles using grid

plot_D <- ggdraw(plot_D) +
  draw_label("iNMDeff", x = 0.04, y = 0.94, hjust = 0.5, vjust = 0.5, angle = 90, size = 12) +
  draw_label("Gene expression (TPM)", x = 0.04, y = 0.5, hjust = 0.5, vjust = 0.05, angle = 90, size = 12)

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot_D, width = 100, height = 150, units = "mm") 

################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig2/panel_E.RData")
input_figureE$ind <- gsub("NMD Consensus", "NMD \n Consensus",input_figureE$ind)
#input_figureD$ind <- factor(input_figureD$ind, levels = c("NMD \n Consensus","RP9P", "GAS5", "SMG5"))
list_plots <- list() 
selected_genes <- c("SMG1","SMG5","SMG6","SMG7","SMG8","SMG9","UPF1","UPF2","UPF3B","UPF3A")

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
                ylim <- c(7,18)
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
                geom_point(alpha=0.1, size = 0.1) +
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

        if (!variable %in% "UPF3A") {
                p <- p + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = "none") +
                                guides(x = "none")     
        }

      list_plots[[length(list_plots) + 1]] <- p
}

plot_E <- ggarrange(list_plots[[1]],list_plots[[2]], list_plots[[3]], list_plots[[4]],
                list_plots[[5]],list_plots[[6]], list_plots[[7]], list_plots[[8]],list_plots[[9]],
                list_plots[[10]],list_plots[[11]],
                labels = c("","","",""), align = "v",
                font.label = list(size = 20), heights = c(0.17,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.42),
                ncol=1, nrow=length(list_plots), common.legend = FALSE)

# Add custom titles using grid

plot_E <- ggdraw(plot_E) +
  draw_label("iNMDeff", x = 0.04, y = 0.94, hjust = 0.5, vjust = 0.05, angle = 90, size = 12) +
  draw_label("Gene expression (TPM)", x = 0.04, y = 0.5, hjust = 0.5, vjust = 0.05, angle = 90, size = 12)

################################
############ PLOT F ############
################################

input_figureF <- input_figureD
selected_genes <- c("CASC3","RBM8A","EIF4A3","MAGOH")
list_plots <- list() 

for (variable in c("NMD \n Consensus",selected_genes)) {

        plot_title <- ""
        face_text <- "italic"
        X_title_size <- 0
        input_figureF_filt <- input_figureF %>%
                                filter(ind == variable)
        size_title <- 11
        if (variable == "CASC3") {
                ylim <- c(3.5,8)
        } else if (variable == "RBM8A") {
                ylim <- c(8,16)
        } else if (variable == "EIF4A3") {
                ylim <- c(5,12)
        } else if (variable == "MAGOH") {
                ylim <- c(5,11)
        } else {
                ylim <- c(-4,4)
                plot_title <- "TCGA"
                face_text <- "plain"
                size_title <- 9
        }

        p <- ggplot(data = input_figureF_filt, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
                geom_smooth(data = subset(input_figureF_filt, ind %in% selected_genes), 
                                aes(y=values, group = factor(NMDeff_type, levels = c("Low","High"))), 
                                method = "gam", se = TRUE, linewidth = 1.5) +
                geom_point(alpha=0.15, size = 0.1) +
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

        if (!variable %in% "MAGOH") {
                p <- p + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = "none") +
                                guides(x = "none")     
        }

      list_plots[[length(list_plots) + 1]] <- p
}

plot_F <- ggarrange(list_plots[[1]],list_plots[[2]], list_plots[[3]], list_plots[[4]], list_plots[[5]],
                labels = c("","","",""), align = "v",
                font.label = list(size = 20), heights = c(0.17,0.22,0.22,0.22,0.30),
                ncol=1, nrow=length(list_plots), common.legend = FALSE)

# Add custom titles using grid
plot_F <- ggdraw(plot_F) +
  draw_label("iNMDeff", x = 0.04, y = 0.92, hjust = 0.5, vjust = 0.05, angle = 90, size = 12) +
  draw_label("Gene expression (TPM)", x = 0.04, y = 0.5, hjust = 0.5, vjust = 0.05, angle = 90, size = 12)

################################
############ PLOT G ############
################################

input_figureG <- input_figureE
list_plots <- list() 

for (variable in c("NMD \n Consensus",selected_genes)) {

        plot_title <- ""
        face_text <- "italic"
        X_title_size <- 0
        input_figureG_filt <- input_figureG %>%
                                filter(ind == variable)
        size_title <- 11
        if (variable == "CASC3") {
                ylim <- c(6,12)
        } else if (variable == "RBM8A") {
                ylim <- c(4,8)
        } else if (variable == "EIF4A3") {
                ylim <- c(5,12)
        } else if (variable == "MAGOH") {
                ylim <- c(3,7)
        } else {
                ylim <- c(-4,4)
                plot_title <- "GTex"
                face_text <- "plain"
                size_title <- 9
        }

        p <- ggplot(data = input_figureG_filt, aes(x = sample, y = values, color = factor(NMDeff_type, levels = c("Low","High")))) +
                geom_smooth(data = subset(input_figureG_filt, ind %in% selected_genes), 
                                aes(y=values, group = factor(NMDeff_type, levels = c("Low","High"))), 
                                method = "gam", se = TRUE, linewidth = 1.5) +
                geom_point(alpha=0.15, size = 0.1) +
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

        if (!variable %in% "MAGOH") {
                p <- p + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = "none") +
                                guides(x = "none")     
        }

      list_plots[[length(list_plots) + 1]] <- p
}

plot_G <- ggarrange(list_plots[[1]],list_plots[[2]], list_plots[[3]], list_plots[[4]], list_plots[[5]], 
                labels = c("","","",""), align = "v",
                font.label = list(size = 20), heights = c(0.17,0.22,0.22,0.22,0.30),
                ncol=1, nrow=length(list_plots), common.legend = FALSE)

# Add custom titles using grid
plot_G <- ggdraw(plot_G) +
  draw_label("iNMDeff", x = 0.04, y = 0.92, hjust = 0.5, vjust = 0.05, angle = 90, size = 12) +
  draw_label("Gene expression (TPM)", x = 0.04, y = 0.5, hjust = 0.5, vjust = 0.05, angle = 90, size = 12)

########################################################
############ PANEL BOTTOM --> SUPP FIG D-G ############
########################################################

fig_bottom <- plot_grid(plot_D, plot_E, plot_F, plot_G, nrow = 1, ncol = 4, labels = c("D","E","F","G"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5), label_y = 1.05)

# fig_bottom <- ggarrange(plot_D, plot_E, plot_F, plot_G, nrow = 1, ncol = 4, labels = c("D","E","F","G"), 
#                     font.label = list(size = 20),  align = "v", common.legend = TRUE)

###################################
############ FINAL FIG ############
###################################

SuppFig2_final <- plot_grid(fig_up, fig_bottom, nrow = 2, rel_heights = c(0.36,0.64))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig2_complete.png"
ggsave(final_figure_path, SuppFig2_final, width = 325, height = 325, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig2_complete.pdf"
ggsave(final_figure_path, SuppFig2_final, width = 325, height = 325, units = "mm")
