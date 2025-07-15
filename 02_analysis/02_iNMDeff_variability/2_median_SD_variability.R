# 1) Obtain median SD variability results for each NMD method and Database

databases <- c("GTEx","TCGA")
NMD_methods <- c("ASE","endogenous")
tissues_percentiles_SD_of_medians <- c()

for (database in databases) {
    if (database == "TCGA") {
        path_str <- "/cancers/pancancer/"
    } else if (database == "GTEx") {
        path_str <- "/GTEx/pantissue/"
    }
    for (NMD_method in NMD_methods) {
        for (random_iteration in 1:3) {
            # file_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project",path_str,NMD_method,"/boxplots/NMD_genesets/",database,"_",NMD_method,"_randomization_median_SD.txt")
            file_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project",path_str,NMD_method,"/boxplots/NMD_genesets/",database,"_",NMD_method,"_randomization_median_SD_random_iteration_",random_iteration,".RData")
            percentiles_SD_of_medians <- readRDS(file_path)
            percentiles_SD_of_medians_df <- percentiles_SD_of_medians$sd_of_medians_diff
            # percentiles_SD <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
            percentiles_SD_of_medians_df$NMD_geneset <- gsub(paste0(NMD_method,"_"),"",rownames(percentiles_SD_of_medians_df))
            percentiles_SD_of_medians_df$database <- database
            percentiles_SD_of_medians_df$NMD_method <- NMD_method
            percentiles_SD_of_medians_df$random_iteration <- random_iteration
            if (length(tissues_percentiles_SD_of_medians) == 0) {
                tissues_percentiles_SD_of_medians <- percentiles_SD_of_medians_df
            } else {
                tissues_percentiles_SD_of_medians <- rbind(tissues_percentiles_SD_of_medians,percentiles_SD_of_medians_df)
            }
        }
    }
}
rownames(tissues_percentiles_SD_of_medians) <- NULL

write.table(tissues_percentiles_SD_of_medians, file = "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/percentiles_SD_of_medians_randomized_test.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# 2) Barplot

library("ggplot2")
library("dplyr")
library("viridis")
library("stringr")

# Some fix
tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("_0.2","",tissues_percentiles_SD_of_medians$NMD_geneset)
tissues_percentiles_SD_of_medians$controls <- "NMD set"
tissues_percentiles_SD_of_medians[tissues_percentiles_SD_of_medians$NMD_geneset %in% c("Synonymous 0.2","PTC NMD-evading 0.2","RandomGenes without NMD features","RandomGenes without NMD features"),"controls"] <- "Control set"
tissues_percentiles_SD_of_medians$NMD_geneset <- factor(tissues_percentiles_SD_of_medians$NMD_geneset, levels = unique(tissues_percentiles_SD_of_medians$NMD_geneset))
# Naming
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("PTC_NMD_triggering","PTC NMD-triggering",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("PTC_NMD_evading","PTC NMD-evading",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("synonymous","Synonymous",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("NMD_all","NMD All",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("NMD_Consensus","NMD Consensus",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("non_NMD_neg_control_with_NMD_features","RandomGenes with NMD features",tissues_percentiles_SD_of_medians$NMD_geneset)
# tissues_percentiles_SD_of_medians$NMD_geneset <- gsub("non_NMD_neg_control","RandomGenes without NMD features",tissues_percentiles_SD_of_medians$NMD_geneset)
tissues_percentiles_SD_of_medians$NMD_method <- gsub("endogenous","ETG",tissues_percentiles_SD_of_medians$NMD_method)
# Filter
tissues_percentiles_SD_of_medians_filt <- tissues_percentiles_SD_of_medians[tissues_percentiles_SD_of_medians$NMD_geneset %in% 
                                c("PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2","NMD All","NMD Consensus",
                                    "RandomGenes with NMD features","RandomGenes without NMD features"),]
# Order genesets by median SD
tissues_percentiles_SD_of_medians_filt$NMD_geneset <- str_wrap(tissues_percentiles_SD_of_medians_filt$NMD_geneset, width = 15) 
tissues_percentiles_SD_of_medians_filt$NMD_geneset <- factor(tissues_percentiles_SD_of_medians_filt$NMD_geneset , levels = unique(tissues_percentiles_SD_of_medians_filt$NMD_geneset[order(tissues_percentiles_SD_of_medians_filt$median_SD_of_medians_diff, decreasing = FALSE)]))

# Average of randomizations
tissues_percentiles_SD_of_medians_means <- tissues_percentiles_SD_of_medians_filt %>%
    group_by(NMD_geneset, database, NMD_method,controls) %>%
    summarise(
        median_SD_of_medians_diff_mean = mean(median_SD_of_medians_diff, na.rm = TRUE),
        perc5_SD_of_medians_diff_mean = mean(perc5_SD_of_medians_diff, na.rm = TRUE),
        perc95_SD_of_medians_diff_mean = mean(perc95_SD_of_medians_diff, na.rm = TRUE),
        p_value_mean = mean(p_value, na.rm = TRUE)
    )

# Plot
p <- ggplot(data = tissues_percentiles_SD_of_medians_means, aes(x = NMD_geneset, y = median_SD_of_medians_diff_mean, color = factor(NMD_method), fill = factor(controls))) +
    geom_bar(stat = "identity", linewidth = 0.7) + 
    # Add error bars for confidence intervals
    geom_errorbar(aes(ymin = perc5_SD_of_medians_diff_mean, ymax = perc95_SD_of_medians_diff_mean), 
                    width = 0.2) +  # Adjust width as needed
    guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
    facet_grid(. ~ database) + coord_flip() + #coord_cartesian(ylim = c(-5,5)) +
    ylab("Variability Deviation") + ggtitle("") + xlab("") +
    theme_bw() + scale_fill_brewer(palette = "Pastel1") + scale_color_viridis(discrete=TRUE, option="inferno", direction = -1) +
    theme(plot.title = element_text(hjust = 0.5, size = 35),
        axis.title.x = element_text(color="black", size=33),
        axis.title.y = element_text(color="black", size=35),
        axis.text.x = element_text(color="black", size=25),# angle = 90, hjust = 1),
        axis.text.y = element_text(color="black", size=18),
        strip.text = element_text(color="black", size = 30),
        legend.position='top', legend.text = element_text(size = 25)) +
    geom_text(aes(label = ifelse(p_value_mean < 0.05, "*", "")), 
                         position = position_dodge(width = 1), size = 15, hjust = -0.1, color = "black", vjust = 0.75)
# Save
plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/variability/percentiles_SD_of_medians_randomization.png")
png(plot_path, width = 3700, height = 2450, res = 300)
ggsave(gsub(".png",".pdf",plot_path), p, width = 4500, height = 3500, units = "px")
print(p)
dev.off()

# Save
write.table(tissues_percentiles_SD_of_medians_means, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(tissues_percentiles_SD_of_medians_means, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig2/Fig2B.RData")

# T-test
x <- tissues_percentiles_SD_of_medians_means[tissues_percentiles_SD_of_medians_means$controls == "NMD set","median_SD_of_medians_diff_mean"]
y <- tissues_percentiles_SD_of_medians_means[tissues_percentiles_SD_of_medians_means$controls != "NMD set","median_SD_of_medians_diff_mean"]

t.test(x, y, alternative = "two.sided", var.equal = FALSE)


## Supplementary
tissues_percentiles_SD_of_medians_filt <- tissues_percentiles_SD_of_medians

# Order genesets by median SD
tissues_percentiles_SD_of_medians_filt$NMD_geneset <- str_wrap(tissues_percentiles_SD_of_medians_filt$NMD_geneset, width = 15) 
tissues_percentiles_SD_of_medians_filt$NMD_geneset <- factor(tissues_percentiles_SD_of_medians_filt$NMD_geneset , levels = unique(tissues_percentiles_SD_of_medians_filt$NMD_geneset[order(tissues_percentiles_SD_of_medians_filt$median_SD_of_medians_diff, decreasing = FALSE)]))

# Average of randomizations
tissues_percentiles_SD_of_medians_means <- tissues_percentiles_SD_of_medians_filt %>%
    group_by(NMD_geneset, database, NMD_method,controls) %>%
    summarise(
        median_SD_of_medians_diff_mean = mean(median_SD_of_medians_diff, na.rm = TRUE),
        perc5_SD_of_medians_diff_mean = mean(perc5_SD_of_medians_diff, na.rm = TRUE),
        perc95_SD_of_medians_diff_mean = mean(perc95_SD_of_medians_diff, na.rm = TRUE),
        p_value_mean = mean(p_value, na.rm = TRUE)
    )

# Save
write.table(tissues_percentiles_SD_of_medians_means, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(tissues_percentiles_SD_of_medians_means, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4C.RData")
