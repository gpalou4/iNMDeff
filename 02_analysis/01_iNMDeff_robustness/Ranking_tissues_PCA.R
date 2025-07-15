library("ggrepel")
library("ggplot2")
library("FactoMineR")
library("factoextra")
library("dplyr")
library("corrplot")
#conda activate /home/gpalou/anaconda3_envs/general

NMDeff_median_tissue_PCA <- function(randomization, sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove) {

  NMD_genesets <- c("endogenous_NMD_global_2_shared","endogenous_non_NMD_neg_control","ASE_stopgain_0.01","ASE_stopgain_0.2",
                    "ASE_synonymous_0.01","ASE_synonymous_0.2", "PTCs_stopgain_NMD_triggering","PTCs_synonymous")

  # Randomize samples before estimating NMDeff median by tissue
  if (randomization == "yes") {
    # Randomized TCGA
    random <- sample(1:nrow(sample_NMD_efficiencies_TCGA))
    cols <- which(colnames(sample_NMD_efficiencies_TCGA) != "cancer_type_strat")
    sample_NMD_efficiencies_TCGA[,cols] <- sample_NMD_efficiencies_TCGA[random,cols]
    # Randomized GTEx
    random <- sample(1:nrow(sample_NMD_efficiencies_GTEx))
    cols <- which(colnames(sample_NMD_efficiencies_GTEx) != "acronyms")
    sample_NMD_efficiencies_GTEx[,cols] <- sample_NMD_efficiencies_GTEx[random,cols]
    plot_char <- "_randomization"
  } else {
    plot_char <- ""
  }

  # 2) Calculate median NMDeff per tissue
  # 2.1) TCGA
  samples_NMDeff_median_TCGA <- data.frame( sample_NMD_efficiencies_TCGA %>%
                              group_by(cancer_type_strat) %>% 
                              summarise(across(all_of(NMD_genesets), ~ median(.x, na.rm = TRUE))) )
  colnames(samples_NMDeff_median_TCGA)[1] <- "tissues"
  samples_NMDeff_median_TCGA_sample_size <- data.frame( sample_NMD_efficiencies_TCGA %>%
                              group_by(cancer_type_strat) %>% 
                              summarise(across(all_of(NMD_genesets), ~ sum(!is.na(.x)))) )
  colnames(samples_NMDeff_median_TCGA_sample_size)[1] <- "tissues"
  # Scale
  samples_NMDeff_median_TCGA[,-1] <- scale(samples_NMDeff_median_TCGA[,-1])
  # 2.2) GTEx
  samples_NMDeff_median_GTEx <- data.frame( sample_NMD_efficiencies_GTEx %>%
                              group_by(acronyms) %>% 
                              summarise(across(all_of(NMD_genesets), ~ median(.x, na.rm = TRUE))) )
  colnames(samples_NMDeff_median_GTEx)[1] <- "tissues"
  samples_NMDeff_median_GTEx_sample_size <- data.frame( sample_NMD_efficiencies_GTEx %>%
                              group_by(acronyms) %>% 
                              summarise(across(all_of(NMD_genesets), ~ sum(!is.na(.x)))) )
  colnames(samples_NMDeff_median_GTEx_sample_size)[1] <- "tissues"
  # Scale
  samples_NMDeff_median_GTEx[,-1] <- scale(samples_NMDeff_median_GTEx[,-1])
  # 2.3) Merge
  samples_NMDeff_median_all_tissues <- rbind(samples_NMDeff_median_GTEx,samples_NMDeff_median_TCGA)
  samples_NMDeff_median_all_tissues_sample_size <- rbind(samples_NMDeff_median_GTEx_sample_size,samples_NMDeff_median_TCGA_sample_size)
  # Filter columns
  samples_NMDeff_median_all_tissues <- samples_NMDeff_median_all_tissues[,-grep("0.01",colnames(samples_NMDeff_median_all_tissues))]
  samples_NMDeff_median_all_tissues_sample_size <- samples_NMDeff_median_all_tissues_sample_size[,c("tissues","ASE_stopgain_0.2")]
  colnames(samples_NMDeff_median_all_tissues_sample_size)[2] <- "sample_size"
  # Remove outliers in PCA
  samples_NMDeff_median_all_tissues <- samples_NMDeff_median_all_tissues[!samples_NMDeff_median_all_tissues$tissues %in% tissues_remove,]
  samples_NMDeff_median_all_tissues_sample_size <- samples_NMDeff_median_all_tissues_sample_size[!samples_NMDeff_median_all_tissues_sample_size$tissues %in% tissues_remove,]
  # Color matched TCGA-GTEx tissues
  samples_NMDeff_median_all_tissues$matched_tissues <- NA
  for (j in 1:nrow(TCGA_GTEx_tissues)) {
    tissue <- TCGA_GTEx_tissues[j,"GTEx_tissues"]
    cancer <- TCGA_GTEx_tissues[j,"TCGA_cancers"]
    index <- grep(tissue,samples_NMDeff_median_all_tissues$tissues)
    samples_NMDeff_median_all_tissues[index,"matched_tissues"] <- paste0(TCGA_GTEx_tissues[j,"merge_tissues"],"_N")
    if (cancer == "ACC") {
      cancer <- "TCGA-ACC"
    }
    index <- grep(cancer,samples_NMDeff_median_all_tissues$tissues)
    samples_NMDeff_median_all_tissues[index,"matched_tissues"] <- paste0(TCGA_GTEx_tissues[j,"merge_tissues"],"_T")
  }
  samples_NMDeff_median_all_tissues$dataset <- "GTEx"
  samples_NMDeff_median_all_tissues[grep("TCGA",samples_NMDeff_median_all_tissues$tissues),"dataset"] <- "TCGA"
  samples_NMDeff_median_all_tissues[samples_NMDeff_median_all_tissues$tissues %in% "TCGA-UCEC_MSS_SBS10ab","matched_tissues"] <- NA
  colnames(samples_NMDeff_median_all_tissues) <- c("tissues","Endogenous","Endogenous_neg_control","ASE","ASE_neg_control","PTC","PTC_neg_control","matched_tissues","dataset")
  # NA
  samples_NMDeff_median_all_tissues <- samples_NMDeff_median_all_tissues[!is.na(samples_NMDeff_median_all_tissues$tissues),]
  samples_NMDeff_median_all_tissues_sample_size <- samples_NMDeff_median_all_tissues_sample_size[!is.na(samples_NMDeff_median_all_tissues_sample_size$tissues),]
  # Color Brain
  samples_NMDeff_median_all_tissues$brain <- "other"
  samples_NMDeff_median_all_tissues[grep("Brain",samples_NMDeff_median_all_tissues$matched_tissues),"brain"] <- "brain"
  # Sample Size
  samples_NMDeff_median_all_tissues <- merge(samples_NMDeff_median_all_tissues,samples_NMDeff_median_all_tissues_sample_size, by = "tissues", all.x = TRUE)
  # 3) PCA
  #pca_iris_rotated <- psych::principal(na.omit(samples_NMDeff_median_all_tissues[,-1]), rotate="varimax", nfactors=5, scores=TRUE)

  for (PTC in c("yes","no")) {

    if (PTC == "yes") {
      samples_NMDeff_median_all_tissues_filt <- samples_NMDeff_median_all_tissues
    } else if (PTC == "no") {
      samples_NMDeff_median_all_tissues_filt <- samples_NMDeff_median_all_tissues[,-grep("PTC",colnames(samples_NMDeff_median_all_tissues))]
    }
    #colnames(samples_NMDeff_median_all_tissues) <- c("tissues","Endogenous_neg_control","Endogenous","ASE","ASE_neg_control","ASE_evading","PTC","PTC_neg_control","PTC_evading")
    #colnames(samples_NMDeff_median_all_tissues) <- c("tissues","Endogenous_neg_control","Endogenous","ASE","ASE_neg_control")

    for (type in c("TCGA","GTEx","TCGA_GTEx")) {

      # Filter
      if (type != "TCGA_GTEx") {
        samples_NMDeff_median_all_tissues_filt2 <- samples_NMDeff_median_all_tissues_filt[samples_NMDeff_median_all_tissues_filt$dataset == type,]
        fill <- "brain"
        labels <- "tissues"
      } else {
        samples_NMDeff_median_all_tissues_filt2 <- samples_NMDeff_median_all_tissues_filt
        fill <- "brain"
        labels <- "matched_tissues"
      }

      # if (type == "TCGA_GTEx") {

      # for (i in 1:nrow(TCGA_GTEx_tissues)) {

      # tissue <- TCGA_GTEx_tissues[i,"merge_tissues"]
      # index <- grep(tissue,samples_NMDeff_median_all_tissues_filt2$matched_tissues)
      # samples_NMDeff_median_all_tissues_filt2$color <- NA
      # samples_NMDeff_median_all_tissues_filt2[index,"color"] <- samples_NMDeff_median_all_tissues_filt2[index,"tissues"]
      

      ###################### TRYING INDIVIDUALS??

      # cols <- c("endogenous_NMD_global_2_shared","endogenous_non_NMD_neg_control","ASE_stopgain_0.2","ASE_synonymous_0.2")
      # sample_NMD_efficiencies_TCGA_filt <- na.omit(sample_NMD_efficiencies_TCGA[,cols])
      # sample_NMD_efficiencies_TCGA_filt <- merge(sample_NMD_efficiencies_TCGA_filt,sample_NMD_efficiencies_TCGA, by = cols, all.x = TRUE)
      # res_pca <- PCA(sample_NMD_efficiencies_TCGA_filt[,cols], 
      #               scale.unit = TRUE, ncp = 1000, graph = FALSE)

      # # Color Brain
      # sample_NMD_efficiencies_TCGA_filt$brain <- "other"
      # sample_NMD_efficiencies_TCGA_filt[grep("LGG|GBM",sample_NMD_efficiencies_TCGA_filt$cancer_type_strat),"brain"] <- "brain"
      
      # fill <- "cancer_type_strat"
      # PCA_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/"
      # png(paste0(PCA_output_path,"/tissue_ranking/PCA/test.png"), width = 2500, height = 2000, res = 300)
      # p <- fviz_pca_biplot(res_pca, repel = FALSE, geom.ind = "point",
      #                       col.var = "#2E9FDF",
      #                       axes = c(1,2),
      #                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
      #                       alpha.var = 0.3, 
      #                       alpha.ind = 0.5,
      #                       pointshape = 21, 
      #                       #addEllipses = TRUE,
      #                       #size ="sample_size",
      #                       #pointsize = sample_NMD_efficiencies_TCGA_filt$cancer_type_strat,
      #                       mean.point = TRUE,
      #                       fill.ind = factor(sample_NMD_efficiencies_TCGA_filt[,fill]),
      #                       legend.title = list(fill = "NMDeff")) +
      #   #geom_label_repel(aes(label=as.factor(samples_NMDeff_median_all_tissues$cancer_type), color = as.factor(samples_NMDeff_median_all_tissues$tissues)), size=3, nudge_y=0.05, max.overlaps = 33) +
      #   #geom_label_repel(aes(label=as.factor(sample_NMD_efficiencies_TCGA_filt[,"cancer_type_strat"])), size=2.5, nudge_y=0.07, max.overlaps = nrow(sample_NMD_efficiencies_TCGA_filt)) +
      #   guides(colour=FALSE) + theme_classic() + ggtitle(type) +
      #   theme(plot.title = element_text(hjust = 0.5, size = 16),
      #     axis.title.x = element_text(color="black", size=16, face="bold"),
      #     axis.title.y = element_text(color="black", size=16, face="bold"),
      #     legend.position='top') + scale_size_continuous(breaks=seq(1, 650, 50), limits = c(1, 553))
      # print(p)
      # dev.off()

      ##########################

      res_pca <- PCA(na.omit(samples_NMDeff_median_all_tissues_filt2[,-which(colnames(samples_NMDeff_median_all_tissues_filt2) %in% c("tissues","matched_tissues","dataset","brain","sample_size"))]), scale.unit = FALSE, ncp = 1000, graph = FALSE)
      PCA_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/"
      png(paste0(PCA_output_path,"/tissue_ranking/PCA/PCA_all_methods_",type,plot_char,"_PTC_",PTC,".png"), width = 2500, height = 2000, res = 300)
      #png(paste0(PCA_output_path,"/tissue_ranking/PCA/tissues/PCA_all_methods_",type,"_",tissue,".png"), width = 2500, height = 2000, res = 300)
      p <- fviz_pca_biplot(res_pca, repel = FALSE, geom.ind = "point",
                            col.var = "#2E9FDF",
                            axes = c(1,2),
                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                            alpha.var = 0.5, 
                            pointshape = 21, 
                            #size ="sample_size",
                            pointsize = samples_NMDeff_median_all_tissues_filt2$sample_size,
                            mean.point = TRUE,
                            fill.ind = factor(samples_NMDeff_median_all_tissues_filt2[,fill]),
                            legend.title = list(fill = "NMDeff")) +
        #geom_label_repel(aes(label=as.factor(samples_NMDeff_median_all_tissues$cancer_type), color = as.factor(samples_NMDeff_median_all_tissues$tissues)), size=3, nudge_y=0.05, max.overlaps = 33) +
        geom_label_repel(aes(label=as.factor(samples_NMDeff_median_all_tissues_filt2[,labels])), size=2.5, nudge_y=0.07, max.overlaps = nrow(samples_NMDeff_median_all_tissues_filt2)) +
        guides(colour=FALSE) + theme_classic() + ggtitle(type) +
        theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.position='top') + scale_size_continuous(breaks=seq(1, 650, 50), limits = c(1, 553))
      print(p)
      dev.off()
      # }
      # }
    }
  }
} # END of function

NMDeff_samples_PCA <- function(test, sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove) {

  NMD_genesets <- c("endogenous_NMD_global_2_shared","endogenous_non_NMD_neg_control","ASE_stopgain_0.01","ASE_stopgain_0.2",
                    "ASE_synonymous_0.01","ASE_synonymous_0.2", "PTCs_stopgain_NMD_triggering","PTCs_synonymous")

  # Randomize samples before estimating samples-NMDeff
  if (test == "randomization") {
    # Randomized TCGA
    random <- sample(1:nrow(sample_NMD_efficiencies_TCGA))
    cols <- which(colnames(sample_NMD_efficiencies_TCGA) != "cancer_type_strat")
    sample_NMD_efficiencies_TCGA[,cols] <- sample_NMD_efficiencies_TCGA[random,cols]
    # Randomized GTEx
    random <- sample(1:nrow(sample_NMD_efficiencies_GTEx))
    cols <- which(colnames(sample_NMD_efficiencies_GTEx) != "acronyms")
    sample_NMD_efficiencies_GTEx[,cols] <- sample_NMD_efficiencies_GTEx[random,cols]
    plot_char <- "_randomization"
  } else if (test == "sample_size") {
    plot_char <- "_sampleSize"
    # Reduce Sample Size for 10 tissues including Brain
    # TCGA
    #random_tissues <- unique(c(sample(unique(sample_NMD_efficiencies_TCGA$cancer_type))[1:10],"TCGA-LGG","TCGA-GBM"))
    sample_NMD_efficiencies_TCGA <- data.frame(sample_NMD_efficiencies_TCGA %>% group_by(cancer_type) %>% slice_sample(n=30))
    # Randomized GTEx
    #random_tissues <- unique(c(sample(unique(sample_NMD_efficiencies_TCGA$cancer_type))[1:10],"TCGA-LGG","TCGA-GBM"))
    sample_NMD_efficiencies_GTEx <- sample_NMD_efficiencies_GTEx[!sample_NMD_efficiencies_GTEx$acronyms %in% c("BLDDER","CVXECT","FLLPNT","KDNMDL","CVSEND"),]
    sample_NMD_efficiencies_GTEx <- data.frame(sample_NMD_efficiencies_GTEx %>% group_by(acronyms) %>% slice_sample(n=30))
  } else {
    plot_char <- ""
  }
  
  # TCGA scale
  sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA[,c("sample","cancer_type_strat",NMD_genesets)]
  # Scale
  sample_NMD_efficiencies_TCGA_filt[,NMD_genesets] <- scale(sample_NMD_efficiencies_TCGA[,NMD_genesets])
  colnames(sample_NMD_efficiencies_TCGA_filt)[2] <- "tissue"
  # GTEx
  sample_NMD_efficiencies_GTEx_filt <- sample_NMD_efficiencies_GTEx[,c("sample_full_barcode","acronyms",NMD_genesets)]
  # Scale
  sample_NMD_efficiencies_GTEx_filt[,NMD_genesets] <- scale(sample_NMD_efficiencies_GTEx_filt[,NMD_genesets])
  colnames(sample_NMD_efficiencies_GTEx_filt)[1:2] <- c("sample","tissue")
  # Merge
  samples_NMDeff_all <- rbind(sample_NMD_efficiencies_GTEx_filt,sample_NMD_efficiencies_TCGA_filt)
  # Filter columns
  samples_NMDeff_all <- samples_NMDeff_all[,-grep("0.01",colnames(samples_NMDeff_all))]
  # Add columns
  samples_NMDeff_all$dataset <- NA
  samples_NMDeff_all[grep("GTEX",samples_NMDeff_all$sample),"dataset"] <- "GTEx"
  samples_NMDeff_all[grep("TCGA",samples_NMDeff_all$sample),"dataset"] <- "TCGA"

  # 3) PCA

  for (PTC in c("yes","no")) {

    if (PTC == "yes") {
      samples_NMDeff_all_filt <- samples_NMDeff_all
    } else if (PTC == "no") {
      samples_NMDeff_all_filt <- samples_NMDeff_all[,-grep("PTC",colnames(samples_NMDeff_all))]
    }

    NMDeff_PCs_list <- list()
    for (type in c("TCGA","GTEx","TCGA_GTEx")) {

      # Filter
      if (type != "TCGA_GTEx") {
        samples_NMDeff_all_filt2 <- samples_NMDeff_all_filt[samples_NMDeff_all_filt$dataset == type,]
        fill <- "brain"
        labels <- "tissue"
      } else {
        samples_NMDeff_all_filt2 <- samples_NMDeff_all_filt
        fill <- "brain"
        labels <- "matched_tissues"
      }

      # 1) PCA of samples-NMDeff
      
      num_cols <- unlist(lapply(samples_NMDeff_all_filt2, is.numeric))
      cols <- colnames(samples_NMDeff_all_filt2)[num_cols]
      samples_NMDeff_all_filt2 <- na.omit(samples_NMDeff_all_filt2)
      #samples_NMDeff_all_filt2 <- merge(samples_NMDeff_all_filt2,sample_NMD_efficiencies_TCGA, by = c("sample"), all.x = TRUE)
      res_pca <- PCA(samples_NMDeff_all_filt2[,cols], 
                    scale.unit = FALSE, ncp = 1000, graph = FALSE)
      # Color Brain
      samples_NMDeff_all_filt2$brain <- "other"
      samples_NMDeff_all_filt2[grep("LGG|GBM|BRN",samples_NMDeff_all_filt2$tissue),"brain"] <- "brain"
      
      # PCA
      PCA_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/tissue_ranking/PCA/samples/"
      png(paste0(PCA_output_path,"/PCA_samples_",type,plot_char,"_PTC_",PTC,".png"), width = 2500, height = 2000, res = 300)
      p <- fviz_pca_biplot(res_pca, repel = FALSE, geom.ind = "point",
                            col.var = "#2E9FDF",
                            axes = c(1,2),
                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                            alpha.var = 1, 
                            alpha.ind = 0.7,
                            pointshape = 21, 
                            addEllipses = TRUE,
                            #size ="sample_size",
                            #pointsize = sample_NMD_efficiencies_TCGA_filt$cancer_type_strat,
                            mean.point = TRUE,
                            fill.ind = factor(data.frame(samples_NMDeff_all_filt2[,fill])[,1]),
                            legend.title = list(fill = "NMDeff")) +
        guides(colour=FALSE) + theme_classic() + ggtitle(type) +
        #geom_label_repel(aes(label=as.factor(samples_NMDeff_all_filt2[,labels])), size=0.5, nudge_y=0.07, max.overlaps = nrow(samples_NMDeff_all_filt2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(color="black", size=16, face="bold"),
          axis.title.y = element_text(color="black", size=16, face="bold"),
          legend.position='top') + scale_size_continuous(breaks=seq(1, 650, 50), limits = c(1, 553)) #+ coord_cartesian(ylim=c(-4,4))
      print(p)
      dev.off()

      # 2) PCA of median samples-NMDeff by tissue
      NMDeff_PCs <- data.frame(res_pca$ind$coord)
      NMDeff_PCs$sample <- samples_NMDeff_all_filt2$sample
      NMDeff_PCs$tissue <- samples_NMDeff_all_filt2$tissue
      NMDeff_PCs$brain <- samples_NMDeff_all_filt2$brain
      # Median by tissue
      NMDeff_PCs_median_tissues <- data.frame( NMDeff_PCs %>%
                                  group_by(tissue, brain) %>% 
                                  summarise(across(all_of(c("Dim.1","Dim.2","Dim.3","Dim.4")), ~ median(.x, na.rm = TRUE))) )
      # Sample Size
      samples_NMDeff_median_GTEx_sample_size <- data.frame( sample_NMD_efficiencies_GTEx %>%
                                  group_by(acronyms) %>% 
                                  summarise(across(all_of("ASE_stopgain_0.2"), ~ sum(!is.na(.x)))) )
      colnames(samples_NMDeff_median_GTEx_sample_size)[1] <- "tissue"
      samples_NMDeff_median_TCGA_sample_size <- data.frame( sample_NMD_efficiencies_TCGA %>%
                                  group_by(cancer_type_strat) %>% 
                                  summarise(across(all_of("ASE_stopgain_0.2"), ~ sum(!is.na(.x)))) )
      colnames(samples_NMDeff_median_TCGA_sample_size)[1] <- "tissue"
      sample_size <- na.omit(rbind(samples_NMDeff_median_TCGA_sample_size,samples_NMDeff_median_GTEx_sample_size))
      colnames(sample_size)[2] <- "sample_size"
      NMDeff_PCs_median_tissues <- merge(NMDeff_PCs_median_tissues,sample_size, by = "tissue", all.x = TRUE)
      # Filter PCA outliers
      NMDeff_PCs_median_tissues <- NMDeff_PCs_median_tissues[!NMDeff_PCs_median_tissues$tissue %in% tissues_remove,]
      # Color matched TCGA-GTEx tissues
      NMDeff_PCs_median_tissues$matched_tissues <- NA
      for (j in 1:nrow(TCGA_GTEx_tissues)) {
        tissue <- TCGA_GTEx_tissues[j,"GTEx_tissues"]
        cancer <- TCGA_GTEx_tissues[j,"TCGA_cancers"]
        index <- grep(tissue,NMDeff_PCs_median_tissues$tissue)
        NMDeff_PCs_median_tissues[index,"matched_tissues"] <- paste0(TCGA_GTEx_tissues[j,"merge_tissues"],"_N")
        if (cancer == "ACC") {
          cancer <- "TCGA-ACC"
        }
        index <- grep(cancer,NMDeff_PCs_median_tissues$tissue)
        NMDeff_PCs_median_tissues[index,"matched_tissues"] <- paste0(TCGA_GTEx_tissues[j,"merge_tissues"],"_T")
      }
      NMDeff_PCs_median_tissues$dataset <- "GTEx"
      NMDeff_PCs_median_tissues[grep("TCGA",NMDeff_PCs_median_tissues$tissue),"dataset"] <- "TCGA"
      NMDeff_PCs_median_tissues[NMDeff_PCs_median_tissues$tissue %in% "TCGA-UCEC_MSS_SBS10ab","matched_tissues"] <- NA
      # PCA
      PCA_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/tissue_ranking/PCA/samples/"
      png(paste0(PCA_output_path,"PCA_tissues_",type,plot_char,"_PTC_",PTC,".png"), width = 2500, height = 2000, res = 300)
      p <- ggplot(data = NMDeff_PCs_median_tissues, aes(x=Dim.1, y = Dim.2, fill = eval(parse(text = fill)))) +
          geom_point(mapping=aes(fill=eval(parse(text = fill)),size = sample_size), shape = 21) + theme_classic() +
          xlab(paste("Dim.1: ",round(res_pca$eig[1,2],2),"%", sep = "")) +
          ylab(paste("Dim.2: ",round(res_pca$eig[2,2],2),"%", sep = "")) +
          ggtitle(paste0(type)) +
          geom_hline(yintercept=0,linetype=2) + 
          geom_vline(xintercept=0,linetype=2) + 
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(color="black", size=18),
                axis.text.y = element_text(color="black", size=18),
                legend.position = "none", axis.text=element_text(size=18)) +
          geom_label_repel(aes(label=as.factor(NMDeff_PCs_median_tissues[,labels])), size=2.5, nudge_y=0.07, max.overlaps = nrow(NMDeff_PCs_median_tissues)) +
          scale_size_continuous(breaks=seq(1, 650, 50), limits = c(min(sample_size$sample_size), max(sample_size$sample_size)))
      print(p)
      dev.off()
      # 3) Save NMDeff-PCs for individuals
      if (PTC == "no") {
        NMDeff_PCs <- data.frame(res_pca$ind$coord)
        NMDeff_PCs$sample <- samples_NMDeff_all_filt2$sample
        NMDeff_PCs$tissue <- samples_NMDeff_all_filt2$tissue
        NMDeff_PCs_list[[type]] <- NMDeff_PCs
      }
    }
  }
  if (!test %in% c("randomization","sample_size")) {
    return(NMDeff_PCs_list)
  }
} # END of function


# 1) Data
# 1.1) sample NMD efficiencies TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert some columns to factors
factor_cols <- c("cancer_type","cancer_type_MSI","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor) 
cancers <- unique(sample_NMD_efficiencies_TCGA$cancer_type_strat)
# Filters
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.2 < 3),c("ASE_stopgain_0.2","ASE_stopgain_NMD_evading_0.2")] <- NA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.01 < 3),c("ASE_stopgain_0.01","ASE_stopgain_NMD_evading_0.01")] <- NA
#sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$PTC_num_PTCs < 5),c("PTCs_stopgain_NMD_triggering","PTCs_stopgain_NMD_evading")] <- NA

# 1.3) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert some columns to factors
factor_cols <- c("tissue","sample")
sample_NMD_efficiencies_GTEx[factor_cols] <- lapply(sample_NMD_efficiencies_GTEx[factor_cols], factor) 
GTEx_tissues <- unique(sample_NMD_efficiencies_GTEx$acronyms)
# Filters for GTEx
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.2 < 3),c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.01 < 3),c("ASE_stopgain_0.01")] <- NA

# 1.4) GTEx - TCGA matched tissues
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_GTEx_match.txt")
TCGA_GTEx_tissues <- read.table(file = output_path, header = TRUE, sep = "\t")
colnames(TCGA_GTEx_tissues)[1] <- "merge_tissues"
# # iCluster metadata
# TCGA_iCluster <- as.data.frame(read_excel(path = "/g/strcombio/fsupek_home/gpalou/data/TCGA_stratification/TCGA_iCluster_2019.xlsx",
#                                                   sheet = "Table S6 - iCluster", skip = 1))
# sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA,TCGA_iCluster, by.x = "sample", by.y="Sample ID", all.x = TRUE)

# 2) PCA
tissues_remove <- c("KDNMDL","TCGA-DLBC","CVSEND","BLDDER","FLLPNT")
NMDeff_median_tissue_PCA(randomization = "no", sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove = tissues_remove)
NMDeff_median_tissue_PCA(randomization = "yes", sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove = tissues_remove)

tissues_remove <- c("KDNMDL","TCGA-DLBC","CVSEND","BLDDER","FLLPNT")
NMDeff_PCs <- NMDeff_samples_PCA(test = "no", sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove = tissues_remove)
NMDeff_samples_PCA(test = "randomization", sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove = tissues_remove)
NMDeff_samples_PCA(test = "sample_size", sample_NMD_efficiencies_TCGA, sample_NMD_efficiencies_GTEx, tissues_remove = tissues_remove)
 
# 3) Add PCs
# 3.1) To TCGA
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMDeff_PCs_TCGA <- NMDeff_PCs$TCGA[,c("sample","tissue","Dim.1")]
colnames(NMDeff_PCs_TCGA)[3] <- "NMDeff_PCA"
sample_NMD_efficiencies_TCGA_filt <- merge(sample_NMD_efficiencies_TCGA,NMDeff_PCs_TCGA[,c("sample","NMDeff_PCA")], by = "sample", all.x = TRUE)
write.table(sample_NMD_efficiencies_TCGA_filt, file = sample_NMD_efficiencies_TCGA_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 3.2) To GTEx
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMDeff_PCs_GTEx <- NMDeff_PCs$GTEx[,c("sample","tissue","Dim.1")]
colnames(NMDeff_PCs_GTEx)[3] <- "NMDeff_PCA"
sample_NMD_efficiencies_GTEx_filt <- merge(sample_NMD_efficiencies_GTEx,NMDeff_PCs_GTEx[,c("sample","NMDeff_PCA")], by.x = "sample_full_barcode", by.y = "sample", all.x = TRUE)
write.table(sample_NMD_efficiencies_GTEx_filt, file = sample_NMD_efficiencies_GTEx_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
