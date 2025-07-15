library("ggplot2")

NMD_eff_correlations_old <- function(TCGA.cancer, num_PTC) {

    correlations.res <- expand.grid(VAFs,methods)
    colnames(correlations.res) <- c("VAF","methods")
    correlations.res[,c("stopgain_PTCs_VS_stopgain_ASE","num_samples")] <- NA

    # 1) PTC NMDeff

    # Current PTC 
    PTCs <- read.table(file = paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/endogenous/",TCGA.cancer,"_PTCs_NB_coeff_gene_id-pca_subtypes-somatic_transcripts_filt-0.01-4.txt"),
                                header = TRUE, sep = "\t", row.names = 1)
    PTCs <- PTCs[,c("stopgain","synonymous")]
    colnames(PTCs) <- c("stopgain_PTCs","synonymous_PTCs")

    # 2) ASE NMDeff

    # Correlation between ASE and PTCs for each different ASE VAF and method parameter

    for (VAF in VAFs) {

        for (method in methods) {

        row <- which(correlations.res$VAF%in%VAF & correlations.res$method%in%method)

        ASE.file <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/ASE/",TCGA.cancer,"_ASE_NB_coeff",num_PTC,"_",VAF,"-",method,".txt")
         
        # Current ASE
        ASE <- read.table(file = ASE.file,header = TRUE, sep = "\t", row.names = 1)
        correlations.res[row,"ASE_diff_stopgain_vs_NMD_evading"] <- abs(median(ASE$stopgain,na.rm = TRUE)-median(ASE$stopgain.NMD.evading,na.rm = TRUE))
        correlations.res[row,"num_ASE"] <- median(ASE$num.PTCs,na.rm=TRUE)
        ASE <- ASE[,c("stopgain","synonymous")]
        colnames(ASE) <- c("stopgain_ASE","synonymous_ASE")
        correlations.res[row,"num_samples"] <- sum(!is.na(ASE$stopgain_ASE))  

        # Merge
        ASE.PTCs <- merge(ASE,PTCs, by.x="row.names",by.y="row.names")
        rownames(ASE.PTCs) <- ASE.PTCs$Row.names
        ASE.PTCs$Row.names <- NULL

        correlations.res[row,"stopgain_PTCs_VS_stopgain_ASE"] <- round(cor.test(ASE.PTCs$stopgain_PTCs,ASE.PTCs$stopgain_ASE, use = "pairwise.complete.obs")$estimate,3)
        correlations.res[row,"stopgain_ASE_VS_synonymous_ASE"] <- round(cor.test(ASE.PTCs$stopgain_ASE,ASE.PTCs$synonymous_ASE, use = "pairwise.complete.obs")$estimate,3)
        }

    }

    return(correlations.res)

}

NMD_eff_correlations <- function(TCGA.cancer, num_PTC) {

    correlations.res <- expand.grid(VAFs,methods)
    colnames(correlations.res) <- c("VAF","methods")
    correlations.res[,c("stopgain_PTCs_VS_stopgain_ASE","endogenous_VS_stopgain_ASE","num_samples")] <- NA

    # 1) Endogenous NMDeff

    # Current endogenous 

    endogenous <- read.table(file = paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/endogenous/",TCGA.cancer,"_endogenous_NB_coeff.txt"),
                                    header = TRUE, sep = "\t", row.names = 1)
    endogenous <- endogenous[,c("NMD.global.2.shared","non.NMD.neg.control")]

    # Correlation between ASE and PTCs for each different VAF and method parameter

    for (VAF in VAFs) {

        # 1) PTCs NMDeff (LOEUF == 4)

        PTCs <- read.table(file = paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/PTCs/",TCGA.cancer,"_PTCs_NB_coeff_gene_id-pca_subtypes-somatic_transcripts_filt-",gsub("e",".E",VAF),"-4.txt"),
                                    header = TRUE, sep = "\t", row.names = 1)
        PTCs <- PTCs[,c("stopgain","synonymous")]
        colnames(PTCs) <- c("stopgain_PTCs","synonymous_PTCs")

        # 2) ASE NMDeff
        for (method in methods) {

            row <- which(correlations.res$VAF%in%VAF & correlations.res$method%in%method)

            ASE.file <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/ASE/",TCGA.cancer,"_ASE_NB_coeff",num_PTC,"_",VAF,"-",method,".txt")
            
            # Current ASE
            ASE <- read.table(file = ASE.file,header = TRUE, sep = "\t", row.names = 1)
            correlations.res[row,"ASE_diff_stopgain_vs_NMD_evading"] <- abs(median(ASE$stopgain,na.rm = TRUE)-median(ASE$stopgain.NMD.evading,na.rm = TRUE))
            correlations.res[row,"num_ASE"] <- median(ASE$num.PTCs,na.rm=TRUE)
            ASE <- ASE[,c("stopgain","synonymous")]
            colnames(ASE) <- c("stopgain_ASE","synonymous_ASE")
            correlations.res[row,"num_samples"] <- sum(!is.na(ASE$stopgain_ASE))  

            # Merge PTCs
            ASE.PTCs <- merge(ASE,PTCs, by.x="row.names",by.y="row.names")
            rownames(ASE.PTCs) <- ASE.PTCs$Row.names
            ASE.PTCs$Row.names <- NULL
            correlations.res[row,"stopgain_PTCs_VS_stopgain_ASE"] <- round(cor.test(ASE.PTCs$stopgain_PTCs,ASE.PTCs$stopgain_ASE, use = "pairwise.complete.obs")$estimate,3)
            correlations.res[row,"synonymous_PTCs_VS_synonymous_ASE"] <- round(cor.test(ASE.PTCs$synonymous_PTCs,ASE.PTCs$synonymous_ASE, use = "pairwise.complete.obs")$estimate,3)
            correlations.res[row,"stopgain_PTCs_VS_synonymous_ASE"] <- round(cor.test(ASE.PTCs$stopgain_PTCs,ASE.PTCs$synonymous_ASE, use = "pairwise.complete.obs")$estimate,3)
            correlations.res[row,"synonymous_PTCs_VS_stopgain_ASE"] <- round(cor.test(ASE.PTCs$synonymous_PTCs,ASE.PTCs$stopgain_ASE, use = "pairwise.complete.obs")$estimate,3)
            correlations.res[row,"stopgain_ASE_VS_synonymous_ASE"] <- round(cor.test(ASE.PTCs$stopgain_ASE,ASE.PTCs$synonymous_ASE, use = "pairwise.complete.obs")$estimate,3)
            # Merge endogenous
            ASE.endogenous <- merge(ASE,endogenous, by.x="row.names",by.y="row.names")
            rownames(ASE.endogenous) <- ASE.endogenous$Row.names
            ASE.endogenous$Row.names <- NULL
            correlations.res[row,"endogenous_VS_stopgain_ASE"] <- round(cor.test(ASE.endogenous$NMD.global.2.shared,ASE.endogenous$stopgain_ASE, use = "pairwise.complete.obs")$estimate,3)
            correlations.res[row,"nonNMD_neg_control_VS_synonymous_ASE"] <- round(cor.test(ASE.endogenous$non.NMD.neg.control,ASE.endogenous$synonymous_ASE, use = "pairwise.complete.obs")$estimate,3)
        }
    }
    return(correlations.res)
}


NMD_eff_ASE_estimates <- function(TCGA.cancer, num_PTC)  {
    NMDeff.res <- expand.grid(VAFs,methods)
    colnames(NMDeff.res) <- c("VAF","methods")
    NMDeff.res[,c("NMDeff_stopgain","NMDeff_synonymous")] <- NA
    # 1) ASE NMDeff
    for (VAF in VAFs) {
        for (method in methods) {
        row <- which(NMDeff.res$VAF%in%VAF & NMDeff.res$method%in%method)

        ASE.file <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/",TCGA.cancer,"/NB_results/ASE/",TCGA.cancer,"_ASE_NB_coeff",num_PTC,"_",VAF,"-",method,".txt")

        # Current ASE
        ASE <- read.table(file = ASE.file,header = TRUE, sep = "\t", row.names = 1)
        # Median NMDeff across samples
        NMDeff.median.ASE <- median(ASE$stopgain, na.rm = TRUE)
        NMDeff.median.synonymous <- median(ASE$synonymous, na.rm = TRUE)
        # Store res
        NMDeff.res[row,"NMDeff_stopgain"] <- NMDeff.median.ASE
        NMDeff.res[row,"NMDeff_synonymous"] <- NMDeff.median.synonymous
        NMDeff.res[row,"num_PTCs"] <- median(ASE$num.PTCs,na.rm=TRUE)
        NMDeff.res[row,"num_samples"] <- sum(!is.na(ASE$stopgain)) 
        }
    }
    return(NMDeff.res)
}

param_opt_corr_plots <- function(cancer, num_PTC) {

    if (cancer == "all") {
        corr_list <- list()
        corr.res.merge <- expand.grid(VAFs,methods)
        colnames(corr.res.merge) <- c("VAF","methods")
        correlations <- c("stopgain_PTCs_VS_stopgain_ASE","endogenous_VS_stopgain_ASE","stopgain_ASE_VS_synonymous_ASE","ASE_diff_stopgain_vs_NMD_evading")
        corr.res.merge[,c(correlations,"num_samples")] <- NA

        for (cncer in c("TCGA-SKCM","TCGA-KIRP","TCGA-GBM","TCGA-PRAD")) {
            corr <- NMD_eff_correlations(TCGA.cancer=cncer,num_PTC = num_PTC)
            corr_list[[cncer]] <- corr
        }
        corr_df <- do.call(cbind,corr_list)
        for (corr in correlations) {
            corr_df_filt <- corr_df[,grep(corr,colnames(corr_df))]
            corr.res.merge[,corr] <- rowMeans(corr_df_filt)
        }
        corr_df_filt <- corr_df[,grep("num_samples",colnames(corr_df))]
        corr.res.merge[,"num_samples"] <- apply(corr_df_filt, 1, FUN = mean)
        corr_df_filt <- corr_df[,grep("num_ASE",colnames(corr_df))]
        corr.res.merge[,"num_ASE"] <- apply(corr_df_filt, 1, FUN = mean)
        corr.df <- corr.res.merge
    } else {
        corr.df <- NMD_eff_correlations(TCGA.cancer=cancer,num_PTC = num_PTC)
    }
    if (num_PTC == "_1PTC") {
        plots.path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Strelka_ASE/",cancer,"/>=1PTCs")
    } else {
        plots.path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Strelka_ASE/",cancer,"/>1PTCs")
    }   

    # Calculate % of missing samples
    corr.df$perc_samples <- round(corr.df$num_samples/max(corr.df$num_samples),2)*100

    # PTC vs ASE correlations
    png(paste0(plots.path,"/corr_stopgain_VAF_method.png"), width = 4500, height = 3000, res = 300)
    tmp <- corr.df[corr.df$methods == "strelka",]
    labels.mut <- paste0(tmp$VAF," (",as.numeric(tmp$num_ASE),")")
    p <- ggplot(data = corr.df, mapping = aes(x = as.factor(VAF), y = eval(parse(text=paste0("stopgain_PTCs_VS_stopgain_ASE"))), 
                group = methods, colour = methods)) +
        scale_fill_hue(name="Methods", labels = c("strelka", "samtools")) +
        scale_x_discrete(labels=labels.mut) +
        geom_text(aes(label=as.character(corr.df$perc_samples)), vjust = 2, size = 10)+
        geom_point(aes(size = corr.df$num_ASE), alpha = 1, group = corr.df$num_ASE) + 
        geom_line(alpha=0.75) + ylim(c(-0.2,0.45)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = paste0("stopgain_PTCs_VS_stopgain_ASE")) +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()

    # Endogenous vs ASE correlations
    png(paste0(plots.path,"/corr_endogenous_VAF_method.png"), width = 4500, height = 3000, res = 300)
    tmp <- corr.df[corr.df$methods == "strelka",]
    labels.mut <- paste0(tmp$VAF," (",as.numeric(tmp$num_ASE),")")
    p <- ggplot(data = corr.df, mapping = aes(x = as.factor(VAF), y = eval(parse(text=paste0("endogenous_VS_stopgain_ASE"))), 
                group = methods, colour = methods)) +
        scale_fill_hue(name="Methods", labels = c("strelka", "samtools")) +
        scale_x_discrete(labels=labels.mut) +
        geom_text(aes(label=as.character(corr.df$perc_samples)), vjust = 2, size = 10)+
        geom_point(aes(size = corr.df$num_ASE), alpha = 1, group = corr.df$num_ASE) + 
        geom_line(alpha=0.75) + ylim(c(-0.2,0.45)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = paste0("endogenous_VS_stopgain_ASE")) +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()

    # PTC vs ASE correlations synonymous
    png(paste0(plots.path,"/corr_synonymous_VAF_method.png"), width = 4500, height = 3000, res = 300)
    tmp <- corr.df[corr.df$methods == "strelka",]
    labels.mut <- paste0(tmp$VAF," (",as.numeric(tmp$num_ASE),")")
    p <- ggplot(data = corr.df, mapping = aes(x = as.factor(VAF), y = eval(parse(text=paste0("stopgain_ASE_VS_synonymous_ASE"))), 
                group = methods, colour = methods)) +
        scale_fill_hue(name="Methods", labels = c("strelka", "samtools")) +
        scale_x_discrete(labels=labels.mut) +
        geom_text(aes(label=as.character(corr.df$perc_samples)), vjust = 2, size = 10, check_overlap = TRUE)+
        geom_point(aes(size = corr.df$num_ASE), alpha = 1, group = corr.df$num_ASE) + 
        geom_line(alpha=0.75) + ylim(c(-0.2,0.45)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = paste0("stopgain_ASE_VS_synonymous_ASE")) +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()

    # ASE stopgain vs PTC NMD-evading difference
    png(paste0(plots.path,"/ASE_diff_stopgain_vs_NMD_evading_VAF_method.png"), width = 4500, height = 3000, res = 300)
    tmp <- corr.df[corr.df$methods == "strelka",]
    labels.mut <- paste0(tmp$VAF," (",as.numeric(tmp$num_ASE),")")
    p <- ggplot(data = corr.df, mapping = aes(x = as.factor(VAF), y = ASE_diff_stopgain_vs_NMD_evading, group = as.factor(methods), colour = as.factor(methods))) +
        geom_point(aes(size = corr.df$num_ASE),alpha = 1) + 
        scale_x_discrete(labels=labels.mut) +
        geom_line(alpha=0.75) + ylim(c(0,3)) +
        geom_text(aes(label=as.character(corr.df$perc_samples)), vjust = 2, size = 10, check_overlap = TRUE)+
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = "PTCs_diff_stopgain_vs_NMD_evading") +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()

}

param_opt_NMDeff_ASE_plots <- function(cancer, num_PTC) {

    if (cancer == "all") {
        NMDeff_list <- list()
        NMDeff.res.merged <- expand.grid(VAFs,methods)
        colnames(NMDeff.res.merged) <- c("VAF","methods")
        NMDeff.res.merged[,c("NMDeff_stopgain","NMDeff_synonymous")] <- NA
        for (cncer in c("TCGA-KIRP","TCGA-SKCM","TCGA-GBM","TCGA-PRAD")) {
            NMDeff.res.cancer <- NMD_eff_ASE_estimates(TCGA.cancer=cncer,num_PTC = num_PTC)
            NMDeff_list[[cncer]] <- NMDeff.res.cancer
        }
        NMDeff_df <- do.call(cbind,NMDeff_list)
        for (NMDeff in c("NMDeff_stopgain","NMDeff_synonymous")) {
            NMDeff_df_filt <- NMDeff_df[,grep(NMDeff,colnames(NMDeff_df))]
            NMDeff.res.merged[,NMDeff] <- rowMeans(NMDeff_df_filt)
        }
        NMDeff_df_filt <- NMDeff_df[,grep("num_samples",colnames(NMDeff_df))]
        NMDeff.res.merged[,"num_samples"] <- apply(NMDeff_df_filt, 1, FUN = mean)
        NMDeff_df_filt <- NMDeff_df[,grep("num_PTCs",colnames(NMDeff_df))]
        NMDeff.res.merged[,"num_PTCs"] <- apply(NMDeff_df_filt, 1, FUN = mean)
        NMDeff.res <- NMDeff.res.merged
    } else {
            NMDeff.res <- NMD_eff_ASE_estimates(TCGA.cancer=cancer,num_PTC = num_PTC)
    }
    if (num_PTC == "_1PTC") {
        plots.path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Strelka_ASE/",cancer,"/>=1PTCs")
    } else {
        plots.path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Strelka_ASE/",cancer,"/>1PTCs")
    }   
   
    # Calculate % of missing samples
    NMDeff.res$perc_samples <- round(NMDeff.res$num_samples/max(NMDeff.res$num_samples),2)*100

    # 1) NMD efficiency using PTCs stopgain
    png(paste0(plots.path,"/NMDeff_ASE_PTCs_by_VAF_methods.png"), width = 4500, height = 3000, res = 300)
    tmp <- NMDeff.res[NMDeff.res$methods == "strelka",]
    labels.mut <- paste0(tmp$VAF," (",as.numeric(tmp$num_PTCs),")")
    p <- ggplot(data = NMDeff.res, mapping = aes(x = as.factor(VAF), y = exp(NMDeff_stopgain), 
                group = as.factor(methods), colour = as.factor(methods))) +
        geom_point(aes(size = NMDeff.res$num_PTCs), alpha = 1, group = NMDeff.res$num_PTCs) + 
        scale_x_discrete(labels=labels.mut) +
        geom_text(aes(label=as.character(NMDeff.res$perc_samples)), vjust = 2, size = 7, check_overlap = TRUE)+
        geom_line(alpha=0.75) + ylim(c(0,1)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = paste0("NMD efficiency values")) +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()
    # 2) NMD efficiency using synonymous (controls)
    png(paste0(plots.path,"/NMDeff_ASE_synonymous_by_VAF_methods.png"), width = 4500, height = 3000, res = 300)
    p <- ggplot(data = NMDeff.res, mapping = aes(x = as.factor(VAF), y = exp(NMDeff_synonymous), 
                group = as.factor(methods), colour = as.factor(methods))) +
        geom_point(size = 6, alpha = 1) +
        geom_text(aes(label=as.character(NMDeff.res$perc_samples)), vjust = 2, size = 7, check_overlap = TRUE)+
        geom_line(alpha=0.75) + ylim(c(0.5,1.2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(color="black", size=15, angle = 45, vjust = 1, hjust=1),
            legend.position='right') +
        labs(y = paste0("NMD efficiency values")) +
        scale_size_continuous(range = c(2, 10), breaks = seq(0,25,by=5))
    print(p)
    dev.off()
}

# Define data
VAFs <- c("1e-05","1e-04","0.001","0.01","0.1","0.2","0.3","0.4","0.5")
methods <- c("strelka","samtools")

#for (cancer in c("all","TCGA-KIRP","TCGA-LIHC","TCGA-GBM","TCGA-SKCM","TCGA-PRAD")) {
for (cancer in c("all","TCGA-PRAD","TCGA-GBM","TCGA-SKCM","TCGA-KIRP")) {
    print(cancer)
    # Correlation plots
    param_opt_corr_plots(cancer = cancer, num_PTC = "_1PTC")
    #param_opt_corr_plots(cancer = cancer, num_PTC = "")
    # NMD efficiency plots
    param_opt_NMDeff_ASE_plots(cancer = cancer, num_PTC = "_1PTC")
    #param_opt_NMDeff_ASE_plots(cancer = cancer, num_PTC = "")
}





