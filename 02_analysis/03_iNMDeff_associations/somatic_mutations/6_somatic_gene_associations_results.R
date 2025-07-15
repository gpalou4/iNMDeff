lambda <- function(p_values) {
    # Lambda calculation
    chisq <- qchisq(1-p_values,1)
    lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
    return(lambda)
}

QC_FDR_and_qvalue_corrections <- function(df_ass_final_mut,mutation) {
    error <<- FALSE
    pval_char <- paste0("som_",mutation,"_pval")
    # Q-value (FDR) correction using "John Storey" method
    tryCatch({
        pvals <- as.numeric(df_ass_final_mut[,pval_char])
        qvals <- qvalue(pvals)$qvalues
        pvals_adjust <- p.adjust(pvals,method="fdr")
        df_ass_final_mut[,paste0("som_",mutation,"_qval")] <- qvals
        df_ass_final_mut[,paste0("som_",mutation,"_pval_adjusted")] <- pvals_adjust
        # Add top hits column
        df_ass_final_mut <- df_ass_final_mut[order(df_ass_final_mut[,paste0("som_",mutation,"_qval")]),]
        rows <- which(df_ass_final_mut[,paste0("som_",mutation,"_qval")] <= 0.1)
        df_ass_final_mut[,c("top_hits","top_hits_type","color")] <- NA
        df_ass_final_mut[rows,"top_hits"] <- df_ass_final_mut[rows,"Gene_symbol"]
        df_ass_final_mut[rows,"top_hits_type"] <- df_ass_final_mut[rows,"Role.in.Cancer"]
        df_ass_final_mut$NMD_type <- ifelse(!is.na(df_ass_final_mut$NMD_type),"NMD",NA)
        df_ass_final_mut[rows,"color"] <- gsub("_?NA_?","",paste0(df_ass_final_mut[rows,"top_hits_type"],"_",df_ass_final_mut[rows,"NMD_type"]))
        }, error = function(e){
            df_ass_final_mut$som_CNV_del_qval <<- df_ass_final_mut[,paste0("som_",mutation,"_coeff")]  
            error <<- TRUE
            print(e)
            return("NA")
        }
    )
    if (isTRUE(error)) {
        return("NA")
    } else {
        return(df_ass_final_mut)
    }   
}
            
CGC_associations <- function(TCGA_cancers, NMD_method, num_PCs = NULL, cancer_types_corrected_CNV_PCs) {

    if (NMD_method == "ASE") {
        NMD_method_VAF <- paste0(NMD_method,"_",VAF)
        grep_char <- "ASE"
    } else {
        NMD_method_VAF <- NMD_method
        grep_char <- "END"
    }
    
    # PCs optimized for each cancer
    best_PC_all_TCGA_cancers <- read.table(file = paste0("~/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/",NMD_method_VAF,"_best_PC_all_TCGA_cancers.txt"), header = TRUE, sep ="\t")

    TCGA_cancers_list <- list()
    for (TCGA_cancer in TCGA_cancers) {TCGA_cancers_list[[TCGA_cancer]] <- NA}
    top_hits <- list(CNV_amp = NA, CNV_del = NA, mut_missense = NA, mut_truncating = NA, mut_synonymous = NA)
    for (mut in names(top_hits)) {top_hits[[mut]] <- TCGA_cancers_list}
    if ("pancancer" %in% TCGA_cancers) {
        char3 <- "pancancer"
    } else {
        char3 <- "by_cancers"
    }
    for (TCGA_cancer in TCGA_cancers) {

        print(TCGA_cancer)
        if (is.null(num_PCs)) {
            # Corrected cancer type or not?
            row1 <- cancer_types_corrected_CNV_PCs$cancer_type %in% gsub("TCGA-","",TCGA_cancer)
            col <- paste0(grep_char,"_CNV")
            cancer_type_correct <- as.character(cancer_types_corrected_CNV_PCs[row1,col])
            if ( cancer_type_correct == 1) {
                PCs <- 1
            } else {
                PCs <- 0
            }
        }

        # CGC results
        CGC_results_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/final_associations/",TCGA_cancer,"_",NMD_method_VAF,"_",char,"_PCs_",PCs,".txt")
        row2 <- best_PC_all_TCGA_cancers$TCGA_cancer %in% TCGA_cancer
        df_ass <- read.csv(file = CGC_results_path, header = TRUE, stringsAsFactors = FALSE, sep ="\t")
        df_ass_final <- merge(df_ass,CGC, by.x = c("Gene_symbol","ENSEMBL_gene"), by.y = c("Gene.Symbol","Synonyms"), all.x = TRUE)
        df_ass_final <- merge(df_ass_final, ensembl_v88_gtf_filt, by.x = "Gene_symbol", by.y ="gene_name", all.x = TRUE)
        df_ass_final$chr <- gsub("(.*)\\:.*","\\1",df_ass_final$genome_location)
        # Fix TSG/OG label
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "oncogene, fusion","oncogene",df_ass_final$Role.in.Cancer)
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "TSG, fusion","TSG",df_ass_final$Role.in.Cancer)
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "oncogene, TSG, fusion","other",df_ass_final$Role.in.Cancer)
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "oncogene, TSG","other",df_ass_final$Role.in.Cancer)
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "",NA,df_ass_final$Role.in.Cancer)
        df_ass_final$Role.in.Cancer <- ifelse(df_ass_final$Role.in.Cancer == "fusion","other",df_ass_final$Role.in.Cancer)
        CGC_genes <- CGC[which(CGC$Role.in.Cancer != "non_CGC_NMD"),"Gene.Symbol"]
        #NMD_genes <- CGC[which(!is.na(CGC$NMD_type)),"Gene.Symbol"]
        random_genes <- CGC[which(CGC$Role.in.Cancer == "non_CGC_NMD"),"Gene.Symbol"]

        for (mutation in c("CNV_amp","CNV_del","mut_missense","mut_truncating","mut_synonymous")) {
            print(mutation)
            if ( sum(na.omit(df_ass_final[,paste0("som_",mutation,"_coeff")])) == 0 ) {
                next
            }
            # Number of PCs for each mutation
            if (mutation == "CNV_amp" || mutation == "CNV_del") {
                CNV_PCs <- best_PC_all_TCGA_cancers[row2,"best_PC_som_CNV"]
            } else if ((mutation == "mut_missense" || mutation == "mut_truncating")) {
                CNV_PCs <- best_PC_all_TCGA_cancers[row2,"best_PC_som_mut_miss_tr"]
            } else if (mutation == "mut_synonymous") {
                CNV_PCs <- best_PC_all_TCGA_cancers[row2,"best_PC_som_mut_synonymous"]
            }

            qqplot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/QC/",char3,"/",TCGA_cancer,"_",NMD_method_VAF,"_qqplot_",mutation,"_PCs_",PCs,".png")
            pvals_histogram <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/QC/",char3,"/",TCGA_cancer,"_",NMD_method_VAF,"_histogram_pvals_",mutation,"_PCs_",PCs,".png")
            volcano_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/QC/",char3,"/",TCGA_cancer,"_",NMD_method_VAF,"_volcano_",mutation,"_PCs_",PCs,".png")
            manhattan_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/QC/",char3,"/",TCGA_cancer,"_",NMD_method_VAF,"_manhattan_",mutation,"_PCs_",PCs,".png")

            # Remove pals with NAs
            pval_char <- paste0("som_",mutation,"_pval")
            pval_exp_char <- paste0("som_",mutation,"_pval_exp")
            df_ass_final_filt <- df_ass_final[-which(is.na(df_ass_final[,pval_char])),]
            print(dim(df_ass_final_filt))

            # 1) QQ-plots
            # Calculate Lambda
            p_values_all <- df_ass_final_filt[,pval_char]
            p_values_cancer_genes <- df_ass_final_filt[df_ass_final_filt$Gene_symbol %in% c(CGC_genes,NMD_genes$gene_symbol),pval_char]
            p_values_random_genes <- df_ass_final_filt[df_ass_final_filt$Gene_symbol %in% random_genes,pval_char]
            lambda_all_genes <- lambda(p_values_all)
            lambda_cancer_genes <- lambda(p_values_cancer_genes)
            lambda_random_genes <- lambda(p_values_random_genes)
            df_ass_final_filt$color_qqplot <- factor(df_ass_final_filt$Role.in.Cancer)

            # 1.1) Cancer + NMD genes
            plot_title <- paste0(mutation," for Cancer + NMD genes ---- Lambda = ",lambda_cancer_genes)
            df_cancer_NMD_genes <- df_ass_final_filt[df_ass_final_filt$Gene_symbol %in% c(CGC_genes,NMD_genes$gene_symbol),]
            df_cancer_NMD_genes <- df_cancer_NMD_genes[order(df_cancer_NMD_genes[,pval_char]),]
            # Expected p-values
            n1 <- length(p_values_cancer_genes)
            if (n1 != 0) {
                a <- 1:n1
                pvals_exp <- a/n1
                df_cancer_NMD_genes[,paste0("som_",mutation,"_pval_exp")] <- pvals_exp
                # QQplot       
                p1 <- ggplot(df_cancer_NMD_genes, aes(x = -log10(eval(parse(text=pval_exp_char))), y = -log10(eval(parse(text=pval_char))),color = color_qqplot)) +
                geom_point(size = 2) +
                geom_abline(intercept = 0, slope = 1) +
                ylab("-log10(P-values") + xlab("-log10(expected P-values") + ggtitle(plot_title) +
                theme_classic(base_size = 16)
            }
            # 1.2) Random genes
            plot_title <- paste0(mutation," for Random genes ---- Lambda = ",lambda_random_genes)
            df_random_genes <- df_ass_final_filt[df_ass_final_filt$Gene_symbol %in% random_genes,]
            df_random_genes <- df_random_genes[order(df_random_genes[,pval_char]),]
            # Expected p-values
            n2 <- length(p_values_random_genes)
            if (n2 != 0) {
                a <- 1:n2
                pvals_exp <- a/n2
                df_random_genes[,paste0("som_",mutation,"_pval_exp")] <- pvals_exp
                # QQplot     
                p2 <- ggplot(df_random_genes, aes(x = -log10(eval(parse(text=pval_exp_char))), y = -log10(eval(parse(text=pval_char))))) +
                geom_point(size = 2) +
                geom_abline(intercept = 0, slope = 1) +
                ylab("-log10(P-values") + xlab("-log10(expected P-values") + ggtitle(plot_title) +
                theme_classic(base_size = 16)
            }
            if ( (n1 != 0) & (n2 != 0) ) {
                # 1.3) Final qqplot
                png(qqplot, width = 5000, height = 2500, res = 300)
                p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
                print(p)
                dev.off()
            }
            #### OLD Plot ####
            # title <- paste0("QQ-plot for --> ",mutation," and Total Lambda = ",lambda_all_genes," \n
            #                 Cancer genes Lambda = ",lambda_cancer_genes," \n
            #                 Random genes Lambda = ",lambda_random_genes)
            # p <- qqPlot(p_values_all, truncate = FALSE, ylim = NULL, thinThreshold = NULL, 
            #             ci=TRUE, main = title, col = as.character(df_ass_final$color_qqplot))
            # legend(x = "topleft", legend=c("cancer","random","NMD"),
            #             fill = 1:3, cex=1.2)        
            ##################

            # 2) Histogram of p-vals
            # 2.1) Cancer + NMD genes
            p1 <- ggplot(df_cancer_NMD_genes, aes(x = eval(parse(text=pval_char)))) +
                geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 15) +
                geom_vline(aes(xintercept=median(eval(parse(text=pval_char)), na.rm = TRUE)),color="red", linetype="dashed", size=1) + 
                xlab("p-values") + ylab("frequency") + ggtitle(paste0(mutation," --- Cancer + NMD Genes")) +
                theme_classic(base_size = 16)
            # 2.2) Random genes
            p2 <- ggplot(df_random_genes, aes(x = eval(parse(text=pval_char)))) +
                geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 15) +
                geom_vline(aes(xintercept=median(eval(parse(text=pval_char)), na.rm = TRUE)),color="red", linetype="dashed", size=1) + 
                xlab("p-values") + ylab("frequency") + ggtitle(paste0(mutation," --- Random Genes")) +
                theme_classic(base_size = 16)
            # 2.3) Final histogram
            png(pvals_histogram, width = 5000, height = 2500, res = 300)
            p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
            print(p)
            dev.off()

            # # Remove missing values
            # NAs <- is.na(df_ass_final[,pval_char])
            # df_ass_final_mut <- df_ass_final[!NAs,]

            # 3) Volcano
            df_cancer_NMD_genes <- QC_FDR_and_qvalue_corrections(df_ass_final_mut = df_cancer_NMD_genes, mutation = mutation)
            if (df_cancer_NMD_genes == "NA") {
                next
            }
            # Save top hits
            sig_hits <- which(df_cancer_NMD_genes[,paste0("som_",mutation,"_qval")] < 1)
            if (is.null(PCs)) {
                cancer_type_correct_mut <- as.character(cancer_types_corrected_CNV_PCs[row1,paste0(grep_char,"_",gsub("_amp|_del","",mutation))])
            } else {
                cancer_type_correct_mut <- "ok"
            }
            if ( (length(sig_hits) != 0) & (cancer_type_correct_mut != "remove") ) {
                df_ass_final_hits <- QC_FDR_and_qvalue_corrections(df_ass_final_mut = df_ass_final_filt, mutation = mutation)
                top_hits[[mutation]][TCGA_cancer][[1]] <- df_ass_final_hits[,grep("coeff|Gene_symbol|pval|qval",colnames(df_ass_final_hits))]
            }
            df_random_genes <- QC_FDR_and_qvalue_corrections(df_ass_final_mut = df_random_genes, mutation = mutation)

            if (!isTRUE(error)) {
                xlim1 <- min(na.omit(df_cancer_NMD_genes[,paste0("som_",mutation,"_coeff")]))-0.05
                xlim2 <- max(na.omit(df_cancer_NMD_genes[,paste0("som_",mutation,"_coeff")]))+0.05        
                # 3) Volcano plot
                # 3.1) Cancer + NMD genes
                p1 <- ggplot(df_cancer_NMD_genes, aes(x = eval(parse(text=paste0("som_",mutation,"_coeff"))), y = -log10(eval(parse(text=paste0("som_",mutation,"_qval")))), fill = color)) +
                    geom_point(size = 2, alpha = 0.5) + geom_label_repel(aes(label=top_hits),hjust=0.5, vjust=0.5, size = 2, max.overlaps = nrow(df_cancer_NMD_genes), alpha = 1) +
                    geom_hline(yintercept=-log10(0.10), linetype="dashed", color = "red") +
                    geom_vline(xintercept=0, linetype="dashed", color = "red") +
                    xlab("LM Coefficients") + ylab("-log10(q-value)") + ggtitle(paste0(mutation,"")) +
                        xlim(c(xlim1,xlim2)) +
                        scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
                        scale_fill_brewer(type = "qual", palette = "Dark2") +
                        scale_color_brewer(type = "qual",palette = "Paired") +
                        theme_classic(base_size = 16) +
                        theme(legend.position="top")
                # 3.2) Random genes
                p2 <- ggplot(df_random_genes, aes(x = eval(parse(text=paste0("som_",mutation,"_coeff"))), y = -log10(eval(parse(text=paste0("som_",mutation,"_qval")))), fill = color)) +
                    geom_point(size = 2, alpha = 0.5) + geom_label_repel(aes(label=top_hits),hjust=0.5, vjust=0.5, size = 2, max.overlaps = nrow(df_random_genes), alpha = 1) +
                    geom_hline(yintercept=-log10(0.10), linetype="dashed", color = "red") +
                    geom_vline(xintercept=0, linetype="dashed", color = "red") +
                    xlab("LM Coefficients") + ylab("-log10(q-value)") + ggtitle(paste0(mutation,"")) +
                        xlim(c(xlim1,xlim2)) +
                        scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
                        scale_fill_brewer(type = "qual", palette = "Dark2") +
                        scale_color_brewer(type = "qual",palette = "Paired") +
                        theme_classic(base_size = 16) +
                        theme(legend.position="top")
                # 3.3) Final volcano
                png(volcano_plot, width = 5000, height = 2500, res = 300)
                p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
                print(p)
                dev.off()

                # 4) Manhattan plot
                for (gene_type in c("cancer_NMD","random")) {
                    if (gene_type == "cancer_NMD") {
                        df_ass_final_mut <- df_cancer_NMD_genes
                    } else if (gene_type == "random") {
                        df_ass_final_mut <- df_random_genes
                    }
                    NAs <- which(is.na(df_ass_final_mut$chr) | (df_ass_final_mut$chr == "M"))
                    if (length(NAs) != 0 ) {
                        df_ass_final_mut <- df_ass_final_mut[-NAs,]
                    }
                    if (nrow(df_ass_final_mut) == 0) {
                        next
                    } else {
                        # Order by chr
                        skip <- c(which(is.na(df_ass_final_mut$chr)),which(df_ass_final_mut$chr %in% c("X","Y")))
                        if (length(skip) != 0) {
                            tmp <- df_ass_final_mut[skip,]
                            tmp <- tmp[order(as.character(tmp$chr)),]
                            df_ass_final_mut <- df_ass_final_mut[-skip,]
                        } else {
                            tmp <- data.frame()
                        }
                        df_ass_final_mut <- df_ass_final_mut[order(as.numeric(as.character(df_ass_final_mut$chr))),]
                        if (nrow(tmp) == 0) {
                            df_ass_final_mut <- df_ass_final_mut
                        } else {
                            df_ass_final_mut <- rbind(df_ass_final_mut,tmp)
                        }
                        df_ass_final_mut$chr <- ifelse(df_ass_final_mut$chr == "X","23",df_ass_final_mut$chr)
                        df_ass_final_mut$chr <- ifelse(df_ass_final_mut$chr == "Y","24",df_ass_final_mut$chr)
                        df_ass_final_mut$chr <- as.numeric(df_ass_final_mut$chr)
                        # Create BP column
                        bp_all <- c()
                        for (chr in unique(df_ass_final_mut$chr)) {
                            bp <- sum(df_ass_final_mut$chr == chr)
                            bp_vec <- seq(1:bp)
                            bp_all <- c(bp_all,bp_vec)
                        }
                        df_ass_final_mut$BP <- bp_all
                        # Label NMD genes
                        NMD_label <- df_ass_final_mut[df_ass_final_mut$Gene_symbol %in% NMD_genes$gene_symbol,"Gene_symbol"]
                        png(gsub("endogenous_",paste0("endogenous_",gene_type,"_"),manhattan_plot), width = 3500, height = 2500, res = 300)
                        manhattan(df_ass_final_mut, chr="chr", bp="BP", snp="Gene_symbol", p=paste0("som_",mutation,"_qval"),
                                    annotatePval = 1, highlight = NMD_label, genomewideline = -log10(1e-01), main = paste0(TCGA_cancer))
                        dev.off()
                    }

                }
            } # VOLCANO + MANHATTAN plot
        } # FOR EACH MUTATION TYPE
    }
    return(top_hits)
}

volcano_plot <- function(all_CGC_assoc, NMD_method_VAF, pval_adjust, FDR, type, PCs) {

    # All hits across cancers
    top_hits_list <- lapply(all_CGC_assoc, function(x) {
        do.call(rbind,x)
    })

    # For each mutation type
    for (mut in c("CNV_amp","CNV_del","mut_missense","mut_truncating","mut_synonymous")) {

        if (NMD_method_VAF == "ASE_0.2") {
            xlim <- c(-2.75,2.75)
            if (mut == "CNV_amp") {
                ylim <- c(0,2.5)
                if (type == "pancancer") {
                    xlim <- c(-0.2,0.2)
                    ylim <- c(0,2)
                }
            } else if (mut == "CNV_del") {
                ylim <- c(0,3)
                if (type == "pancancer") {
                    xlim <- c(-0.2,0.2)
                    ylim <- c(0,1.75)
                }
            } else if (mut == "mut_missense") {
                ylim <- c(0,2)
                if (type == "pancancer") {
                    xlim <- c(-0.5,0.5)
                    ylim <- c(0,1.25)
                }
            } else if (mut == "mut_truncating") {
                ylim <- c(0,2)
                if (type == "pancancer") {
                    xlim <- c(-1,1)
                    ylim <- c(0,1.25)
                }
            } else if (mut == "mut_synonymous") {
                ylim <- c(0,5)
                if (type == "pancancer") {
                    xlim <- c(-1,1)
                    ylim <- c(0,1.25)
                }
            }
        } else if (NMD_method_VAF == "endogenous") {
            xlim <- c(-2.75,2.75)
            if (mut == "CNV_amp") {
                ylim <- c(0,10)
                ybreak <- c(8,18)
                if (type == "pancancer") {
                    xlim <- c(-0.2,0.2)
                    ylim <- c(0,3)
                }
            } else if (mut == "CNV_del") {
                ylim <- c(0,12)
                ybreak <- c(13,20)
                if (type == "pancancer") {
                    xlim <- c(-0.2,0.2)
                    ylim <- c(0,3)
                }
            } else if (mut == "mut_missense") {
                xlim <- c(-1,1)
                ylim <- c(0,5)
                ybreak <- c(7,10)
                if (type == "pancancer") {
                    xlim <- c(-0.25,0.25)
                    ylim <- c(0,2.5)
                }
            } else if (mut == "mut_truncating") {
                ylim <- c(0,5)
                xlim <- c(-1,1)
                ybreak <- c(3,4.5)
                if (type == "pancancer") {
                    xlim <- c(-0.5,0.5)
                    ylim <- c(0,2.5)
                }
            } else if (mut == "mut_synonymous") { 
                ylim <- c(0,5)
                xlim <- c(-1,1)
                ybreak <- c(7,10)
                if (type == "pancancer") {
                    xlim <- c(-0.75,0.75)
                    ylim <- c(0,2.5)
                }
            }  
        } 
        if (type == "by_cancers") {
            xlimFisher <- -2
            ylimFisher <- 1.5
        } else if (type == "pancancer") {
            xlimFisher <- -0.15
            ylimFisher <- 1.5
        }

        top_hits_df <- top_hits_list[[mut]]
        top_hits_df <- top_hits_df[!is.na(top_hits_df$Gene_symbol),]
        top_hits_df$cancer <- gsub("(.*)\\..*","\\1",gsub(".*\\.(.*)\\..*","\\1",rownames(top_hits_df)))
        top_hits_df$gene_cancer <- paste0(top_hits_df$Gene_symbol,"_",top_hits_df$cancer)
        # NMD
        top_hits_df <- merge(top_hits_df,NMD_genes, by.x="Gene_symbol", by.y = "gene_symbol", all.x = TRUE)

        # Filter by p-value
        #pval_adjust <- "qval"
        #pval_adjust <- "pval_adjusted"
        sig_hits <- which(top_hits_df[,paste0("som_",mut,"_",pval_adjust)] < FDR)

        # Fisher test
        top_hits_df$sig <- FALSE
        top_hits_df[sig_hits,"sig"] <- TRUE
        top_hits_df$NMD_factor <- FALSE
        top_hits_df[which(top_hits_df$NMD_type %in% c("NMD_factor","EJC","DECID_complex")),"NMD_factor"] <- TRUE
        tryCatch({
                fisher_test <<- fisher.test(table(top_hits_df$NMD_factor, top_hits_df$sig))
            }, error = function(e){
                fisher_test <<- data.frame(p.value = 0, estimate = 0)
            print(e)
            }
        )
        # Volcano plot
        top_hits_df <- top_hits_df[order(top_hits_df[,paste0("som_",mut,"_",pval_adjust)]),]
        top_hits_df$top_hits <- NA
        #filter <- ( top_hits_df$sig ) & ( top_hits_df$NMD_type %in% c("NMD_factor","EJC","DECID_complex") )
        filter <- ( top_hits_df$sig )
        top_hits_df[filter,"top_hits"] <- top_hits_df[filter,"gene_cancer"] 
        top_hits_df$NMD_type <- ifelse(is.na(top_hits_df$NMD_type),"TSG/OG",top_hits_df$NMD_type)
        top_hits_df$top_hits <- factor(top_hits_df$top_hits)
        top_hits_df <- top_hits_df[!duplicated(top_hits_df),]
        top_hits_df <- top_hits_df[!abs(top_hits_df[,paste0("som_",mut,"_coeff")]) > 100,]

        volcano_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/volcanos/",type,"/",type,"_",NMD_method_VAF,"_volcano_",mut,"_",pval_adjust,"_",FDR,"_PCs_",PCs,".png")
        png(volcano_plot, width = 3500, height = 2500, res = 300)
        p <- ggplot(top_hits_df, aes(x = eval(parse(text=paste0("som_",mut,"_coeff"))), y = -log10(eval(parse(text=paste0("som_",mut,"_",pval_adjust)))), color = NMD_type)) +
            geom_point(size = 1, alpha = 0.5) + geom_label_repel(aes(label=top_hits),hjust=0.5, vjust=0.5, size = 3, max.overlaps = nrow(top_hits_df), alpha = 1) +
            geom_hline(yintercept=-log10(FDR), linetype="dashed", color = "red") +
            geom_vline(xintercept=0, linetype="dashed", color = "red") +
            #geom_vline(xintercept=0.5, linetype="dashed", color = "red") +
            xlab("LM Coefficients") + ylab("-log10(adjusted p-value)") + ggtitle(paste0("Somatic variant type --> ",mut)) +
            theme(legend.position="top", plot.title = element_text(hjust = 0.5, size = 18),
                axis.title.x = element_text(color="black", size=18, face="bold"),
                axis.title.y = element_text(color="black", size=18, face="bold"),
                legend.title = element_text(colour="black", size=18, face="bold")) +
                theme_classic() + #xlim(xlim) + ylim(ylim) + #facet_zoom(ylim = c(0, 5)) + #scale_y_break(ybreak) + # + 
                geom_text(x=xlimFisher, y=ylimFisher, color = "black",
                            label=paste0("NMD enrichment (Fisher test) \n p-value: ",round(fisher_test$p.value,2),"\n OR: ",round(fisher_test$estimate,2))) +
                #scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
                scale_fill_brewer(type = "qual", palette = "Dark2")+#, name = "Cancer Gene Census", labels = unique(df_ass_final$top_hits_type)) +
                scale_color_brewer(type = "qual",palette = "Paired")#, name = "NMD genes", labels = unique(df_ass_final$NMD_type))
        print(p)
        #plot_grid(plotlist=plots_list, labels = "AUTO", align = "v", ncol = 3, nrow = 1)
        dev.off()
    }
}

freq_genes_plots <- function(all_CGC_assoc_end, all_CGC_assoc_ASE, pval_adjust, FDR, type) {

    top_hits_list_end <- lapply(all_CGC_assoc_end, function(x) {
        do.call(rbind,x)
    })
    top_hits_list_ASE <- lapply(all_CGC_assoc_ASE, function(x) {
        do.call(rbind,x)
    })

    for (mut in c("CNV_amp","CNV_del","mut_missense","mut_truncating","mut_synonymous")) {

        if (type == "by_cancers") {      
            if (mut == "CNV_amp") {
                xlimNMD <- c(-2,3)
            } else if (mut == "mut_synonymous") {
                xlimNMD <- c(-2,2)
            } else {
                xlimNMD <- c(-1,1)
            }
            xlimAll <- c(-2.5,2.5)
            xlimFisher <- 15
            ylimFisher <- -2
        } else if (type == "pancancer") {
            xlimNMD <- c(-0.4,0.4)
            xlimAll <- c(-0.4,0.4)
            xlimFisher <- 15
            ylimFisher <- -0.15
        }
        top_hits_df_end <- top_hits_list_end[[mut]]
        top_hits_df_end$NMD_method <- "endogenous"
        top_hits_df_ASE <- top_hits_list_ASE[[mut]]
        top_hits_df_ASE$NMD_method <- "ASE"
        # Check sig hits in ASE (discovery)
        sig_hits_genes <- top_hits_df_ASE[top_hits_df_ASE[,paste0("som_",mut,"_",pval_adjust)] < FDR,"Gene_symbol"]
        # Validate sig hits in Endogenous (validation)
        top_hits_df_end_tmp <- top_hits_df_end[top_hits_df_end$Gene_symbol %in% sig_hits_genes,]
        sig_hits_genes_validated <- top_hits_df_end_tmp[top_hits_df_end_tmp[,paste0("som_",mut,"_",pval_adjust)] < FDR,"Gene_symbol"]
        # Merge
        top_hits_df <- rbind(top_hits_df_end,top_hits_df_ASE)
        top_hits_df <- top_hits_df[!is.na(top_hits_df$Gene_symbol),]
        top_hits_df$cancer <- gsub("(.*)\\..*","\\1",gsub(".*\\.(.*)\\..*","\\1",rownames(top_hits_df)))
        top_hits_df$gene_cancer <- paste0(top_hits_df$Gene_symbol,"_",top_hits_df$cancer)
        # NMD
        top_hits_df <- merge(top_hits_df,NMD_genes, by.x="Gene_symbol", by.y = "gene_symbol", all.x = TRUE)
        # Fisher test
        sig_hits <- which(top_hits_df$Gene_symbol %in% sig_hits_genes_validated)
        top_hits_df$sig <- FALSE
        top_hits_df[sig_hits,"sig"] <- TRUE
        top_hits_df$NMD_factor <- FALSE
        top_hits_df[which(top_hits_df$NMD_type %in% c("NMD_factor","EJC","DECID_complex")),"NMD_factor"] <- TRUE
        tryCatch({
                fisher_test <<- fisher.test(table(top_hits_df$NMD_factor, top_hits_df$sig))
            }, error = function(e){
                fisher_test <<- data.frame(p.value = 0, estimate = 0)
            print(e)
            }
        )
        # Sign of coefficient
        top_hits_df$coeff_sign <- NA
        top_hits_df[top_hits_df[,paste0("som_",mut,"_coeff")] > 0,"coeff_sign"] <- "pos"
        top_hits_df[top_hits_df[,paste0("som_",mut,"_coeff")] < 0,"coeff_sign"] <- "neg"
        # Filter by significant validated genes
        top_hits_df_filt <- top_hits_df[sig_hits,]
        if (nrow(top_hits_df_filt) == 0 ){
            next
        }
        medians <- top_hits_df_filt %>% group_by(Gene_symbol) %>% summarize(median = median(eval(parse(text=paste0("som_",mut,"_coeff")))))
        top_hits_df_filt <- left_join(top_hits_df_filt, medians, by = "Gene_symbol")
        top_hits_df_filt <- top_hits_df_filt[!duplicated(top_hits_df_filt),]
        # Add frequency to the plot
        tmp <- data.frame(Gene_symbol = table(top_hits_df_filt$Gene_symbol))
        colnames(tmp) <- c("Gene_symbol","frequency")
        top_hits_df_filt <- left_join(top_hits_df_filt, tmp, by = "Gene_symbol")
        top_hits_df_filt$Gene_symbol_ss <- paste0(top_hits_df_filt$Gene_symbol," (",top_hits_df_filt$frequency,")")
        # Save results
        output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/freq_hits/",type,"/",type,"_NMD_factors_results_",mut,"_",pval_adjust,"_",FDR,".txt")
        res <- top_hits_df_filt[top_hits_df_filt$NMD_factor,]
        write.table(res, file = output_path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

        # 3.1) Show only NMD factors, regardless if replicated and/or high freq
        top_hits_df_NMD <- top_hits_df_filt[top_hits_df_filt$NMD_factor,]
        freq_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/freq_hits/",type,"/",type,"_freq_boxplot_NMD_",mut,"_",pval_adjust,"_",FDR,".png")
        png(freq_plot, width = 3500, height = 2500, res = 300)
        p <- ggplot(data = top_hits_df_NMD, aes(x=reorder(Gene_symbol_ss, median), y = eval(parse(text=paste0("som_",mut,"_coeff"))), color = NMD_method)) +
            coord_flip() + #coord_flip(ylim = xlimNMD) +
            geom_point(size = 1) +
            geom_boxplot(width=1, fill="white")+
            geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
            ylab("Somatic variant LM coefficient") + xlab("Gene Symbol") + ggtitle(paste0(mut)) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(color="black", size=18),
                axis.text.y = element_text(color="black", size=18),
                legend.position = "top", axis.text=element_text(size=18))
        print(p)
        dev.off()

        # 3.2) Show only most frequent genes replicated in ASE and End
        top_hits_df_ASE <- top_hits_df_filt[top_hits_df_filt$NMD_method %in% "ASE",]
        top_hits_df_End <- top_hits_df_filt[top_hits_df_filt$NMD_method %in% "endogenous",]
        replicated_genes <- intersect(top_hits_df_ASE$Gene_symbol,top_hits_df_End$Gene_symbol)
        top_hits_df_filt$replicated_gene <- FALSE
        top_hits_df_filt[top_hits_df_filt$Gene_symbol %in% replicated_genes,"replicated_gene"] <- TRUE
        # Are there more NMD factor genes replicated between ASE and End than not replicated?
        tryCatch({
                fisher_test <- fisher.test(table(top_hits_df_filt$NMD_factor, top_hits_df_filt$replicated_gene))
            }, error = function(e){
                fisher_test <<- data.frame(p.value = 0, estimate = 0)
            print(e)
            }
        )
        top_hits_df_all <- top_hits_df_filt[top_hits_df_filt$replicated_gene,]
        top_hits_df_all <- top_hits_df_all[!duplicated(top_hits_df_all),]
        
        if (nrow(top_hits_df_all) == 0 ) {
            next
        }
        # Order by frequency (by cancers) or by p-value (pancancer)
        # if (type == "by_cancers") {
        #     top_hits_df_all <- top_hits_df_all[order(top_hits_df_all$frequency, decreasing = TRUE),] 
        # } else if (type == "pancancer") {
        #     top_hits_df_all <- top_hits_df_all[order(top_hits_df_all[,paste0("som_",mut,"_",pval_adjust)]),]
        # }  
        top_hits_df_all <- top_hits_df_all[1:200,]
        na_rows <- grep("NA",rownames(top_hits_df_all))
        if (length(na_rows) != 0) {
            top_hits_df_all <- top_hits_df_all[-grep("NA",rownames(top_hits_df_all)),]
        }
        top_hits_df_all$NMD_type <- ifelse(is.na(top_hits_df_all$NMD_type),"other",top_hits_df_all$NMD_type)
        # Order by median coefficient
        top_hits_df_all <- top_hits_df_all[order(top_hits_df_all$median),]
        # Color NMD factors in X axis using an external vector 
        odd_indexes <- seq(1, length(top_hits_df_all$Gene_symbol), by = 2) 
        odd_elements <- factor(top_hits_df_all$NMD_type[odd_indexes]) 
        #x = reorder(Gene_symbol_ss,median)
        freq_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/",char2,"/analysis_plots/",char3,"freq_hits/",type,"/",type,"_freq_boxplot_all_",mut,"_",pval_adjust,"_",FDR,".png")
        png(freq_plot, width = 3000, height = 3500, res = 300)
        p <- ggplot(data = top_hits_df_all, aes(x=reorder(Gene_symbol,median), y = eval(parse(text=paste0("som_",mut,"_coeff"))), color = NMD_method)) +
            coord_flip() + #coord_flip(ylim = xlimAll)
            geom_point(size = 1) +
            geom_boxplot(width=1, fill="white")+
            geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) +
            ylab("Somatic variant LM coefficient") + xlab("Gene Symbol") + ggtitle(paste0(mut)) +
            theme_classic() + 
            geom_text(x=xlimFisher, y=ylimFisher, color = "black",
            label=paste0("NMD enrichment (Fisher test) \n p-value: ",round(fisher_test$p.value,2),"\n OR: ",round(fisher_test$estimate,2))) +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(size=16, colour = "black"),
                axis.text.y = element_text(size=16, colour = odd_elements),
                legend.position = "top", axis.text=element_text(size=18))
        print(p)
        dev.off()
    }
}

library("qvalue")
library("ggplot2")
library("ggrepel")
library("RColorBrewer")
library("GWASTools")
# Manhattan plot
library("qqman")
library("pheatmap")
library("dplyr")
library("cowplot")
library("ggbreak")
library("ggforce")

# 1) Arguments and data
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)

# NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
# Cancers with no somatic VCFs
TCGA_cancers <- TCGA_cancers[!TCGA_cancers%in%c("TCGA-LAML","TCGA-MESO","TCGA-SKCM")]

# args <- commandArgs(trailingOnly=TRUE)
# print(paste0("ARGUMENTS USED --> ",args))
# pancancer <- args[1]
# NMD_method <- args[2]
# test <- args[3]
# VAF <- as.numeric(args[4])

pancancer <- "yes"
test <- "CGC"
NMD_method <- "ASE"
TCGA_cancer <- "pancancer"
VAF <- 0.2

if ( test == "CGC") {
    char <- "CGC_somatic_mut_CNV"
    char2 <- "somatic_associations"
} else if ( test == "subtypes") {
    char <- "subtypes"
    char2 <- "subtypes"
}
char3 <- ""

# ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
# # Remove versions from IDs
# ensembl_v88_gtf[,c("gene_id","transcript_id")] <- sapply(ensembl_v88_gtf[,c("gene_id","transcript_id")], function(col) { 
# gsub("(.*)\\..*","\\1", col)
# })
ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$type == "gene",]
ensembl_v88_gtf_filt$genome_location <- gsub("chr","",paste0(ensembl_v88_gtf_filt$seqnames,":",ensembl_v88_gtf_filt$start,"-",ensembl_v88_gtf_filt$end))
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[,c("gene_name","genome_location")]

# 1.2) Cancer types to be corrected (1) or not (0) by CNV arm-level PCs
cancer_types_corrected_CNV_PCs <- read.table(file = "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/cancer_types_CNV_PCs_correction.csv",
        sep = ",", header = TRUE)

# 2) QQplots and Volcano plots
# List of top hits per cancer type

# 2.1) Endogenous
# PCs == 1 for adding PCs to adjust. PCs == 0 without PCs. PCs == NULL to automatically select if PC is 1/0 depending on the cancer type
NMD_method <- "endogenous"
NMD_method_VAF <- "endogenous"
# Pancancer
all_CGC_assoc_end_pancancer <- CGC_associations(TCGA_cancers = "pancancer", NMD_method = NMD_method, num_PCs = "1",
                                    cancer_types_corrected_CNV_PCs = cancer_types_corrected_CNV_PCs) 
for (FDR in c(0.05,0.1,0.2,0.3,0.4,0.5)) {
    volcano_plot(all_CGC_assoc_end_pancancer, NMD_method_VAF, pval_adjust = "pval_adjusted", FDR = FDR, type = "pancancer", PCs = "1")
    volcano_plot(all_CGC_assoc_end_pancancer, NMD_method_VAF, pval_adjust = "qval", FDR = FDR, type = "pancancer", PCs = "1")
}
# By cancer type
all_CGC_assoc_end_bycancer <- CGC_associations(TCGA_cancers = TCGA_cancers, NMD_method = NMD_method, 
                                        num_PCs = NULL, cancer_types_corrected_CNV_PCs = cancer_types_corrected_CNV_PCs) 
for (pval_adjust in c("qval","pval_adjusted")) {
    print(pval_adjust)
    #for (FDR in c(0.01,0.05,0.1,0.25,0.5)) {
    for (FDR in c(0.01,0.05,0.1)) {
        print(FDR)
        volcano_plot(all_CGC_assoc_end_bycancer,NMD_method_VAF, pval_adjust = pval_adjust, FDR = FDR, type = "by_cancers", PCs = "1")
    }
}

# 2.2) ASE
NMD_method <- "ASE"
NMD_method_VAF <- paste0(NMD_method,"_",VAF)
# Pancancer
PCs <- "0"
all_CGC_assoc_ASE_pancancer <- CGC_associations(TCGA_cancers = "pancancer", NMD_method = NMD_method, num_PCs = PCs,
                                            cancer_types_corrected_CNV_PCs = cancer_types_corrected_CNV_PCs) 
for (FDR in c(0.05,0.1,0.2,0.3,0.4,0.5)) {
    volcano_plot(all_CGC_assoc_ASE_pancancer, NMD_method_VAF, pval_adjust = "pval_adjusted", FDR = FDR, type = "pancancer", PCs = PCs)
    volcano_plot(all_CGC_assoc_ASE_pancancer, NMD_method_VAF, pval_adjust = "qval", FDR = FDR, type = "pancancer", PCs = PCs)
}

# By cancer type
all_CGC_assoc_ASE_bycancer <- CGC_associations(TCGA_cancers = TCGA_cancers, NMD_method = NMD_method, num_PCs = NULL, 
                                         cancer_types_corrected_CNV_PCs = cancer_types_corrected_CNV_PCs) 
for (pval_adjust in c("qval","pval_adjusted")) {
    print(pval_adjust)
    for (FDR in c(0.01,0.05,0.1,0.25,0.5)) {
        print(FDR)
        volcano_plot(all_CGC_assoc_ASE_bycancer, NMD_method_VAF, pval_adjust = pval_adjust, FDR = FDR, type = "by_cancers", PCs = "1")
    }
}

# df <- all_CGC_assoc_ASE_pancancer$CNV_amp$pancancer
# df[df$Gene_symbol %in% "SMG5",]
# all_CGC_assoc_ASE_pancancer$CNV_amp$pancancer[all_CGC_assoc_ASE_pancancer$CNV_amp$pancancer$som_CNV_amp_pval_adjusted < 0.4,]
# all_CGC_assoc_ASE_pancancer$CNV_amp$pancancer[all_CGC_assoc_ASE_pancancer$CNV_amp$pancancer$Gene_symbol %in% c("SMG9","SMG1","SMG5","SMG6","SMG7","UPF1","UPF3B","UPF2","UPF3A"),]

# 3) Frequency and coefficient of replicated genes in ASE (discovery) + Endogenous (validation)
# Pancancer
for (pval_adjust in c("qval","pval_adjusted")) {
    print(pval_adjust)
    for (FDR in c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)) {
        print(FDR)
        freq_genes_plots(all_CGC_assoc_end_pancancer, all_CGC_assoc_ASE_pancancer, pval_adjust = pval_adjust, FDR = FDR, type = "pancancer")
    }
}
# All hits across cancers
for (pval_adjust in c("qval","pval_adjusted")) {
    print(pval_adjust)
    for (FDR in c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)) {
        print(FDR)
        freq_genes_plots(all_CGC_assoc_end_bycancer, all_CGC_assoc_ASE_bycancer, pval_adjust = pval_adjust, FDR = FDR, type = "by_cancers")
    }
}


########################################################################################

# Frequency
top_hits_df_filt$NMD_method <- factor(top_hits_df_filt$NMD_method)
df_counts <- top_hits_df_filt %>% group_by(NMD_method) %>% count(Gene_symbol)
replicated_genes <- df_counts$Gene_symbol[table(df_counts$Gene_symbol)>1]


subset(df_counts, Gene_symbol,replicated_genes)
df_counts_filt <- df_counts %>% filter(Gene_symbol %in% replicated_genes) %>% arrange(n)


# Take top 100
if (nrow(top_hits_df_freq) > 100) {
    top_hits_df_freq <- top_hits_df_freq[1:100,]
}

top_hits_df_freq <- top_hits_df_freq[order(top_hits_df_freq$Freq, decreasing = FALSE),]
top_hits_df_freq$Var3 <- paste0(top_hits_df_freq$Var1,"_",top_hits_df_freq$Var2)
top_hits_df_freq$Var3 <- factor(top_hits_df_freq$Var3, levels = top_hits_df_freq$Var3)
#top_hits_df_freq$Var1 <- factor(top_hits_df_freq$Var1, levels = top_hits_df_freq$Var1)

# NMD type
top_hits_df_freq <- merge(top_hits_df_freq,CGC, by.x = "Var1", by.y = "Gene.Symbol", all.x = TRUE)

barplot_freq_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/",NMD_method,"_barplot_freq_sig_hits.png")
freq_df_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/",NMD_method,"_freq_sig_hits.txt")

####
top_hits_df_freq <- read.table(file = freq_df_path, sep = "\t", header = TRUE)
top_hits_df_freq <- top_hits_df_freq[1:100,]
top_hits_df_freq <- top_hits_df_freq[order(top_hits_df_freq$frequency),]
top_hits_df_freq$Var3 <- paste0(top_hits_df_freq$gene_symbol,"_",top_hits_df_freq$coeff_direction)
top_hits_df_freq[which(is.na(top_hits_df_freq$NMD_type)),"Var3"] <- length(which(is.na(top_hits_df_freq$NMD_type))):1
top_hits_df_freq$Var3 <- factor(top_hits_df_freq$Var3, levels = top_hits_df_freq$Var3)
top_hits_df_freq$NMD_type <- as.character(top_hits_df_freq$NMD_type)
top_hits_df_freq$NMD_type <- ifelse(top_hits_df_freq$NMD_type != "NMD_factor","NMD_related",top_hits_df_freq$NMD_type)
top_hits_df_freq[which(is.na(top_hits_df_freq$NMD_type)),"NMD_type"] <- "TSG/OG"
####

png(barplot_freq_plot, width = 1500, height = 1500, res = 300)
p <- ggplot(data = top_hits_df_freq, aes(x = frequency, y = Var3, fill = NMD_type)) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_discrete(guide = guide_axis(check.overlap=TRUE)) +
    geom_bar(stat = "identity", position = "dodge") + #coord_flip(ylim = c(-1,1)) +
    ylab("") + xlab("Frequency") + ggtitle(paste0("Significant hits across cancer types")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text.x = element_text(color="black", size=20),
            axis.text.y = element_text(color="black", size=7)) +
    labs(color = "Direction of coefficient") + theme_classic()
print(p)
dev.off()

# Save df results
colnames(top_hits_df_freq_save) <- c("gene_symbol","coeff_direction","frequency")
top_hits_df_freq_save <- merge(top_hits_df_freq_save,CGC, by.x = "gene_symbol", by.y = "Gene.Symbol", all.x = TRUE)
top_hits_df_freq_save <- top_hits_df_freq_save[order(top_hits_df_freq_save$frequency, decreasing = TRUE),]

write.table(top_hits_df_freq_save, file = freq_df_path, sep = "\t", quote = FALSE,
                    col.names = TRUE, row.names = FALSE)






