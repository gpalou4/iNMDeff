rm(list=ls())

################################ FUNCTIONS ################################
################################ FUNCTIONS ################################
################################ FUNCTIONS ################################

negative_control_match <- function(NMD_geneset) {

    # Filter genes:
    # Use only NMD_control genes (they can have NMD features but at least they are MANE/Karousis)
    NMD_geneset_filt <- NMD_geneset %>% 
                        filter(!grepl(".*no stop.*|.*no start.*",comment)) %>%
                        filter(TCGA_pancancer_gene_exp >= 1 & GTEx_pantissue_gene_exp >= 1) %>%
                        filter(NMD_type == "NMD_control")
    NMD_geneset_ensembl_genes <- unique(NMD_geneset_filt$ensembl_gene_id)
    NMD_targets_controls <- data.frame(NMD_target = NA, NMD_control = NA)

    for (ensembl_gene in NMD_geneset_ensembl_genes) {
        NMD_nonNMD_pairs <- data.frame(NMD_target = NA, NMD_control = NA)
        NMD_geneset_gene <- NMD_geneset_filt %>%
                            filter(ensembl_gene_id %in% ensembl_gene)
        transcripts_num <- nrow(NMD_geneset_gene)
        # Keep genes with at least 1 NMD_target and 1 NMD_control
        if ( transcripts_num <= 1) {
            #print("Gene with no 1-1 NMD_target-NMD_control...")
            next
        } else {
            random_transcripts <- sample(transcripts_num)
            NMD_target <- NMD_geneset_gene[random_transcripts[1],"ensembl_transcript_id"]
            NMD_control <- NMD_geneset_gene[random_transcripts[2],"ensembl_transcript_id"]
            NMD_targets_controls <- rbind(NMD_targets_controls,c(NMD_target,NMD_control))
        }
    }
    NMD_targets_controls <- na.omit(NMD_targets_controls)
    return(NMD_targets_controls)
}

NMD_target_control_combinations <- function(NMD_geneset) {

    NMD_geneset_ensembl_genes <- unique(NMD_geneset$ensembl_gene_id)

    NMD_combinations_list <- list()
    for (ensembl_gene in NMD_geneset_ensembl_genes) {

        NMD_nonNMD_pairs <- data.frame(NMD_target = NA, NMD_control = NA)
        NMD_geneset_gene <- NMD_geneset %>%
                            filter(ensembl_gene_id %in% ensembl_gene)
        # Keep genes with at least 1 NMD_target and 1 NMD_control
        if (length(unique(NMD_geneset_gene$NMD_type)) == 1) {
            #print("Gene with no 1-1 NMD_target-NMD_control...")
            next
        }
        # For each NMD_target-NMD_control combination
        NMD_targets <- NMD_geneset_gene$ensembl_transcript_id[NMD_geneset_gene$NMD_type == "NMD_target"]
        NMD_controls <- NMD_geneset_gene$ensembl_transcript_id[NMD_geneset_gene$NMD_type == "NMD_control"]
        NMD_combinations <- expand.grid(NMD_target = NMD_targets, NMD_control = NMD_controls)
        NMD_combinations[,c("ensembl_gene_id","NMD_event","CDS_diff","TSS_diff","TES_diff","Karousis_NMD_target","Karousis_NMD_control", 
                        "NMD_target_stop_codon_to_3UTR_splice_site","NMD_control_stop_codon_to_3UTR_splice_site", "NMD_target_CNV_PCs_adj_r_squared", "NMD_control_CNV_PCs_adj_r_squared",
                        "MANE_control","CL_UPF1_KD_NMD_nonNMD_ratio_median","CL_WT_NMD_nonNMD_ratio_median","CL_NMD_nonNMD_ratio_median_diff",
                        "transcript_type_NMD_target","transcript_type_NMD_control","TSL_NMD_control","TSL_NMD_target","NMD_nonNMD_ratio_TCGA", "NMD_nonNMD_ratio_GTEx",
                        "NMD_control_TCGA_pancancer_gene_exp", "NMD_target_TCGA_pancancer_gene_exp","NMD_control_GTEx_pantissue_gene_exp", "NMD_target_GTEx_pantissue_gene_exp",
                        "NMD_target_DGE","NMD_control_DGE","NMD_target_PC2_outliers","NMD_control_PC2_outliers", "NMD_control_NMD_Colombo", "NMD_target_NMD_Colombo",
                        "NMD_target_purity_corr","NMD_control_purity_corr","NMD_target_CNV_burden_corr", "NMD_control_CNV_burden_corr",
                        "GSE152435_KD1_ratio","GSE86148_KD2_ratio","GSE88466_KD3_ratio","GSE88140_KD4_ratio",
                        "GSE152435_WT1_ratio","GSE86148_WT2_ratio","GSE88148_WT3_ratio","GSE88266_WT4_ratio",
                        "GSE152435_GSE152435_KD_WT_OR","GSE86148_GSE86148_KD_WT_OR","GSE88466_GSE88148_KD_WT_OR","GSE88140_GSE88266_KD_WT_OR",
                        "GSE152435_GSE152435_log2FC_NMD_target","GSE86148_GSE86148_log2FC_NMD_target","GSE88466_GSE88148_log2FC_NMD_target","GSE88140_GSE88266_log2FC_NMD_target",
                        "GSE152435_GSE152435_log2FC_NMD_control","GSE86148_GSE86148_log2FC_NMD_control","GSE88466_GSE88148_log2FC_NMD_control","GSE88140_GSE88266_log2FC_NMD_control")] <- NA
        # Obtain info
        for (i in 1:nrow(NMD_combinations)) {
            NMD_combination <- NMD_combinations[i,]
            NMD_target <- NMD_geneset_gene %>%
                                filter(ensembl_transcript_id == NMD_combination$NMD_target)
            NMD_control <- NMD_geneset_gene %>%
                                filter(ensembl_transcript_id == NMD_combination$NMD_control)
            NMD_combinations[i,"ensembl_gene_id"] <- unique(NMD_target$ensembl_gene_id)
            # Karousis, MANE, NMD event type
            NMD_combinations[i,"NMD_event"] <- NMD_target$NMD_event_type
            NMD_combinations[i,"NMD_target_uORFs"] <- NMD_target$uORFs
            NMD_combinations[i,"NMD_control_uORFs"] <- NMD_control$uORFs
            NMD_combinations[i,"NMD_target_UTR3_GC_content"] <- NMD_target$UTR3_GC_content
            NMD_combinations[i,"NMD_control_UTR3_GC_content"] <- NMD_control$UTR3_GC_content
            NMD_combinations[i,"Karousis_NMD_target"] <- NMD_target$Karousis_paper
            NMD_combinations[i,"Karousis_NMD_control"] <- NMD_control$Karousis_paper
            NMD_combinations[i,"MANE_control"] <- NMD_control$MANE_Select
            NMD_combinations[i,"transcript_type_NMD_target"] <- NMD_target$transcript_type
            NMD_combinations[i,"transcript_type_NMD_control"] <- NMD_control$transcript_type
            NMD_combinations[i,"TSL_NMD_target"] <- NMD_target$transcript_support_level
            NMD_combinations[i,"TSL_NMD_control"] <- NMD_control$transcript_support_level
            NMD_combinations[i,"NMD_control_TCGA_pancancer_gene_exp"] <- NMD_control$TCGA_pancancer_gene_exp
            NMD_combinations[i,"NMD_target_TCGA_pancancer_gene_exp"] <- NMD_target$TCGA_pancancer_gene_exp
            NMD_combinations[i,"NMD_control_GTEx_pantissue_gene_exp"] <- NMD_control$GTEx_pantissue_gene_exp
            NMD_combinations[i,"NMD_target_GTEx_pantissue_gene_exp"] <- NMD_target$GTEx_pantissue_gene_exp
            NMD_combinations[i,"NMD_target_stop_codon_to_3UTR_splice_site"] <- NMD_target$stop_codon_to_3UTR_splice_site
            NMD_combinations[i,"NMD_control_stop_codon_to_3UTR_splice_site"] <- NMD_control$stop_codon_to_3UTR_splice_site
            NMD_combinations[i,"NMD_target_CNV_PCs_adj_r_squared"] <- NMD_target$CNV_PCs_adj_r_squared
            NMD_combinations[i,"NMD_control_CNV_PCs_adj_r_squared"] <- NMD_control$CNV_PCs_adj_r_squared
            NMD_combinations[i,"NMD_target_DGE"] <- NMD_target$DGE_transcript
            NMD_combinations[i,"NMD_control_DGE"] <- NMD_control$DGE_transcript
            NMD_combinations[i,"NMD_target_PC2_outliers"] <- NMD_target$PC2_outliers
            NMD_combinations[i,"NMD_control_PC2_outliers"] <- NMD_control$PC2_outliers
            NMD_combinations[i,"NMD_target_NMD_Colombo"] <- NMD_target$NMD_Colombo
            NMD_combinations[i,"NMD_control_NMD_Colombo"] <- NMD_control$NMD_Colombo
            NMD_combinations[i,"NMD_target_purity_corr"] <- NMD_target$purity_corr
            NMD_combinations[i,"NMD_control_purity_corr"] <- NMD_control$purity_corr
            NMD_combinations[i,"NMD_target_CNV_burden_corr"] <- NMD_target$CNV_burden_corr
            NMD_combinations[i,"NMD_control_CNV_burden_corr"] <- NMD_control$CNV_burden_corr
            # CDS/TSS/TES distance
            CDS_diff <- abs(as.numeric(NMD_control$CDS) - as.numeric(NMD_target$CDS))
            NMD_combinations[i,"CDS_diff"] <- CDS_diff
            TSS_diff <- abs(as.numeric(NMD_control$transcription_start_site) - as.numeric(NMD_target$transcription_start_site))
            NMD_combinations[i,"TSS_diff"] <- TSS_diff
            TES_diff <- abs(as.numeric(NMD_control$transcription_end_site) - as.numeric(NMD_target$transcription_end_site))
            NMD_combinations[i,"TES_diff"] <- TES_diff
            # NMD/nonMMD ratio in TCGA and GTEx
            NMD_combinations[i,"NMD_nonNMD_ratio_TCGA"] <- NMD_target$TCGA_pancancer_gene_exp / NMD_control$TCGA_pancancer_gene_exp
            NMD_combinations[i,"NMD_nonNMD_ratio_GTEx"] <- NMD_target$GTEx_pantissue_gene_exp / NMD_control$GTEx_pantissue_gene_exp
            # NMD/nonNMD ratio in CL UPF1 KD
            CL_UPF1KD_NMD_filt <- CL_UPF1KD_NMD[rownames(CL_UPF1KD_NMD) %in% c(NMD_target$ensembl_transcript_id,NMD_control$ensembl_transcript_id),]
            CL_UPF1KD_NMD_nonNMD_ratios <- CL_UPF1KD_NMD_filt[NMD_target$ensembl_transcript_id,] / CL_UPF1KD_NMD_filt[NMD_control$ensembl_transcript_id,]
            colnames(CL_UPF1KD_NMD_nonNMD_ratios) <- paste0(colnames(CL_UPF1KD_NMD_nonNMD_ratios),"_ratio")
            #CL_UPF1KD_NMD_nonNMD_ratios <- ifelse(CL_UPF1KD_NMD_nonNMD_ratios %in% c("Inf","-Inf"), 0, CL_UPF1KD_NMD_nonNMD_ratios)
            CL_UPF1_KD_NMD_nonNMD_ratio_median <- median(as.numeric(CL_UPF1KD_NMD_nonNMD_ratios), na.rm = TRUE)
            NMD_combinations[i,"CL_UPF1_KD_NMD_nonNMD_ratio_median"] <- CL_UPF1_KD_NMD_nonNMD_ratio_median
            NMD_combinations[i,colnames(CL_UPF1KD_NMD_nonNMD_ratios)] <- CL_UPF1KD_NMD_nonNMD_ratios
            # NMD/nonNMD ratio in CL WT
            CL_WT_NMD_filt <- CL_UPF1KD_controls[rownames(CL_UPF1KD_controls) %in% c(NMD_target$ensembl_transcript_id,NMD_control$ensembl_transcript_id),]
            CL_WT_NMD_nonNMD_ratios <- CL_WT_NMD_filt[NMD_target$ensembl_transcript_id,] / CL_WT_NMD_filt[NMD_control$ensembl_transcript_id,]
            colnames(CL_WT_NMD_nonNMD_ratios) <- paste0(colnames(CL_WT_NMD_nonNMD_ratios),"_ratio")
            # CL_WT_NMD_nonNMD_ratios <- ifelse(CL_WT_NMD_nonNMD_ratios %in% c("Inf","-Inf"), 0, CL_WT_NMD_nonNMD_ratios)
            CL_WT_NMD_nonNMD_ratio_median <- median(as.numeric(CL_WT_NMD_nonNMD_ratios), na.rm = TRUE)
            NMD_combinations[i,"CL_WT_NMD_nonNMD_ratio_median"] <- CL_WT_NMD_nonNMD_ratio_median
            NMD_combinations[i,colnames(CL_WT_NMD_nonNMD_ratios)] <- CL_WT_NMD_nonNMD_ratios
            # NMD/nonNMD ratio diff
            NMD_combinations[i,"CL_NMD_nonNMD_ratio_median_diff"] <- CL_WT_NMD_nonNMD_ratio_median - CL_UPF1_KD_NMD_nonNMD_ratio_median  
            CL_NMD_nonNMD_ratio_diff <- CL_WT_NMD_nonNMD_ratios - CL_UPF1KD_NMD_nonNMD_ratios
            # colnames(CL_NMD_nonNMD_ratio_diff) <- paste0(colnames(CL_NMD_nonNMD_ratio_diff))
            # threshold <- -0.1
            # NMD_combinations[i,"num_CL_ratios_above_threshold"] <- sum(CL_NMD_nonNMD_ratio_diff <= threshold)
            # NMD_combinations[i,"num_CL_ratios_same_direction"] <- sum(CL_NMD_nonNMD_ratio_diff <= 0)
            # 1) NMD_target/NMD_control ratio in KD vs WT
            CL_NMD_OR_ratios <- log2(CL_UPF1KD_NMD_nonNMD_ratios / CL_WT_NMD_nonNMD_ratios)
            NMD_combinations[i,grep("KD_WT_OR",colnames(NMD_combinations))] <- as.numeric(CL_NMD_OR_ratios)
            OR_threshold <- 0.5
            NMD_combinations[i,"num_CL_OR_above_threshold"] <- sum(CL_NMD_OR_ratios >= OR_threshold)
            NMD_combinations[i,"num_CL_OR_same_direction"] <- sum(CL_NMD_OR_ratios >= 0)
            # 2) Log2FC between KD and WT for NMD_target and NMD_control separately
            CL_UPF1KD_log2FC <- log2(CL_UPF1KD_NMD_filt[NMD_target$ensembl_transcript_id,]) - log2(CL_WT_NMD_filt[NMD_target$ensembl_transcript_id,])
            NMD_combinations[i,grep("log2FC_NMD_target",colnames(NMD_combinations))] <- as.numeric(CL_UPF1KD_log2FC)
            CL_WT_log2FC <- log2(CL_UPF1KD_NMD_filt[NMD_control$ensembl_transcript_id,]) - log2(CL_WT_NMD_filt[NMD_control$ensembl_transcript_id,])
            NMD_combinations[i,grep("log2FC_NMD_control",colnames(NMD_combinations))] <- as.numeric(CL_WT_log2FC)
            log2FC_threshold <- 0.5
            NMD_combinations[i,"num_CL_KD_log2FC_above_threshold"] <- sum(CL_UPF1KD_log2FC >= log2FC_threshold)
            NMD_combinations[i,"num_CL_WT_log2FC_above_threshold"] <- sum(CL_WT_log2FC >= log2FC_threshold)
            NMD_combinations[i,"num_CL_KD_log2FC_same_direction"] <- sum(CL_UPF1KD_log2FC >= 0)
            NMD_combinations[i,"num_CL_WT_log2FC_same_direction"] <- sum(CL_WT_log2FC >= 0)
            # Without Weird Cell Line
            NMD_combinations[i,"num_CL_KD_log2FC_above_threshold_outlier"] <- sum(CL_UPF1KD_log2FC[,-1] >= log2FC_threshold)
            NMD_combinations[i,"num_CL_WT_log2FC_above_threshold_outlier"] <- sum(CL_WT_log2FC[,-1] >= log2FC_threshold)
            NMD_combinations[i,"num_CL_KD_log2FC_same_direction_outlier"] <- sum(CL_UPF1KD_log2FC[,-1] >= 0)
            NMD_combinations[i,"num_CL_WT_log2FC_same_direction_outlier"] <- sum(CL_WT_log2FC[,-1] >= 0)
        }

        NMD_combinations_filt <- NMD_combinations
        if ( nrow(NMD_combinations_filt) == 0) {
            next
        }
        # 2) Order transcripts as follows:
        # Karousis NMD target > MANE control > Karousis Control > NMD_events (2 features > 3EJC > uORF) >
        # UTR3 GC NMD_target > number of uORFs NMD_target (if possible) > gene expression of NMD_control 
        NMD_events <- c(" | >=2 uORFs | Intronic-spliced 3UTR"," | Intronic-spliced 3UTR"," | >=2 uORFs")
        NMD_combinations_filt <- arrange(NMD_combinations_filt,desc(Karousis_NMD_target),desc(Karousis_NMD_control),desc(MANE_control),match(NMD_event, NMD_events),desc(NMD_target_UTR3_GC_content),desc(NMD_target_uORFs),desc(NMD_control_TCGA_pancancer_gene_exp),desc(NMD_control_GTEx_pantissue_gene_exp))
        # Save
        NMD_combinations_list[[ensembl_gene]] <- NMD_combinations_filt
    }
    return(NMD_combinations_list)
}

NMD_target_control_final_matching <- function(NMD_combinations_list, dataset) {
    NMD_combinations_list_filt <- lapply(NMD_combinations_list, function( gene ) {
        if (dataset == "TCGA") {
            gene_filt <- gene %>%
                        filter( ( NMD_event %in% c(" | Intronic-spliced 3UTR"," | >=2 uORFs | Intronic-spliced 3UTR") 
                                & NMD_target_stop_codon_to_3UTR_splice_site >= 50) | ( NMD_event == " | >=2 uORFs" ) ) %>%        
                        filter(NMD_nonNMD_ratio_TCGA < 0.9) %>%
                        filter(NMD_target_TCGA_pancancer_gene_exp >= 1 & NMD_control_TCGA_pancancer_gene_exp >= 3)
                        # filter(NMD_nonNMD_ratio_TCGA < 1) %>%
                        # filter(NMD_target_TCGA_pancancer_gene_exp >= 1 & NMD_control_TCGA_pancancer_gene_exp >= 1) %>%
                        # filter(NMD_target_TCGA_pancancer_gene_exp + NMD_control_TCGA_pancancer_gene_exp >= 5) %>%
                        # filter(! (num_CL_KD_log2FC_above_threshold_outlier <=1 & NMD_target_PC2_outliers == "yes" & TSL_NMD_target == 5)) %>%
                        # filter( ( (abs(NMD_target_CNV_burden_corr) < 0.25) | ( abs(NMD_control_CNV_burden_corr) < 0.25 ) ) ) %>%
                        # filter( ( (NMD_target_purity_corr > -0.10) & (NMD_control_purity_corr > -0.10 ) ) )
        } else if (dataset == "GTEx") {
            gene_filt <- gene %>%
                        filter( ( NMD_event %in% c(" | Intronic-spliced 3UTR"," | >=2 uORFs | Intronic-spliced 3UTR") 
                                & NMD_target_stop_codon_to_3UTR_splice_site >= 50) | ( NMD_event == " | >=2 uORFs" ) ) %>%        
                        filter(NMD_nonNMD_ratio_GTEx < 0.9) %>%
                        filter(NMD_target_GTEx_pantissue_gene_exp >= 0.2 & NMD_control_GTEx_pantissue_gene_exp >= 2)
                        #filter(NMD_target_GTEx_pantissue_gene_exp + NMD_control_GTEx_pantissue_gene_exp >= 5) %>%
                        # filter(! (num_CL_KD_log2FC_above_threshold_outlier <=1 & NMD_target_PC2_outliers == "yes" & TSL_NMD_target == 5)) %>%
                        # filter( ( (abs(NMD_target_CNV_burden_corr) < 0.25) | ( abs(NMD_control_CNV_burden_corr) < 0.25 ) ) ) %>%
                        # filter( ( (NMD_target_purity_corr > -0.10) & (NMD_control_purity_corr > -0.10 ) ) )
        }
        if (nrow(gene_filt) != 0 ) {
            #gene_filt <- arrange(gene_filt,desc(num_CL_KD_log2FC_above_threshold), num_CL_WT_log2FC_above_threshold)
            gene_filt[1,]
        } else {
            NULL
        }
    })

    NMD_combinations_df <- do.call(rbind,NMD_combinations_list_filt)
    # Add NMD_type_final label (final NMD_target-NMD_control match)
    NMD_geneset_final <- NMD_geneset
    NMD_geneset_final[NMD_geneset_final$ensembl_transcript_id %in% NMD_combinations_df$NMD_target,"NMD_type_final"] <- "NMD_target"
    NMD_geneset_final[NMD_geneset_final$ensembl_transcript_id %in% NMD_combinations_df$NMD_control,"NMD_type_final"] <- "NMD_control"
    NMD_geneset_final <- NMD_geneset_final[!is.na(NMD_geneset_final$NMD_type_final),]
    NMD_combinations_df$ensembl_gene_id <- NULL
    # Add gene info
    cols <- c("NMD_control","NMD_event","Karousis_NMD_target","Karousis_NMD_control","MANE_control","ensembl_gene_id",
                "TSL_NMD_control","TSL_NMD_target","transcript_type_NMD_target","transcript_type_NMD_control")
    NMD_geneset_final <- merge(NMD_geneset_final,NMD_combinations_df[,!colnames(NMD_combinations_df) %in% cols], by.x = "ensembl_transcript_id", by.y = "NMD_target", all.x = TRUE)
    # Add CL gene expression from KD and WT
    CL_UPF1KD_NMD_filt <- CL_UPF1KD_NMD[rownames(CL_UPF1KD_NMD) %in% c(as.character(NMD_combinations_df$NMD_target),as.character(NMD_combinations_df$NMD_control)),]
    CL_WT_NMD_filt <- CL_UPF1KD_controls[rownames(CL_UPF1KD_controls) %in% c(as.character(NMD_combinations_df$NMD_target),as.character(NMD_combinations_df$NMD_control)),]
    NMD_geneset_final <- merge(NMD_geneset_final,CL_UPF1KD_NMD_filt, by.x = "ensembl_transcript_id", by.y = "row.names", all.x = TRUE)
    NMD_geneset_final <- merge(NMD_geneset_final,CL_WT_NMD_filt, by.x = "ensembl_transcript_id", by.y = "row.names", all.x = TRUE)
    return(NMD_geneset_final)
}

median_quartile <- function(x) {
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
}

check_NMD_features <- function(ensembl_NMD_features, set = NULL, control, gene_exp_NMD_target = 1, gene_exp_control = 3) {

  if (!is.null(set)) {
    ensembl_NMD_features_set <- ensembl_NMD_features[ensembl_NMD_features$ensembl_transcript_id%in%set,]
  } else {
    ensembl_NMD_features_set <- ensembl_NMD_features
  }
  if (control == "yes") {
    # no Start or stop codons
    ensembl_NMD_features_set <- ensembl_NMD_features_set[grep(".*no stop.*|.*no start.*",ensembl_NMD_features_set$comment, invert = TRUE),]
    filt <- (ensembl_NMD_features_set$NMD_event_type == "") & (ensembl_NMD_features_set$TCGA_pancancer_gene_exp >= gene_exp_control) & (ensembl_NMD_features_set$uORFs==0 | ensembl_NMD_features_set$uORFs == "")
    controls.index <- which(filt)
    # If no controls
    if (length(controls.index) == 0 ) {
      #print("Gene with no control(s) or low expressed (<= 3), skipping")
      return("NA")
    }
    ensembl_NMD_features_set_filt <- ensembl_NMD_features_set[controls.index,]
    # Choose the control with highest expression
    ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt[which(ensembl_NMD_features_set_filt$TCGA_pancancer_gene_exp == max(ensembl_NMD_features_set_filt$TCGA_pancancer_gene_exp)),]
    # If more than 1 transcript is chosen (same gene expression), then take one random
    ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_gene_exp[sample(1:nrow(ensembl_NMD_features_set_gene_exp))[1],]
    control.transcript <- ensembl_NMD_features_set_gene_exp$ensembl_transcript_id
    return(control.transcript)
  } else if (control == "no") {
    # no Start or stop codons
    ensembl_NMD_features_set <- ensembl_NMD_features_set[grep(".*no stop.*|.*no start.*",ensembl_NMD_features_set$comment, invert = TRUE),]
    NMD.targets.index <- which( (ensembl_NMD_features_set$NMD_event_type != "") & (ensembl_NMD_features_set$TCGA_pancancer_gene_exp >= gene_exp_NMD_target) )
    # If there aren't NMD features or low expression, skip
    if (length(NMD.targets.index) == 0 ) {
      #print("Gene with no NMD target(s) or low expressed (< 1 TPM), skipping")
      return("NA")
  }
  # Criteria
  # 1) 3UTR EJC 50nt > uORF 
  # 2) highest 3'UTR GC content
  ensembl_NMD_features_set_filt <- ensembl_NMD_features_set[NMD.targets.index ,]
  ensembl_NMD_features_set_filt_3UTR_EJC <- ensembl_NMD_features_set_filt[grep("Intronic",ensembl_NMD_features_set_filt$NMD_event_type),]
  ensembl_NMD_features_set_filt_uORFs <- ensembl_NMD_features_set_filt[ensembl_NMD_features_set_filt$NMD_event_type == " | >=2 uORFs",]
  # 50nt
  ensembl_NMD_features_set_filt_3UTR_EJC <- ensembl_NMD_features_set_filt_3UTR_EJC[which(as.numeric(ensembl_NMD_features_set_filt_3UTR_EJC$stop_codon_to_3UTR_splice_site) >= 50),]

  if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) == 0 && nrow(ensembl_NMD_features_set_filt_uORFs) == 0) {
    return("NA")
  } else if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) >= 1 && nrow(ensembl_NMD_features_set_filt_uORFs) >=0 ) { # 3UTR EJC
    # Take the isoform with highest 3'UTR GC content
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_filt_3UTR_EJC[ensembl_NMD_features_set_filt_3UTR_EJC$UTR3_GC_content == max(ensembl_NMD_features_set_filt_3UTR_EJC$UTR3_GC_content),]
    #ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt_3UTR_EJC[which(ensembl_NMD_features_set_filt_3UTR_EJC$TCGA_pancancer_gene_exp == min(ensembl_NMD_features_set_filt_3UTR_EJC$TCGA_pancancer_gene_exp)),]
  } else if (nrow(ensembl_NMD_features_set_filt_3UTR_EJC) == 0 && nrow(ensembl_NMD_features_set_filt_uORFs) >= 1 ) { # uORFs
    # Take the isoform with highest 3'UTR GC content
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_filt_uORFs[ensembl_NMD_features_set_filt_uORFs$UTR3_GC_content == max(ensembl_NMD_features_set_filt_uORFs$UTR3_GC_content),]
    #ensembl_NMD_features_set_gene_exp <- ensembl_NMD_features_set_filt_uORFs[which(ensembl_NMD_features_set_filt_uORFs$TCGA_pancancer_gene_exp == min(ensembl_NMD_features_set_filt_uORFs$TCGA_pancancer_gene_exp)),]
  }
  # If more than 1 transcript is kept (same 3'UTR GC content)
  # then take the one with highest number of uORFs
  if (nrow(ensembl_NMD_features_set_3utr_gc) >1) {
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_3utr_gc[ensembl_NMD_features_set_3utr_gc$uORFs == max(ensembl_NMD_features_set_3utr_gc$uORFs),]
  }
  # If more than 1 transcript is kept (same #uORFs)
  # then take the one random
  if (nrow(ensembl_NMD_features_set_3utr_gc) >1) {
    ensembl_NMD_features_set_3utr_gc <- ensembl_NMD_features_set_3utr_gc[sample(1:nrow(ensembl_NMD_features_set_3utr_gc))[1],]
  }
  NMD.target <- ensembl_NMD_features_set_3utr_gc$ensembl_transcript_id
  return(NMD.target)
  }
}

NMD_geneset_info_features <- function(ensembl_NMD_features, NMD_geneset, all_NMD_genes = all.NMD.genes, non_NMD_genes = non.NMD.genes) {

  # Keep info
  n <- nrow(ensembl_NMD_features)
  n.uORFs <- as.numeric(table((as.numeric(ensembl_NMD_features$uORFs) > 0))[2])
  n.genes <- length(unique(ensembl_NMD_features$ensembl_gene_id))
  NMD.features.df <- data.frame(NMD_features = c(NMD.features,"fisher_pval","fisher_OR"), value = NA, row.names = 1)
  colnames(NMD.features.df)[1] <- NMD_geneset

  NMD.features.df["transcripts_size",] <- n
  NMD.features.df["genes_size",] <- n.genes
  NMD.features.df["uORFs",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs" & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["uORFs_1_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_1_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["uORFs_2_30_bp",] <- round(sum(ensembl_NMD_features$uORFs_2_30_bp & ensembl_NMD_features$UTR3 <= 1) / n,2)
  NMD.features.df["UTR3_EJC",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR") / n,2)
  NMD.features.df["UTR3_EJC_50nt",] <- round(length(which(ensembl_NMD_features$stop_codon_to_3UTR_splice_site >= 50 & ensembl_NMD_features$NMD_event_type == " | Intronic-spliced 3UTR")) / n,2)
  #NMD.features.df["uORFs_UTR3_intron",] <- round(sum(ensembl_NMD_features$NMD_event_type == " | >=2 uORFs | Intronic-spliced 3UTR") / n,2)
  NMD.features.df["non_NMD",] <- round(sum(ensembl_NMD_features$NMD_event_type == "") / n,2)
  NMD.features.df["UTR3_GC_content",] <- median(as.numeric(ensembl_NMD_features$UTR3_GC_content), na.rm = TRUE)
  NMD.features.df["uORFs_num_translated",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_translated), na.rm = TRUE),2)
  NMD.features.df["uORFs_num_trans_reinit",] <- round(mean(as.numeric(ensembl_NMD_features$uORFs_num_trans_reinit), na.rm = TRUE),2)
  NMD.features.df["penultime_codon_GCG",] <- round(length(grep("GCG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Ala
  #NMD_geneset_names.res[filt,"penultime_codon_GTG_control"] <- round(length(grep("GTG",NMD_geneset$uORF_penultimate_codon)) / n,6) # Val
  NMD.features.df["penultime_codon_ACA_control",] <- round(length(grep("ACA",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Thr
  NMD.features.df["penultime_codon_AGG",] <- round(length(grep("AGG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Arg
  #NMD_geneset_names.res[filt,"penultime_codon_AAG_control"] <- round(length(grep("AAG",NMD_geneset$uORF_penultimate_codon)) / n,6) # Lys
  NMD.features.df["penultime_codon_CCA_control",] <- round(length(grep("CCA",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Pro
  NMD.features.df["penultime_codon_CTG",] <- round(length(grep("CTG",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Leu
  #NMD_geneset_names.res[filt,"penultime_codon_CAT_control"] <- round(length(grep("CAT",NMD_geneset$uORF_penultimate_codon)) / n,6) # His
  NMD.features.df["penultime_codon_TAC_control",] <- round(length(grep("TAC",ensembl_NMD_features$uORF_penultimate_codon)) / n.uORFs,6) # Tyr
  
  # Fisher test
  ensembl.NMD.features[,"NMD_geneset"] <- NA
  ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%eval(parse(text=paste0(NMD_geneset))),"NMD_geneset"] <- NMD_geneset
  if (NMD_geneset == "non.NMD.genes") {
    ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%all_NMD_genes,"NMD_geneset"] <- "controls"
  } else {
    ensembl.NMD.features[ensembl.NMD.features$ensembl_gene_id%in%non_NMD_genes,"NMD_geneset"] <- "controls"
  }

  print(table(ensembl.NMD.features$NMD_geneset))
  # NMD.feature <- names(table(ensembl.NMD.features$NMD_event_type))[4]
  NMD.noNMDfeat <- length(which(ensembl.NMD.features$NMD_event_type=="" & ensembl.NMD.features$NMD_geneset == NMD_geneset))
  NMD.NMDfeat <- length(which(ensembl.NMD.features$NMD_event_type!="" & ensembl.NMD.features$NMD_geneset == NMD_geneset))
  noNMD.noNMDfeat <- length(which(ensembl.NMD.features$NMD_event_type=="" & ensembl.NMD.features$NMD_geneset == "controls"))
  noNMD.NMDfeat <- length(which(ensembl.NMD.features$NMD_event_type!="" & ensembl.NMD.features$NMD_geneset == "controls"))
  fisher.table <- data.frame(NMD_feature = c(NMD.NMDfeat,noNMD.NMDfeat), nonNMD_feature = c(NMD.noNMDfeat,noNMD.noNMDfeat))
  rownames(fisher.table) <- c(NMD_geneset,"controls")
  fisher.res <- fisher.test(fisher.table)
  print(fisher.table)
  print(fisher.res)
  # chisq.test(fisher.table)$expected
  NMD.features.df["fisher_pval",]  <- fisher.res$p.value
  NMD.features.df["fisher_OR",]  <- fisher.res$estimate
  return(t(NMD.features.df))

}

karousis_test <- function(ensembl.gene, ensembl.NMD.features.ensembl.gene) {
  # Check NMD features in all the transcript of the given gene and select the best NMD target
  NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
  # Check NMD features in all the transcript of the given gene and select the best nonNMD control one
  nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
  # Add our criteria to Karousis
  row <- which(Karousis_original_dataset$gene_id%in%ensembl.gene)
  Karousis_original_dataset[row,"NMD_target"] <<- NMD.target
  Karousis_original_dataset[row,"control"] <<- nonNMD.control
}


################################ LIBRARIES ################################
################################ LIBRARIES ################################
################################ LIBRARIES ################################

library("dplyr")
library("factoextra")
library("FactoMineR")
library("matrixStats")
library("tidyr")

################################ SCRIPT ################################
################################ SCRIPT ################################
################################ SCRIPT ################################

# Matching of NMD target and NMD control for each gene
# Criteria for matching NMD-target / NMD-control

# 1) Data

# 1.1) Associations between CNV and NMD transcripts expression
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/TCGA_NMD_transcripts_CNV_PCs_correlation.txt")
TCGA_NMD_transcripts_CNV_PCs_correlation <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1.2) DGE NMD genes in Cell Lines
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/DGE/top_hits_overlap_cell_lines.txt")
DGE_CL_top_hits <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(DGE_CL_top_hits)[1] <- "ensembl_transcript_id"
DGE_CL_top_hits$DGE_transcript <- "yes"

# 1.3) PC2 outliers
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/plots/PC2_outliers_TCGA_NMD_geneset_test_15.txt")
PC2_transcript_outliers <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(PC2_transcript_outliers)[1] <- "ensembl_transcript_id"
PC2_transcript_outliers$PC2_outliers <- "yes"

# # 1.4) NMD_Colombo old genes (NMDeff associated to CNV, R ~ 0.3)
# input_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/old/NMD_Colombo_ensembl_filt.txt"
# NMD_Colombo_old <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# NMD_Colombo_old <- NMD_Colombo_old[,c("ensembl_transcript_id"), drop = FALSE]
# NMD_Colombo_old$NMD_Colombo <- "yes"

# 1.4) NMD_Colombo all genes (old NMDeff geneset associated to CNV, R ~ 0.3)
input_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_Colombo_ensembl.txt"
NMD_Colombo_all <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMD_Colombo_all <- NMD_Colombo_all[,c("ensembl_transcript_id"), drop = FALSE]
NMD_Colombo_all$NMD_Colombo <- "yes"     

# 1.5) Transcripts associated with purity
input_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_pancancer_transcripts_purity_correlation.txt"
transcripts_purity_corr <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(transcripts_purity_corr) <- "purity_corr"

# 1.6) Transcripts associated with CNV burden
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_pancancer_transcripts_CNV_burden_correlation.txt")
CNV_burden_corr <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(CNV_burden_corr) <- "CNV_burden_corr"

# 1.7) Karousis original NMD geneset
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/raw/Karousis_ensembl_gene_symbol.txt")
Karousis_original_dataset <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Karousis_original_dataset[,c("NMD_target","control")] <- NA

# 1.7) NMD targets genesets
NMD_genesets <- c("NMD_Courtney","NMD_Karousis","NMD_Colombo","NMD_Tani","non_NMD_genes",
                "NMD_global","NMD_global_all_shared","NMD_global_2_shared","NMD_ensembl","SMG6","SMG7")
NMD_targets_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts"

for (NMD_geneset_name in NMD_genesets) {
    if (NMD_geneset_name %in% "non_NMD_genes") {
        char <- "ensembl.txt"
    } else {
        char <- "ensembl_filt.txt"
    }
    NMD_geneset <- read.table(file = paste0(NMD_targets_path,"/",NMD_geneset_name,"_",char), 
                                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    remove <- which(NMD_geneset$tag == "PAR")
    if (length(remove != 0)) {
        NMD_geneset <- NMD_geneset[-remove,] 
    }
    # Add CNV associations
    NMD_geneset <- merge(NMD_geneset,TCGA_NMD_transcripts_CNV_PCs_correlation, by.x = "ensembl_transcript_id", by.y = "NMD_transcript", all.x = TRUE)
    # Add DGE genes
    NMD_geneset <- merge(NMD_geneset,DGE_CL_top_hits, by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all.x = TRUE)
    NMD_geneset$DGE_transcript <- ifelse(is.na(NMD_geneset$DGE_transcript),"no","yes")
    # Add PC2 outliers
    NMD_geneset <- merge(NMD_geneset,PC2_transcript_outliers, by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all.x = TRUE)
    NMD_geneset$PC2_outliers <- ifelse(is.na(NMD_geneset$PC2_outliers),"no","yes")
    # Add Colombo genes
    NMD_geneset <- merge(NMD_geneset,NMD_Colombo_all, by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all.x = TRUE)
    NMD_geneset$NMD_Colombo <- ifelse(is.na(NMD_geneset$NMD_Colombo),"no","yes")  
    # Add transcripts purity correlations
    NMD_geneset <- merge(NMD_geneset,transcripts_purity_corr, by.x = "ensembl_transcript_id", by.y = "row.names", all.x = TRUE)
    # Add CNV_burden correlations
    NMD_geneset <- merge(NMD_geneset,CNV_burden_corr, by.x = "ensembl_transcript_id", by.y = "row.names", all.x = TRUE)
    # Remove duplicates
    NMD_geneset <- NMD_geneset[!duplicated(NMD_geneset),]
    # Create variable
    assign(NMD_geneset_name, NMD_geneset, envir = parent.frame())
}

# 1.2) Cell Line UPF1 KD
CL_UPF1KD_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/ENSEMBL/rsem.merged.transcript_tpm_ENSEMBL.tsv"
CL_metadata_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/CL_metadata_subset_and_controls.xlsx"
CL_UPF1KD <- read.table(file = CL_UPF1KD_path, header = TRUE, sep = "\t", row.names = 1)
CL_UPF1KD <- CL_UPF1KD[,-1]
CL_UPF1KD <- CL_UPF1KD[grep(".*PAR.*",rownames(CL_UPF1KD), invert = TRUE),]
rownames(CL_UPF1KD) <- sub("\\..*","", rownames(CL_UPF1KD))

# Merge replicates together
CL_UPF1KD_merged <- as.data.frame(matrix(nrow = nrow(CL_UPF1KD)), row.names = rownames(CL_UPF1KD) )
CL_UPF1KD_merged[,"GSE152435_C"] <- rowMeans(CL_UPF1KD[,c(1:5)])
CL_UPF1KD_merged[,"GSE152435"] <- rowMeans(CL_UPF1KD[,c(6:10)])
CL_UPF1KD_merged[,"GSE86148_C"] <- rowMeans(CL_UPF1KD[,c(11:13)])
CL_UPF1KD_merged[,"GSE86148"] <- rowMeans(CL_UPF1KD[,c(14:16)])
CL_UPF1KD_merged[,"GSE88140"] <- rowMeans(CL_UPF1KD[,c(17:18)])
CL_UPF1KD_merged[,"GSE88148"] <- rowMeans(CL_UPF1KD[,c(19:20)])
CL_UPF1KD_merged[,"GSE88266"] <- rowMeans(CL_UPF1KD[,c(21:22)])
CL_UPF1KD_merged[,"GSE88466"] <- rowMeans(CL_UPF1KD[,c(23:24)])
CL_UPF1KD_merged <- CL_UPF1KD_merged[,-1]

# Split controls from KD
CL_controls <- c("GSE152435_C","GSE86148_C","GSE88148","GSE88266")
CL_UPF1KD_controls <- CL_UPF1KD_merged[,colnames(CL_UPF1KD_merged)%in%CL_controls]
colnames(CL_UPF1KD_controls) <- c("GSE152435_WT1","GSE86148_WT2","GSE88148_WT3","GSE88266_WT4")
CL_UPF1KD_NMD <- CL_UPF1KD_merged[,!colnames(CL_UPF1KD_merged)%in%CL_controls]
colnames(CL_UPF1KD_NMD) <- c("GSE152435_KD1","GSE86148_KD2","GSE88140_KD4","GSE88466_KD3")
CL_UPF1KD_NMD <- CL_UPF1KD_NMD[,c(1,2,4,3)]

# Add pseudocounts for avoiding Inf ratios
CL_UPF1KD_NMD <- CL_UPF1KD_NMD+0.01
CL_UPF1KD_controls <- CL_UPF1KD_controls+0.01

# 2) Final NMD_target - NMD_control matching

for (NMD_geneset_name in c(NMD_genesets,"non_NMD_genes_with_NMD_features")) {
 
    print(paste0("........................ NMD geneset ---> ", NMD_geneset_name, "........................"))
    # Filter for NMD geneset
    # We create a new negative control
    if ( NMD_geneset_name == "non_NMD_genes_with_NMD_features" ) {
        NMD_geneset <- eval(parse(text=paste0("non_NMD_genes")))
    } else {
        NMD_geneset <- eval(parse(text=paste0(NMD_geneset_name)))
    }

    # Check for each gene which isoform is NMD-target and which one could be used as control
    NMD_geneset$final_consensus <- NA

    for (ensembl.gene in unique(NMD_geneset$ensembl_gene_id)) {

        NMD.nonNMD.pairs <- data.frame(NMD_target=NA,nonNMD_control=NA)
        ensembl.NMD.features.ensembl.gene <- NMD_geneset[NMD_geneset$ensembl_gene_id%in%ensembl.gene,]

        # Skip if the gene has only 1 isoform
        # For the negative control geneset, just take two random isoforms from the gene 
        if (length(unique(ensembl.NMD.features.ensembl.gene$ensembl_transcript_id)) == 1) {
            print("Gene with only 1 isoform, skipping...")
            next
        } else if (NMD_geneset_name == "non.NMD.genes") {
            isoforms.index <- which( (ensembl.NMD.features.ensembl.gene$TCGA_pancancer_gene_exp >= 3) )
            # If no controls
            if (length(isoforms.index) <= 1 ) {
                print("Gene with less than 2 isoforms expressed (<= 3), skipping")
                next
            } else {
                ensembl.NMD.features.ensembl.gene.filt <- ensembl.NMD.features.ensembl.gene[isoforms.index,]
                isoforms <- sample(ensembl.NMD.features.ensembl.gene.filt$ensembl_transcript_id)[1:2]
                NMD.target <- isoforms[1]
                nonNMD.control <- isoforms[2]
            }
        } else {
            # Check first if the gene is among the Karousis dataset, then obtain its corresponding NMD-nonNMD pair (if possible)
            if (ensembl.gene %in% unique(Karousis_original_dataset$gene_id)) {
                Karousis.ensembl.transcript <- Karousis_original_dataset[Karousis_original_dataset$gene_id%in%ensembl.gene,]
                # Check NMD isoform
                NMD.isoforms <- unique(Karousis.ensembl.transcript[Karousis.ensembl.transcript$NMD_isoforms%in%ensembl.NMD.features.ensembl.gene$ensembl_transcript_id,"NMD_isoforms"])
                # Check non_NMD_isoform
                nonNMD.isoforms <- unique(Karousis.ensembl.transcript[Karousis.ensembl.transcript$non_NMD_isoforms%in%ensembl.NMD.features.ensembl.gene$ensembl_transcript_id,"non_NMD_isoforms"])
                if (length(NMD.isoforms) >= 1) {
                    # Check NMD features in a specific set of transcripts and select the best NMD target
                    NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NMD.isoforms, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
                } else if (length(NMD.isoforms) == 0) {
                    # Check NMD features in all the transcripts of the given gene and select the best NMD target
                    NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
                }
                if (length(nonNMD.isoforms) >= 1) {
                    # Check NMD features in a specific set of transcripts and select the best nonNMD control one
                    nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = nonNMD.isoforms, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
                } else if (length(nonNMD.isoforms) == 0) {
                    # Check NMD features in all the transcripts of the given gene and select the best nonNMD control one
                    nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
                    #nonNMD.control <- "NA"
                }
                # Check NMD features in all the transcripts of the given gene and select the best NMD target when the Karousis target has no NMD features
                if (NMD.target == "NA") {
                    NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
                }
                # Karousis test
                karousis_test(ensembl.gene = ensembl.gene, ensembl.NMD.features.ensembl.gene = ensembl.NMD.features.ensembl.gene)
            } else {
                # Check NMD features in all the transcript of the given gene and select the best NMD target
                NMD.target <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "no", gene_exp_NMD_target = 1, gene_exp_control = 3)
                # Check NMD features in all the transcript of the given gene and select the best nonNMD control one
                nonNMD.control <- check_NMD_features(ensembl_NMD_features = ensembl.NMD.features.ensembl.gene, set = NULL, control = "yes", gene_exp_NMD_target = 1, gene_exp_control = 3)
            }
        }
        # Final consensus
        NMD.nonNMD.pairs$NMD_target <- NMD.target
        NMD.nonNMD.pairs$nonNMD_control <- nonNMD.control

        # Take gene exp and make ratio
        if (NMD.target == "NA" | nonNMD.control == "NA") {
            next
        } else {
            NMD.target.gene.exp <- as.numeric(ensembl.NMD.features.ensembl.gene[ensembl.NMD.features.ensembl.gene$ensembl_transcript_id%in%NMD.target,"TCGA_pancancer_gene_exp"])
            control.gene.exp <- as.numeric(ensembl.NMD.features.ensembl.gene[ensembl.NMD.features.ensembl.gene$ensembl_transcript_id%in%nonNMD.control,"TCGA_pancancer_gene_exp"])
            NMD.nonNMD.ratio <- NMD.target.gene.exp/control.gene.exp 
            if ( ((NMD.nonNMD.ratio <= 0.9) && (NMD_geneset_name != "non.NMD.genes")) || (NMD_geneset_name == "non.NMD.genes") ) {
                print(NMD.nonNMD.pairs)
                NMD_geneset[NMD_geneset$ensembl_transcript_id%in%as.character(NMD.target),"final_consensus"] <- "NMD_target"
                NMD_geneset[NMD_geneset$ensembl_transcript_id%in%as.character(nonNMD.control ),"final_consensus"] <- "control"
            } else {
                print("NMD target has more gene expression than control transcript")
                next
            }
        }

    } # END of Gene Loop

    ###### MANUAL INSPECTION ######

    NMDeff_old <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/old/NMD_global_2_shared_ensembl_filt.txt",
                            header = TRUE, sep = "\t")

    NMD_geneset_final <- NMD_geneset[!is.na(NMD_geneset$final_consensus),]

    NMDeff_old[!NMDeff_old$ensembl_gene_id %in% NMD_geneset_final$ensembl_gene_id,"ensembl_gene_id"]
    dim(NMD_geneset_final[!NMD_geneset_final$ensembl_gene_id %in% NMDeff_old$ensembl_gene_id,])

    NMD_geneset[NMD_geneset$ensembl_gene_id %in% "ENSG00000136271",c("ensembl_transcript_id","TCGA_pancancer_gene_exp")]
    NMDeff_old[NMDeff_old$ensembl_gene_id %in% "ENSG00000136271",c("final_consensus","ensembl_transcript_id","gene_exp")]
    NMDeff_old[NMDeff_old$ensembl_gene_id %in% "ENSG00000136271",]
    NMD_geneset[NMD_geneset$ensembl_gene_id %in% "ENSG00000136271",]


    new <- NMD_geneset_final[!NMD_geneset_final$ensembl_transcript_id %in% NMDeff_old$ensembl_transcript_id,"ensembl_transcript_id"]
    NMD_geneset_final[NMD_geneset_final$ensembl_transcript_id %in% new,c("ensembl_transcript_id","TCGA_pancancer_gene_exp")]
    NMDeff_old[NMDeff_old$ensembl_transcript_id %in% new,c("ensembl_transcript_id","gene_exp")]

    ###### MANUAL INSPECTION ######

    # 1) Check for each gene all NMD_targets-NMD_controls combinations
    # Obtain certain information for making the matching later

    if (NMD_geneset_name == "non_NMD_genes") {
        NMD_targets_controls_match <- negative_control_match(NMD_geneset)
    } else {
        if (NMD_geneset_name == "non_NMD_genes_with_NMD_features") {
            NMD_geneset <- NMD_geneset[!is.na(NMD_geneset$NMD_type),]
        }
        NMD_combinations_list <- NMD_target_control_combinations(NMD_geneset)
        # Save the NMD_target vs NMD_control combinations
        output_list_file = paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_target_controls_combinations/",NMD_geneset_name,"_NMD_target_controls_combinations.RData")
        saveRDS(NMD_combinations_list, file = output_list_file)
    }
    # 2) Matching NMD_target-NMD_control

    # 2.2) Final matching
    for (dataset in c("TCGA","GTEx")) {
        if (NMD_geneset_name == "non_NMD_genes") {
            NMD_geneset$NMD_type_final <- NA
            NMD_targets_index <- NMD_geneset$ensembl_transcript_id %in% NMD_targets_controls_match$NMD_target
            NMD_geneset[NMD_targets_index,"NMD_type_final"] <- "NMD_target"
            NMD_controls_index <- NMD_geneset$ensembl_transcript_id %in% NMD_targets_controls_match$NMD_control
            NMD_geneset[NMD_controls_index,"NMD_type_final"] <- "NMD_control"
            NMD_geneset_final <- NMD_geneset[!is.na(NMD_geneset$NMD_type_final),]
            # Select randomly 250 genes
            random_genes <- sample(unique(NMD_geneset_final$ensembl_gene_id))[1:250]
            NMD_geneset_final <- NMD_geneset_final[NMD_geneset_final$ensembl_gene_id %in% random_genes, ]
        } else {
            NMD_geneset_final <- NMD_target_control_final_matching(NMD_combinations_list, dataset = dataset)
        }
        # Save
        NMD_geneset_final <- NMD_geneset_final[order(NMD_geneset_final$ensembl_gene_id,NMD_geneset_final$NMD_type_final),] 
        print(length(unique(NMD_geneset_final$ensembl_gene_id)))
        output_path <- paste0(NMD_targets_path,"/",NMD_geneset_name,"_ensembl_final_",dataset,".txt")
        write.table(NMD_geneset_final, file = output_path, quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
        # Plot gene expression of NMD targets and controls
        df_long <- NMD_geneset_final[,c("NMD_type_final","TCGA_pancancer_gene_exp","GTEx_pantissue_gene_exp")] %>% 
                    pivot_longer(cols = c(TCGA_pancancer_gene_exp, GTEx_pantissue_gene_exp), names_to = "dataset", values_to = "gene_exp")
        plot_path <- paste0(NMD_targets_path,"/",NMD_geneset_name,"_final_",dataset,"_gene_exp_boxplot.png")
        png(file = plot_path, width = 4000, height = 3000, res = 300)
        p <- ggplot(df_long, aes(x = NMD_type_final, y = log10(gene_exp), color = NMD_type_final)) +
            geom_violin(alpha = 0.5,draw_quantiles = c(0.25, 0.5, 0.75)) +
            stat_summary(fun = median_quartile, geom = 'point') +
            geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
            facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
            theme_classic() +
            geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
            ylab(paste0("log10(TPM)")) + ggtitle(paste0("Pantissue gene expression of --> ",NMD_geneset_name, " - fitted in ",dataset)) +
            theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=13, face="bold"),
            axis.title.y = element_text(color="black", size=13, face="bold"))
        print(p)
        dev.off()
    }
} # END for each NMD_geneset

#####################

TCGA_NMDeff_main <- read.table(file ="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_TCGA.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GTEx_NMDeff_main <- read.table(file ="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_GTEx.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMDeff_old <- read.table(file ="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/old/NMD_global_2_shared_ensembl_filt.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
outliers <- read.table(file ="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/old/NMD_global_2_shared_outliers.txt",
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
# Remove outliers 
NMDeff_old <- NMDeff_old[-which(NMDeff_old$ensembl_gene_id %in% outliers$V1)]

# Check outiers
table(outliers$V1 %in% TCGA_NMDeff_main$ensembl_gene_id)
table(outliers$V1 %in% GTEx_NMDeff_main$ensembl_gene_id)

# Check overlap
table(unique(NMDeff_old$ensembl_gene_id) %in% unique(TCGA_NMDeff_main$ensembl_gene_id))
table(unique(NMDeff_old$ensembl_gene_id) %in% unique(GTEx_NMDeff_main$ensembl_gene_id))

df <- NMDeff_old[-which(NMDeff_old$ensembl_gene_id %in% GTEx_NMDeff_main$ensembl_gene_id),]
