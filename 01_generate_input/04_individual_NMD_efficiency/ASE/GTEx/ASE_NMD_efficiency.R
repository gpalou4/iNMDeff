rm(list=ls())

ASE_germline_sample_file <- function(GTEx.tissue, GTEx.sample) {
  ASE.germline.path.sample <- gsub("\\[X\\]",GTEx.sample, paste0(GTEx.ASE.files.path,paths[paths$folder_or_object=="ASE_files","path_or_filename"]))
  if (file.exists(ASE.germline.path.sample)) {
    ASE.germline.sample <- read.table(file = ASE.germline.path.sample, header = TRUE, sep = "\t")
  } else {
    return(NULL)
  }
  return(ASE.germline.sample)
}

ASE_germline_sample_filtering <- function(ASE.germline.sample, mutation, VAF, filter = NULL, coverage = NULL) {
    # Filter variants by: nonsense (PTCs) and VAF or synonymous (negative control)
    if (mutation == "synonymous") {
      variant.type <- "synonymous_variant"
      ASE.germline.sample.filt <- ASE.germline.sample[ASE.germline.sample$VARIANT_ANNOTATION%in%variant.type,]
      # Remove weird ASE ratios (peaks near 0 and 1)
      #ASE.ratios <- ASE.germline.sample.filt$altCount/(ASE.germline.sample.filt$altCount+ASE.germline.sample.filt$refCount)
      #ASE.germline.sample.filt <- ASE.germline.sample.filt[-which(ASE.ratios <= 0.25 | ASE.ratios >= 0.75),]      
    } else if (mutation == "stopgain") {
      variant.type <- "stop_gained"
      ASE.germline.sample.filt <- ASE.germline.sample[ASE.germline.sample$VARIANT_ANNOTATION%in%variant.type,]
    }
    # PASS variants
    if (!is.null(filter)) {
      ASE.germline.sample.filt <- ASE.germline.sample.filt[ASE.germline.sample.filt$LOW_MAPABILITY == 0 & ASE.germline.sample.filt$GENOTYPE_WARNING == 0,]
    }
    # Minimum ASE coverage
    if (!is.null(coverage)) {
      ASE.germline.sample.filt <- ASE.germline.sample.filt[which(ASE.germline.sample.filt$TOTAL_COUNT >= coverage),]
    }
    #Subset in ASE synonymous
    # if (mutation == "synonymous") { 
    #   ASE.germline.sample.filt <- ASE.germline.sample.filt[sample(1:nrow(ASE.germline.sample.filt))[1:50],]
    # }
    # Obtain AF (NOT DONE YET)    
    return(ASE.germline.sample.filt)
}

VCF_germline_sample_file <- function(GTEx.tissue, GTEx.sample) {
  VCF.germline.path.sample <- gsub("\\[X\\]",GTEx.sample, paste0(GTEx.VCF.files.path,paths[paths$folder_or_object=="VCF_files","path_or_filename"]))
  if (file.exists(VCF.germline.path.sample)) {
    VCF.tmp <- paste0(GTEx.VCF.files.path,GTEx.sample,"_",GTEx.tissue,"_VCF_tmp.txt")
    system(paste0("zcat ",VCF.germline.path.sample," | head -1 > ",VCF.tmp), intern = TRUE)
    system(paste0("zcat ",VCF.germline.path.sample," | grep '[^non]synonymous SNV' >> ",VCF.tmp), intern = TRUE)
    VCF.germline.sample <- read.table(file = VCF.tmp, header = TRUE, sep = "\t", comment.char = "")
    VCF.germline.sample$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",VCF.germline.sample$Gene.refGene)
    system(paste0("rm ",VCF.tmp))
  } else {
    return(NULL)
  }
  return(VCF.germline.sample)
}

glm_nb_regression <- function(NMD.final.df) {
  
  # Create Dataframe
  df1 <- data.frame(gene.exp = NMD.final.df$refCount, gene.id=NMD.final.df$gene_id, NMD.target = 0)
  df2 <- data.frame(gene.exp = NMD.final.df$altCount, gene.id=NMD.final.df$gene_id, NMD.target = 1)
  df <- rbind(df1,df2)

  # Bayesian NB regression model
  if (length(unique(df$gene.id))==1) {
    glm.nb.model <- "stan_glm(gene.exp ~ as.factor(NMD.target), data = df, family = neg_binomial_2, verbose = FALSE)"
  } else {
    glm.nb.model <- "stan_glm(gene.exp ~ as.factor(NMD.target) + as.factor(gene.id), data = df, family = neg_binomial_2, verbose = FALSE)"
  }

  # Perform the regression
  tryCatch( { 
    nb.res <- NA
    nb.res <- eval(parse(text = glm.nb.model))
  }, error = function(e) {
    print("################################################################################################")
    print("ERROR in NB regression. Model:")
    print(e)
    nb.error.3 <<- TRUE
    print(glm.nb.model)
    print(head(df))
    print("################################################################################################")
  })
  return(nb.res)
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

# Library home

#install.packages("",lib="/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6")
#.libPaths("/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6")
.libPaths( rev( .libPaths() ) )
# Bayesian for NB regression
library("V8")
# library("shinystan")
library("rstan")
library("rstanarm")
#.libPaths( c( "/g/strcombio/fsupek_home/gpalou/R/x86_64-redhat-linux-gnu-library/3.6" , .libPaths() ) )
# Library Size and VST
library("DBI")
library("DESeq2")
# PCA
library("factoextra")
library("FactoMineR")
#library("RColorBrewer")
# Split
library("stringr")
# Read GTF
library("rtracklayer")
# To read fasta
library("seqinr")
# ORFs
library("ORFik")
# Random strings
library("stringi")

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

#paths <- read.table(file = "/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/ASE/GTEx/ASE_NMD_efficiency_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#GTEx.tissue <- "Brain_Spinal_cord_cervical_c1"

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.tissue <- args[2]
nextflow.params <- args[3:length(args)]
#nextflow.params <- c("0.10","0")
VAF <- as.numeric(nextflow.params[1])
LOEUF.score <- nextflow.params[2]
if (LOEUF.score != "no") {
  LOEUF.score <- as.numeric(LOEUF.score)
}

GTEx.names.path <- paths[paths$folder_or_object=="GTEx_names_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
RNAseq.path <- paths[paths$folder_or_object=="RNAseq","path_or_filename"]
GTEx.VCF.files.path <- paths[paths$folder_or_object=="VCF_files_path","path_or_filename"]
GTEx.ASE.nb.res.path <- paths[paths$folder_or_object=="ASE_nb_res_path","path_or_filename"]
#transcripts.germ.mut.path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
germline.PTCs.transcripts.path <- paths[paths$folder_or_object=="germline_PTCs_transcripts_path","path_or_filename"]
nb.coeff.res.path <- paths[paths$folder_or_object=="ASE_nb_res_path","path_or_filename"]
GTEx.ASE.files.path <- paths[paths$folder_or_object=="ASE_files_path","path_or_filename"]

print(paste0("Tissue --> ",GTEx.tissue))

# 1.1) Transcripts length
ensembl.v88.transcripts.length <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_transcripts_length","path_or_filename"]), 
                                header = TRUE, sep = "\t")
# 1.2) Protein coding transcripts
ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t")    
# 1.3) Positive selected genes (PSG)
PSG.ensembl.gene <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="PSG_331","path_or_filename"]),
                                   header = TRUE, sep = "\t")
# 1.4) NMD targets (for negative control)
NMD.global.ensembl <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_ensembl","path_or_filename"]),
                                   header = TRUE, sep = "\t")
# 1.5) GTEx tissues
GTEx.tissues <- read.table(file = paste0(GTEx.names.path,paths[paths$folder_or_object=="GTEx_names","path_or_filename"]))$V1
# Acronyms of tissues
acronyms <- c("ADPSBQ","ADPVSC","ADRNLG","ARTAORT","ARTCRN","ARTTBL","BLDDER","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB",
              "BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG","BREAST","FIBRBLS","LCL","CML","CVXECT","CVSEND","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL",
              "FLLPNT","HRTAA","HRTLV","KDNCTX","KDNMDL","LIVER","LUNG","SLVRYG","MSCLSK","NERVET","OVARY","PNCREAS","PTTARY","PRSTTE","SKINNS","SKINS",
              "SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")
tissues.acronyms <- data.frame(tissues = GTEx.tissues, acronyms = acronyms )
GTEx.tissue.acronym <- tissues.acronyms[tissues.acronyms$tissues%in%GTEx.tissue,"acronyms"]
# 1.6) LOEUF score (negatively selected genes)
# This is to remove constrained transcripts, i.e. those with potentially negative selection. Why? Because if not, dosage compensation from the other wild-type parental allele might happen and can mask the NMD effect
# If the transcript is not negatively selected, no need for dosage compensation!
LOEUF.score.table <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="LOEUF_score","path_or_filename"]), header = TRUE, sep = "\t")
# oe_lof_upper --> LOEUF: upper bound of 90% confidence interval for o/e ratio for pLoF variants (lower values indicate more constrained)
# oe_lof_upper_bin --> Decile bin of LOEUF for given transcript (lower values indicate more constrained)
LOEUF.score.ensembl <- LOEUF.score.table[,c("gene_id","oe_lof_upper","oe_lof_upper_bin")]
LOEUF.genes.to.remove <- unique(na.omit(LOEUF.score.ensembl[which(LOEUF.score.ensembl$oe_lof_upper_bin<=LOEUF.score),"gene_id"]))

# 3) Negative Binomial regression on ASE data
columns <- c("sample","subtype","num.PTCs","stopgain","stopgain.NMD.evading",
            "stopgain.pval","synonymous","synonymous.pval","error")
nb.coeff.res <- data.frame(matrix(ncol=length(columns),nrow=0, dimnames=list(NULL, columns)))
PTC.set.names <- c("stopgain","synonymous")
total.num.PTC.sample <- list()

# Obtain samples names with available ASE germline data
GTEx.ASE.files <- list.files(GTEx.ASE.files.path, recursive = FALSE, full.names = FALSE)
GTEx.ASE.files <- GTEx.ASE.files[grep("wasp",GTEx.ASE.files)]
GTEx.ASE.sample.names <- gsub("(GTEX-\\w{4,5})\\..*gz$","\\1",GTEx.ASE.files)

for (sample in seq(1,length(GTEx.ASE.sample.names))) {

  start_time_total <- Sys.time()
  GTEx.sample <- GTEx.ASE.sample.names[sample]
  nb.coeff.res[sample,"sample"] <- GTEx.sample
  print(paste("SAMPLE", sample, " ----> ",GTEx.sample),sep = "")
  
  # 1) Obtain Germline PTCs (from nonsense + FS) with predicted NMD-rules for the sample
  germline.PTCs.transcripts.sample.path <- paste0(germline.PTCs.transcripts.path,paths[paths$folder_or_object=="germline_PTCs_transcripts","path_or_filename"])
  germline.PTCs.transcripts.sample.path <- gsub("\\[X1\\]",GTEx.tissue, germline.PTCs.transcripts.sample.path)
  germline.PTCs.transcripts.sample.path <- gsub("\\[X2\\]",GTEx.sample, germline.PTCs.transcripts.sample.path)

  if (file.exists(germline.PTCs.transcripts.sample.path)) {
    germline.PTCs.transcripts.sample <- read.table(file = germline.PTCs.transcripts.sample.path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    print("germline PTCs file do not exist, skipping...")
    next
  }

  # 2) Obtain ASE germline data file
  # Check if file exists
  ASE.germline.sample <- NULL
  ASE.germline.sample <- ASE_germline_sample_file(GTEx.sample = GTEx.sample, GTEx.tissue = GTEx.tissue)
  if(is.null(ASE.germline.sample)){
    nb.coeff.res[sample,"error"] <- "ASE germline file not found"
    print("ASE germline file not found, skipping...")
    next
  }
  # If the variant only appears in one tissue that is a somatic variant not germline
  somatic_variants <- names(which(table(ASE.germline.sample$VARIANT_ID) == 1))
  # Filter variants for the specific tissue
  ASE.germline.sample.tissue <- ASE.germline.sample[ASE.germline.sample$TISSUE_ID%in%GTEx.tissue.acronym,]
  if(nrow(ASE.germline.sample)==0){
    print("Sample with no mutations in this tissue")
    next
  }
  ASE.germline.sample.tissue <- ASE.germline.sample.tissue[-which(ASE.germline.sample.tissue$VARIANT_ID %in% somatic_variants),]

  # 3) NB regression per set
  for (PTC.set.name in PTC.set.names) {
    for (type in c("NMD-triggering","NMD-evading")) {
      if (PTC.set.name == "synonymous" & type == "NMD-evading") {next}
      # Init
      start_time <- Sys.time()
      nb.res <- NA
      print(paste0("------------------------",toupper(PTC.set.name),"------------------------------"))
      print(paste0("------------------------",toupper(type),"------------------------------"))
      # 4) Filter VCF for variants
      if (PTC.set.name == "stopgain") {
        ASE.germline.sample.filt <- ASE_germline_sample_filtering(ASE.germline.sample = ASE.germline.sample.tissue,
                                    mutation = PTC.set.name, VAF = 1, filter = "PASS", coverage = 0)
        if (nrow(ASE.germline.sample.filt) == 0 ){
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" )  {
            nb.coeff.res[sample,"error"] <- "<1 total genes"
          }
          next
        }
        # Change ASE columns
        ASE.germline.sample.save <- ASE.germline.sample.filt[,c("CHR","POS","REF_ALLELE","ALT_ALLELE","REF_COUNT","ALT_COUNT","TOTAL_COUNT")]
        colnames(ASE.germline.sample.save) <- c("CHROM","START","REF_DNA","ALT_DNA","refCount_ASE","altCount_ASE","totalCount_ASE")
        # Merge with PTCs metadata
        germline.PTCs.transcripts.sample.ASE <- merge(germline.PTCs.transcripts.sample,ASE.germline.sample.save, by.x = c("chr","start_pos","Ref","Alt"), 
                by.y = c("CHROM","START","REF_DNA","ALT_DNA"), all.x = TRUE)

        # Save raw PTCs_ASE df for sample
        germline.PTCs.ASE.transcripts.tissue.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(germline.PTCs.transcripts.path,paths[paths$folder_or_object=="germline_PTCs_ASE_transcripts","path_or_filename"]))
        germline.PTCs.ASE.transcripts.sample.path <- gsub("\\[X2\\]",GTEx.sample, germline.PTCs.ASE.transcripts.tissue.path)
        #if ( (!file.exists(germline.PTCs.ASE.transcripts.sample.path)) & ( PTC.set.name == "stopgain" & type == "NMD-triggering"  & VAF == 0.01) ) {
        if ( ( PTC.set.name == "stopgain" & type == "NMD-triggering" & VAF == 0.2) ) {
          write.table(germline.PTCs.transcripts.sample.ASE, file = germline.PTCs.ASE.transcripts.sample.path, 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
        }
        # 5) Filterings for the NB regression
        # 5.1) Coverage 5 for nonsense and 2 for FS (filter for ASE PTCs)
        filter.coverage.nonsense <- which(germline.PTCs.transcripts.sample.ASE$totalCount_ASE >= 5 & germline.PTCs.transcripts.sample.ASE$stopgain == "nonsense")
        filter.coverage.FS <- which(germline.PTCs.transcripts.sample.ASE$totalCount_ASE >= 2 & germline.PTCs.transcripts.sample.ASE$stopgain != "nonsense")
        germline.PTCs.transcripts.sample.ASE.filt <- germline.PTCs.transcripts.sample.ASE[c(filter.coverage.nonsense,filter.coverage.FS),]
        # 5.2) Remove ASE-PTCs from FS with no predicted PTC
        if (length(filter.coverage.FS) != 0) {
          filter.FS.no.PTC <- which(germline.PTCs.transcripts.sample.ASE.filt$stopgain != "nonsense" & !is.na(germline.PTCs.transcripts.sample.ASE.filt$PTC_CDS_pos)) 
        } else {
          filter.FS.no.PTC <- NULL
        }
        # 5.3) Remove ASE-PTCs with specified VAF
        filter.VAF <- which(germline.PTCs.transcripts.sample.ASE.filt$VAF > VAF)
        # 5.4) Remove transcripts with 1 exon
        filter.single.exon <- which(germline.PTCs.transcripts.sample.ASE.filt$transcript_CDS_exon_num == 1 & germline.PTCs.transcripts.sample.ASE.filt$splice_site_3UTR == "no")
        # 5.5) LOEUF score 1
        filter.LOEUF <- which(germline.PTCs.transcripts.sample.ASE.filt$LOEUF_decile <= LOEUF.score)
        # 5.6) PSG genes
        filter.PSG <- which(germline.PTCs.transcripts.sample.ASE.filt$PSG == "yes")
        # 5.6) Coefficient variation
        # Coefficient of variation (mean/SD) > 0.5 --> This reduces a lot the NMDeff values (from 1.6 raw to 0.3)
        #filter.coeff.var <- which(germline.PTCs.transcripts.sample.ASE.filt$coeff_var > 0.5)
        # 5.7) NMD-triggering vs NMD-evading
        if (type == "NMD-triggering") {
          filter.NMDtype <- which(
                germline.PTCs.transcripts.sample.ASE.filt$X55_nt_last_exon == "NMD-evading" |
                germline.PTCs.transcripts.sample.ASE.filt$TSS_PTC_dist <= 200
                  )
        } else if (type == "NMD-evading") {
          # filter.NMDtype <- which(germline.PTCs.transcripts.sample.ASE.filt$X55_nt_last_exon == "NMD-triggering")
          filter.NMDtype <- which(
                germline.PTCs.transcripts.sample.ASE.filt$X55_nt_last_exon == "NMD-triggering" &
                germline.PTCs.transcripts.sample.ASE.filt$TSS_PTC_dist > 200
                  )
        }
        # 5.9) Non protein-coding
        filter.nonCoding <- which(germline.PTCs.transcripts.sample.ASE.filt$protein_coding == "no")
        # 5.10) Overlapping germline mutations...
        # 5.11) Take only Het variants
        filter.homVariants <- which(!germline.PTCs.transcripts.sample.ASE.filt$genotype %in% c("0/1","1|0","0|1"))
        # Add up filters
        all.filters <- unique(c(filter.FS.no.PTC,filter.VAF,filter.single.exon,filter.LOEUF,filter.NMDtype,filter.nonCoding,filter.homVariants))
        germline.PTCs.transcripts.sample.ASE.filt <- germline.PTCs.transcripts.sample.ASE.filt[-all.filters,]
        if (nrow(germline.PTCs.transcripts.sample.ASE.filt) == 0 ){
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" )  {
            nb.coeff.res[sample,"error"] <- "<1 total genes"
          }
          next
        } else {
          NMD.final.df <- germline.PTCs.transcripts.sample.ASE.filt[,c("start_pos","gene_id","refCount_ASE","altCount_ASE","VAF")]
          colnames(NMD.final.df) <- c("start_pos","gene_id","refCount","altCount","AF")
          # Remove redundant transcripts with same ASE counts
          NMD.final.df <- NMD.final.df[which(!duplicated(NMD.final.df)),]
        }
      } else if (PTC.set.name == "synonymous") {
        ASE.germline.sample.filt <- ASE_germline_sample_filtering(ASE.germline.sample = ASE.germline.sample.tissue, 
                                    mutation = PTC.set.name, VAF = VAF, filter = "PASS", coverage = 5)
        if(nrow(ASE.germline.sample.filt)==0){
          print("No variants left, skipping...")
          next
        }
        # VCF (this is only to obtain the VAF...)
        VCF.germline.sample <- VCF_germline_sample_file(GTEx.tissue, GTEx.sample)
        # Change ASE columns
        ASE.germline.sample.filt <- ASE.germline.sample.filt[,c("CHR","POS","REF_ALLELE","ALT_ALLELE","REF_COUNT","ALT_COUNT","TOTAL_COUNT")]
        colnames(ASE.germline.sample.filt) <- c("CHROM","START","REF_DNA","ALT_DNA","refCount_ASE","altCount_ASE","totalCount_ASE")
        # Merge with PTCs metadata
        germline.PTCs.transcripts.sample.ASE <- merge(VCF.germline.sample,ASE.germline.sample.filt, by.x = c("Chr","Start","Ref","Alt"), 
                by.y = c("CHROM","START","REF_DNA","ALT_DNA"), all.x = TRUE)
        filter <- colnames(germline.PTCs.transcripts.sample.ASE) %in% c("Start","Gene.refGene","refCount_ASE","altCount_ASE")
        colnames(germline.PTCs.transcripts.sample.ASE)[filter] <- c("start_pos","gene_id","refCount","altCount")
        # Remove NA's
        germline.PTCs.transcripts.sample.ASE <- germline.PTCs.transcripts.sample.ASE[which(!is.na(germline.PTCs.transcripts.sample.ASE$refCount)),]
        # 5) Filterings for the NB regression 
        # 5.1) LOEUF score
        filter.LOEUF <- which(germline.PTCs.transcripts.sample.ASE$gene_id %in% LOEUF.genes.to.remove)
        # 5.2) Non protein-coding
        genes.id <- strsplit(germline.PTCs.transcripts.sample.ASE$gene_id,";")
        genes.id.df <- as.data.frame(do.call(rbind,genes.id))
        filter.nonCoding <- which( !as.character(genes.id.df$V1)%in%ensembl.v88.coding.transcripts$V2 | !as.character(genes.id.df$V2)%in%ensembl.v88.coding.transcripts$V2 )
        # 5.2) PSG genes
        filter.PSG <- which( as.character(genes.id.df$V1)%in%PSG.ensembl.gene$Ensembl.Gene.ID | as.character(genes.id.df$V2)%in%PSG.ensembl.gene$Ensembl.Gene.ID )
        # 5.3) NMD genes
        NMD.genes <- unique(germline.PTCs.transcripts.sample.ASE[as.character(germline.PTCs.transcripts.sample.ASE$gene_id)%in%NMD.global.ensembl$ensembl_gene_id,"gene_id"])
        filter.NMD <- which(as.character(germline.PTCs.transcripts.sample.ASE$gene_id)%in%NMD.genes)
        # 5.4) Remove ASE-PTCs with specified VAF
        filter.VAF <- which(germline.PTCs.transcripts.sample.ASE$AF > VAF)
        # Add up filters
        all.filters <- unique(c(filter.LOEUF,filter.NMD,filter.nonCoding,filter.VAF,filter.PSG))
        germline.PTCs.transcripts.sample.ASE.filt <- germline.PTCs.transcripts.sample.ASE[-all.filters,]
        # 5.5) Remove redundant transcripts with same ASE counts
        germline.PTCs.transcripts.sample.ASE.filt <- germline.PTCs.transcripts.sample.ASE.filt[,c("start_pos","gene_id","refCount","altCount","AF")]
        NMD.final.df <- germline.PTCs.transcripts.sample.ASE.filt[which(!duplicated(germline.PTCs.transcripts.sample.ASE.filt)),]
        if (nrow(NMD.final.df) == 0 ){
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" )  {
            nb.coeff.res[sample,"error"] <- "<1 total genes"
          }
          next
        }
        # 5.6) Sample 100 variants maximum
        if ( nrow(NMD.final.df) > 100) { 
          NMD.final.df <- NMD.final.df[sample(1:nrow(NMD.final.df))[1:100],]
        }
      }
      # Number of PTC genes for the regression
      if ( PTC.set.name == "stopgain" & type == "NMD-triggering" ) {
        num.PTC <- nrow(NMD.final.df)
        nb.coeff.res[sample,"num.PTCs"] <- num.PTC
        print(paste0("Final number of PTCs --> ", num.PTC))
        total.num.PTC.sample[[GTEx.sample]] <- num.PTC
      }
      # NB regression
      nb.res <- glm_nb_regression(NMD.final.df = NMD.final.df)
      # Save results
      end_time <- Sys.time()
      print(paste0(PTC.set.name," time taken --> ",end_time - start_time))
      if (is.na(nb.res)) {
        nb.res.pval <- NA
        nb.res.coeff <- NA
        nb.coeff.res[sample,paste0(PTC.set.name,".sd")] <- NA
        nb.coeff.res[sample,paste0(PTC.set.name,".CI_2.5")] <- NA
        nb.coeff.res[sample,paste0(PTC.set.name,".CI_97.5")] <- NA
        if (PTC.set.name == "stopgain") {
          nb.coeff.res[sample,"error"] <- "NB regression error" 
        }
      } else if (!is.na(nb.res)) {
          nb.res.coeff <- as.numeric(coef(nb.res)[2])
          print(paste0("NMD efficiency --> ", nb.res.coeff))
          # Store NB coefficient
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" ) {
            nb.coeff.res[sample,PTC.set.name] <- nb.res.coeff
            col <- PTC.set.name
          } else if ( PTC.set.name == "stopgain" & type == "NMD-evading" ) {
            nb.coeff.res[sample,paste0(PTC.set.name,".",gsub("-","\\.",type))] <- nb.res.coeff
            col <- paste0(PTC.set.name,".",gsub("-","\\.",type))
          } else if ( PTC.set.name == "synonymous" & type == "NMD-triggering" ) { 
            nb.coeff.res[sample,PTC.set.name] <- nb.res.coeff
            col <- PTC.set.name
          } 
          nb.coeff.res[sample,paste0(col,".sd")] <- as.numeric(nb.res$stan_summary[2,"sd"])
          nb.coeff.res[sample,paste0(col,".CI_2.5")] <- as.numeric(nb.res$stan_summary[2,"2.5%"])
          nb.coeff.res[sample,paste0(col,".CI_97.5")] <- as.numeric(nb.res$stan_summary[2,"97.5%"])
      }
    }
  }
  end_time_total <- Sys.time()
  print(paste0("TOTAL time take for the sample -->",end_time_total - start_time_total))
}

file.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(nb.coeff.res.path,paths[paths$folder_or_object=="NB_coeff_res","path_or_filename"]))
file.path <- gsub("\\[X2\\]",as.character(VAF), file.path)
# write.table(nb.coeff.res, file = file.path , sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = FALSE)

# Obtain length for each PTC set

PTC.sets.lengths <- "c("
for (PTC.set.name in PTC.set.names) {
  PTC.set.length <- nrow(nb.coeff.res[!is.na(nb.coeff.res[,PTC.set.name]),])
  PTC.sets.lengths <- paste0(PTC.sets.lengths,",\'",PTC.set.name," (",PTC.set.length,")\'")
}
PTC.sets.lengths <- sub(",","",PTC.sets.lengths)
PTC.sets.lengths <- paste0(PTC.sets.lengths,")")
PTC.sets.lengths <- eval(parse(text = PTC.sets.lengths ))

plot.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(nb.coeff.res.path,paths[paths$folder_or_object=="NB_coeff_boxplot","path_or_filename"]))
plot.path <- gsub("\\[X2\\]",as.character(VAF), plot.path)

png(file = plot.path, width = 4000, height = 3000, res = 300)
nb.coeff.res.filt <- stack(nb.coeff.res[,PTC.set.names])
ggplot(data = nb.coeff.res.filt, aes(x=ind, y = exp(values), color = ind)) +
  geom_violin(width=0.75, trim=FALSE) + #coord_flip(ylim = c(0,2)) +  
  geom_boxplot(width=0.1, fill="white")+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5) +
  ylab("Targets/controls ratio") + xlab("PTC sets") + ggtitle(paste0("NMD efficiency by PTC set --> ",GTEx.tissue)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        legend.position = "none", axis.text=element_text(size=15)) +
  scale_x_discrete(labels=PTC.sets.lengths) + ylim(c(0,2))
dev.off()

# # Gene expression boxplot PTCs NMD-triggering vs NMD-evading
# plot.path <- gsub("\\[X\\]",gsub("TCGA-","",GTEx.tissue), paste0(nb.coeff.res.path,paths[paths$folder_or_object=="gene_exp_boxplot","path_or_filename"]))
# plot.path <- gsub("\\[X2\\]",paste0(as.character(VAF),"-",method),plot.path)
# # file.path <- gsub("\\[X3\\]",LOEUF.score,file.path)
# png(file = plot.path, width = 4000, height = 3000, res = 300)
# boxplot(log2(nb.coeff.res$PTCs.gene.exp.NMD.triggering),log2(nb.coeff.res$PTCs.gene.exp.NMD.evading), names = c("PTCs triggering","PTCs evading"), ylab = "Gene exp TPM medians", ylim = c(-5,5), main = GTEx.tissue)
# dev.off()












