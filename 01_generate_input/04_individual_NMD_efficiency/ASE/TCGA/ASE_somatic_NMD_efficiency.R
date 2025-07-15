rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

ASE_extract_info <- function(ASE) {
    variants.list <- list()
    for (i in 1:nrow(ASE)) {
      # Transcripts list from same variant
      transcripts.list <- strsplit(as.character(ASE[i,"AAChange.refGene"]), split=",")
      # Obtain the ID name and other info (exon location, PTC CDS location)
      transcripts.exons.info <- lapply(transcripts.list, function(transcript) {
        variant.CDS.exon.num <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon([0-9]{1,3}).*","\\1",transcript)
        transcript.id <- gsub(".*(ENST[0-9]{11})\\.[0-9]{1,3}:exon.*","\\1",transcript)
        #variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)\\D*:p\\..*","\\1",transcript)
        variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)_?.*:p\\..*","\\1",transcript)
        data.frame(transcript_id=transcript.id,exon_num=variant.CDS.exon.num, variant_CDS_pos = variant.CDS.pos)
      })
      transcripts.exons.info <- transcripts.exons.info[[1]]
      transcripts.exons.info$start_pos <- ASE[i,"POS"]
      transcripts.exons.info$gene_id <- ASE[i,"Gene.refGene"]
      transcripts.exons.info$refCount <- ASE[i,"refCount_RNA"]
      transcripts.exons.info$altCount <- ASE[i,"altCount_RNA"]
      transcripts.exons.info$AF <- ASE[i,"AF"]
      transcripts.exons.info$chr <- ASE[i,"CHROM"]
      # Type of stopgain (nonsense or frameshift)
      ref.allele <- as.character(ASE[i,"REF"])
      mut.allele <- as.character(ASE[i,"ALT"])
      transcripts.exons.info$Ref <- ref.allele
      transcripts.exons.info$Alt <- mut.allele
      if (mut.allele == "-") {
        mut.allele <- ""
      }
      if (ref.allele == "-"){
        ref.allele <- ""
      }
      ref.allele.num <- nchar(ref.allele)
      mut.allele.num <- nchar(mut.allele)
      if ( ref.allele.num < mut.allele.num ) {
        transcripts.exons.info$stopgain <- "frameshift_insertion"
      } else if ( ref.allele.num > mut.allele.num ) {
        transcripts.exons.info$stopgain <- "frameshift_deletion"
      } else {
        transcripts.exons.info$stopgain <- "nonsense"
      }
      variants.list[[i]] <- transcripts.exons.info
  }
  return(variants.list)
}

ASE_somatic_sample_file <- function(TCGA.cancer, TCGA.barcode, method) {
  if (method == "strelka") {
    ASE.somatic.path.sample <- gsub("\\[X2\\]",gsub("\\.","-", TCGA.barcode), paste0(ASE.somatic.path.cancer,paths[paths$folder_or_object=="ASE_somatic","path_or_filename"]))
  } else if (method == "samtools") {
    ASE.somatic.path.sample <- gsub("\\[X2\\]",gsub("\\.","-", TCGA.barcode), paste0(ASE.somatic.path.cancer,paths[paths$folder_or_object=="ASE_somatic","path_or_filename"]))
    ASE.somatic.path.sample <- gsub("ASE_somatic_strelka_calls.txt","ase_somatic_counts_VCF_added.txt",ASE.somatic.path.sample)
  }
  if (file.exists(ASE.somatic.path.sample)) {
    error <- FALSE
    tryCatch( { 
      ASE.somatic.sample <- read.table(file = ASE.somatic.path.sample, header = TRUE, sep = "\t")      
      #ASE.somatic.sample <- read.csv(file = ASE.somatic.path.sample, header = TRUE, sep = "\t")
    }, error = function(e) {
      print(e)
      error <<- TRUE
    })
  } else {
    return(NULL)
  }
  if (isTRUE(error)) {
    return(NULL)
  }
  # Fix VAF numeric
  VAFs <- ASE.somatic.sample$AF
  VAFs <- ifelse(VAFs == "" | VAFs == ".",0,format(VAFs,scientific = TRUE))
  VAFs <- as.numeric(VAFs)
  ASE.somatic.sample$AF <- VAFs
  # Fix ENSEMBL Genes IDs
  # Some transcript swill have different Genes IDs (for whatever reason), but bc it's only used on the regression as gene.id it's not problem
  ASE.somatic.sample$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",ASE.somatic.sample$Gene.refGene)
  return(ASE.somatic.sample)
}

ASE_somatic_sample_filtering <- function(ASE.somatic.sample, mutation, VAF, filter = NULL, coverage = NULL, method) {

    # Filter variants by: nonsense/FSs (PTCs) and VAF or synonymous (negative control)
    if (mutation == "synonymous") {
      variant.type <- "synonymous SNV"
      ASE.somatic.sample.filt <- ASE.somatic.sample[ASE.somatic.sample$ExonicFunc.refGene==variant.type & ASE.somatic.sample$AF <= VAF,]
    } else if (mutation == "stopgain") {
      variant.type <- c("stopgain","frameshift insertion","frameshift deletion")
      ASE.somatic.sample.filt <- ASE.somatic.sample[ASE.somatic.sample$ExonicFunc.refGene%in%variant.type &  ASE.somatic.sample$AF <= VAF,]
    }
    # Fix column names
    if (method == "samtools") {
      colnames(ASE.somatic.sample.filt)[colnames(ASE.somatic.sample.filt) %in% c("position","refCount","altCount","Otherinfo10")] <- c("POS","refCount_RNA","altCount_RNA","Otherinfo10_DNA")
      # Remove variants with ALT == 0 (we don't know if these are false positives or true NMD-ed)
      ASE.somatic.sample.filt <- ASE.somatic.sample.filt[ASE.somatic.sample.filt$altCount_RNA != 0,]
    }
    # PASS variants
    if (!is.null(filter)) {
      ASE.somatic.sample.filt <- ASE.somatic.sample.filt[ASE.somatic.sample.filt$Otherinfo10_DNA == filter,]
    }
    # Minimum ASE coverage
    if (!is.null(coverage)) {
      ASE.coverage <- (ASE.somatic.sample.filt$altCount_RNA+ASE.somatic.sample.filt$refCount_RNA)
      ASE.somatic.sample.filt <- ASE.somatic.sample.filt[which(ASE.coverage >= coverage),]
    }
    #Subset in ASE synonymous
    # if ( (mutation == "synonymous") && (nrow(ASE.somatic.sample.filt) > 100) ) { 
    #   ASE.somatic.sample.filt <- ASE.somatic.sample.filt[sample(1:nrow(ASE.somatic.sample.filt))[1:100],]
    # }
    return(ASE.somatic.sample.filt)
}

somatic_transcripts_filtering <- function(ASE.somatic, transcripts.filtering.matrix, TCGA.barcode) {
  if (TCGA.barcode%in%colnames(transcripts.filtering.matrix)) {
    sample.transcripts.to.remove <- rownames(transcripts.filtering.matrix[transcripts.filtering.matrix[,TCGA.barcode]==1,])
    genes.to.remove <- unique(ASE.somatic[as.character(ASE.somatic$transcript_id)%in%sample.transcripts.to.remove,"gene_id"])
    filter.somatic.mut.CNV <- which(ASE.somatic$gene_id%in%genes.to.remove)
    return(filter.somatic.mut.CNV)
  } else {
    return(NULL)
  }
}

germline_transcripts_filtering <- function(ASE.germline, transcripts.filtering.matrix, sample) {
  germ.mut <- "splicing|nonframeshift deletion|nonframeshift insertion|startloss|stoploss"
  transcripts.germ.mut.sample <- transcripts.filtering.matrix[,sample, drop = FALSE]
  transcripts.germ.mut.sample <- transcripts.germ.mut.sample[grep(germ.mut, transcripts.germ.mut.sample[,1]),, drop = FALSE]
  # Don't remove transcripts with PTC mutations
  PTC.mut <- "stopgain|[^non]frameshift deletion|[^non]frameshift insertion"
  transcripts.germ.mut.sample <- transcripts.germ.mut.sample[grep(PTC.mut, transcripts.germ.mut.sample[,1], invert = TRUE),, drop = FALSE]
  sample.transcripts.to.remove <- rownames(transcripts.germ.mut.sample)
  genes.to.remove <- unique(ASE.somatic[!as.character(ASE.somatic$transcript_id)%in%sample.transcripts.to.remove,"gene_id"])
  ASE.somatic <- ASE.somatic[!ASE.somatic$gene_id%in%genes.to.remove,]
  return(ASE.somatic)
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
    #print(paste0("purity: ",offset3))
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

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
TCGA.cancer <- args[2]
nextflow.params <- args[3:length(args)]

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/ASE/TCGA/ASE_NMD_efficiency_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#TCGA.cancer <- "TCGA-KIRC"
##nextflow.params <- c("somatic_transcripts_filt","somatic_transcripts_filt","0.01","strelka")
#nextflow.params <- c("somatic_transcripts_filt","0.2","strelka","0")
VAF <- as.numeric(nextflow.params[2])
method <- as.character(nextflow.params[3])
LOEUF.score <- nextflow.params[4]
if (LOEUF.score != "no") {
  LOEUF.score <- as.numeric(LOEUF.score)
}

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
#firehose.subtypes.path <- paths[paths$folder_or_object=="firehose_subtypes_path","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
transcripts.to.remove.path <- paths[paths$folder_or_object=="transcripts_to_remove_path","path_or_filename"]
nb.coeff.res.path <- paths[paths$folder_or_object=="NB_results","path_or_filename"]
CNV.path <- paths[paths$folder_or_object=="CNV_file_path","path_or_filename"]
ASE.somatic.path <- paths[paths$folder_or_object=="ASE_somatic_path","path_or_filename"]
transcripts.germ.mut.path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
somatic.PTCs.transcripts.path <- paths[paths$folder_or_object=="somatic_PTCs_transcripts_path","path_or_filename"]
TCGA.cancer.samples.metadata.path <- paths[paths$folder_or_object=="TCGA_cancer_samples_metadata_path","path_or_filename"]

print(paste0("TCGA cancer --> ",TCGA.cancer))
print("nextflow parameters:")
print(nextflow.params)

# 1.1) Protein coding transcripts
ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                  header = FALSE, sep = "\t")   
# 1.2) LOEUF score (negative selected genes)
# This is to remove constrained transcripts, i.e. those with potentially negative selection. Why? Because if not, dosage compensation from the other wild-type parental allele might happen and can mask the NMD effect
# If the transcript is not negatively selected, no need for dosage compensation!
LOEUF.score.table <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="LOEUF_score","path_or_filename"]), header = TRUE, sep = "\t")
# oe_lof_upper --> LOEUF: upper bound of 90% confidence interval for o/e ratio for pLoF variants (lower values indicate more constrained)
# oe_lof_upper_bin --> Decile bin of LOEUF for given transcript (lower values indicate more constrained)
LOEUF.score.ensembl <- LOEUF.score.table[,c("transcript","oe_lof_upper","oe_lof_upper_bin")]
LOEUF.transcripts.to.remove <- na.omit(LOEUF.score.ensembl[which(LOEUF.score.ensembl$oe_lof_upper_bin<=LOEUF.score),"transcript"])
# 1.3) Positive selected genes (PSG)
PSG.ensembl.gene <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="PSG_331","path_or_filename"]),
                                   header = TRUE, sep = "\t")
# 1.4) NMD targets (for negative control)
NMD.global.ensembl <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_ensembl","path_or_filename"]),
                                   header = TRUE, sep = "\t")
# 1.5) transcripts filtering matrix with somatic AND CNV
transcripts.to.remove <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(transcripts.to.remove.path,paths[paths$folder_or_object=="transcripts_to_remove","path_or_filename"])),
                                    header = TRUE, sep = "\t", row.names = 1)

# 2) Filterings
# Samples with bad purity/LF or without CNV file
TCGA.cancer.samples.metadata <- read.table(file = gsub("\\[X1\\]",TCGA.cancer, paste0(TCGA.cancer.samples.metadata.path,paths[paths$folder_or_object=="TCGA_cancer_samples_metadata","path_or_filename"])), sep ="\t", header = TRUE)
filter <- (TCGA.cancer.samples.metadata$LF_remove == "yes") | (TCGA.cancer.samples.metadata$purity_remove == "yes")
samples.to.remove <- TCGA.cancer.samples.metadata[filter,"sample"]

# 3) Negative Binomial regression on ASE data
columns <- c("sample","subtype","purity","LF","num.PTCs","stopgain","stopgain.NMD.evading",
            "stopgain.pval","synonymous","synonymous.pval","error")
nb.coeff.res <- data.frame(matrix(ncol=length(columns),nrow=0, dimnames=list(NULL, columns)))
PTC.set.names <- c("stopgain","synonymous")
total.num.PTC.sample <- list()

# Obtain samples names with available ASE somatic data
ASE.somatic.path.cancer <- gsub("\\[X1\\]",TCGA.cancer,ASE.somatic.path)
samples <- list.dirs(ASE.somatic.path.cancer, full.names = FALSE, recursive = FALSE)
samples <- gsub("-",".", samples)

for (sample in seq(1,length(samples))) {

  start_time_total <- Sys.time()
  TCGA.barcode <- samples[sample]
  nb.coeff.res[sample,"sample"] <- TCGA.barcode
  print(paste("SAMPLE", sample, " ----> ",TCGA.barcode),sep = "")
  purity <- TCGA.cancer.samples.metadata[TCGA.cancer.samples.metadata$sample%in%TCGA.barcode,"purity"]
  LF <- TCGA.cancer.samples.metadata[TCGA.cancer.samples.metadata$sample%in%TCGA.barcode,"LF"]
  # Check LF and purity
  if (TCGA.barcode%in%samples.to.remove) {
    print("High Leukocyte Fraction / Low tumor purity sample / no CNA file, skipping...")
    nb.coeff.res[sample,"error"] <- "High Leukocyte Fraction / Low tumor purity / no CNA file"
    #next
  }
  if (TCGA.cancer == "TCGA-THYM") {
    nb.coeff.res[sample,"purity"] <- NA
    purity <- NULL
    LF <- NA
  } else {
    nb.coeff.res[sample,"purity"] <- mean(purity)
    nb.coeff.res[sample,"LF"] <- mean(LF)
  }

  # 1) Obtain somatic PTCs (from nonsense + FS) with predicted NMD-rules for the sample
  somatic.PTCs.transcripts.cancer.path <- gsub("\\[X1\\]",TCGA.cancer, paste0(somatic.PTCs.transcripts.path,paths[paths$folder_or_object=="somatic_PTCs_transcripts","path_or_filename"]))
  somatic.PTCs.transcripts.sample.path <- gsub("\\[X2\\]",gsub("\\.","-",TCGA.barcode), somatic.PTCs.transcripts.cancer.path)
  if (file.exists(somatic.PTCs.transcripts.sample.path)) {
    somatic.PTCs.transcripts.sample <- read.table(file = somatic.PTCs.transcripts.sample.path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    print("somatic PTCs file do not exist, skipping...")
    next
  }
  # 2) Obtain ASE somatic data file
  # Check if file exists
  ASE.somatic.sample <- NULL
  ASE.somatic.sample <- ASE_somatic_sample_file(TCGA.barcode = TCGA.barcode, TCGA.cancer = TCGA.cancer, method = method)
  if(is.null(ASE.somatic.sample)){
    nb.coeff.res[sample,"error"] <- "ASE somatic file not found"
    print("ASE somatic file not found, skipping...")
    next
  }

  # 3) NB regression per NMD set
  for (PTC.set.name in PTC.set.names) {
    for (type in c("NMD-triggering","NMD-evading")) {
    #for (type in c("NMD-triggering")) {
      if (PTC.set.name == "synonymous" & type == "NMD-evading") {next}
      # Init
      start_time <- Sys.time()
      nb.res <- NA
      print(paste0("------------------------",toupper(PTC.set.name),"------------------------------"))
      print(paste0("------------------------",toupper(type),"------------------------------"))
      # 4) Filter VCF for variants
      if (PTC.set.name == "stopgain") {
        ASE.somatic.sample.filt <- ASE_somatic_sample_filtering(ASE.somatic.sample = ASE.somatic.sample, method = method, 
                                    mutation = PTC.set.name, VAF = 1, filter = "PASS", coverage = 0)
        if(nrow(ASE.somatic.sample.filt)==0){
          print("No variants left, skipping...")
          next
        }
        # Change ASE columns
        ASE.somatic.sample.save <- ASE.somatic.sample.filt[,c("CHROM","START","REF_DNA","ALT_DNA","refCount_RNA","altCount_RNA","totalCount_RNA")]
        colnames(ASE.somatic.sample.save) <- c("CHROM","START","REF_DNA","ALT_DNA","refCount_ASE","altCount_ASE","totalCount_ASE")
        # Merge with PTCs metadata
        somatic.PTCs.transcripts.sample.ASE <- merge(somatic.PTCs.transcripts.sample,ASE.somatic.sample.save, by.x = c("chr","start_pos","Ref","Alt"), 
                by.y = c("CHROM","START","REF_DNA","ALT_DNA"), all.x = TRUE)
        # Save raw PTCs_ASE df for sample
        somatic.PTCs.ASE.transcripts.cancer.path <- gsub("\\[X1\\]",TCGA.cancer, paste0(somatic.PTCs.transcripts.path,paths[paths$folder_or_object=="somatic_PTCs_ASE_transcripts","path_or_filename"]))
        somatic.PTCs.ASE.transcripts.sample.path <- gsub("\\[X2\\]",gsub("\\.","-",TCGA.barcode), somatic.PTCs.ASE.transcripts.cancer.path)
        #if ( (!file.exists(somatic.PTCs.ASE.transcripts.sample.path)) & ( PTC.set.name == "stopgain" & type == "NMD-triggering"  & VAF == 0.01) ) {
        if ( ( PTC.set.name == "stopgain" & type == "NMD-triggering" & VAF == 0.2) ) {
          write.table(somatic.PTCs.transcripts.sample.ASE, file = somatic.PTCs.ASE.transcripts.sample.path , 
                      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
        }
        
        # 5) Filterings for the NB regression
        # 5.1) Coverage 5 for nonsense and 2 for FS (filter for ASE PTCs)
        filter.coverage.nonsense <- which(somatic.PTCs.transcripts.sample.ASE$totalCount_ASE >= 5 & somatic.PTCs.transcripts.sample.ASE$stopgain == "nonsense")
        filter.coverage.FS <- which(somatic.PTCs.transcripts.sample.ASE$totalCount_ASE >= 2 & somatic.PTCs.transcripts.sample.ASE$stopgain != "nonsense")
        somatic.PTCs.transcripts.sample.ASE.filt <- somatic.PTCs.transcripts.sample.ASE[c(filter.coverage.nonsense,filter.coverage.FS),]
        # 5.2) Remove ASE-PTCs from FS with no predicted PTC
        if (length(filter.coverage.FS) != 0) {
          filter.FS.no.PTC <- which(somatic.PTCs.transcripts.sample.ASE.filt$stopgain != "nonsense" & !is.na(somatic.PTCs.transcripts.sample.ASE.filt$PTC_CDS_pos)) 
        } else {
          filter.FS.no.PTC <- NULL
        }
        # 5.3) Remove ASE-PTCs with specified VAF
        filter.VAF <- which(somatic.PTCs.transcripts.sample.ASE.filt$VAF > VAF)
        # 5.4) Remove transcripts with 1 exon
        filter.single.exon <- which(somatic.PTCs.transcripts.sample.ASE.filt$transcript_CDS_exon_num == 1 & somatic.PTCs.transcripts.sample.ASE.filt$splice_site_3UTR == "no")
        # 5.5) LOEUF score 1
        filter.LOEUF <- which(somatic.PTCs.transcripts.sample.ASE.filt$LOEUF_decile <= LOEUF.score)
        # 5.6) PSG genes
        filter.PSG <- which(somatic.PTCs.transcripts.sample.ASE.filt$PSG == "yes")
        # 5.6) Transcripts overlapping somatic mut/CNV
        #filter.somatic.mut.CNV <- which(somatic.PTCs.transcripts.sample.ASE.filt$somatic_CNV_SNV == "yes")
        # 5.7) Coefficient variation
        # Coefficient of variation (mean/SD) > 0.5 --> This reduces a lot the NMDeff values (from 1.6 raw to 0.3)
        #filter.coeff.var <- which(somatic.PTCs.transcripts.sample.ASE.filt$coeff_var > 0.5)
        # 5.8) NMD-triggering vs NMD-evading
        if (type == "NMD-triggering") {
          filter.NMDtype <- which(somatic.PTCs.transcripts.sample.ASE.filt$X55_nt_last_exon == "NMD-evading")
        } else if (type == "NMD-evading") {
          filter.NMDtype <- which(somatic.PTCs.transcripts.sample.ASE.filt$X55_nt_last_exon == "NMD-triggering")
        }
        # 5.9) Non protein-coding
        filter.nonCoding <- which(somatic.PTCs.transcripts.sample.ASE.filt$protein_coding == "no")
        # 5.10) Overlapping somatic mutations...
        # Add up filters
        all.filters <- unique(c(filter.FS.no.PTC,filter.VAF,filter.single.exon,filter.LOEUF,filter.NMDtype,filter.nonCoding,filter.PSG))
        somatic.PTCs.transcripts.sample.ASE.filt <- somatic.PTCs.transcripts.sample.ASE.filt[-all.filters,]
        if (nrow(somatic.PTCs.transcripts.sample.ASE.filt) == 0 ){
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" )  {
            nb.coeff.res[sample,"error"] <- "<1 total genes"
          }
          next
        } else {
          NMD.final.df <- somatic.PTCs.transcripts.sample.ASE.filt[,c("start_pos","gene_id","refCount_ASE","altCount_ASE","VAF")]
          colnames(NMD.final.df) <- c("start_pos","gene_id","refCount","altCount","AF")
          # Remove redundant transcripts with same ASE counts
          NMD.final.df <- NMD.final.df[which(!duplicated(NMD.final.df)),]
        }
      } else if ( PTC.set.name == "synonymous") {
        ASE.somatic.sample.filt <- ASE_somatic_sample_filtering(ASE.somatic.sample = ASE.somatic.sample, method = method, 
                                    mutation = PTC.set.name, VAF = VAF, filter = "PASS", coverage = 5)
        if(nrow(ASE.somatic.sample.filt)==0){
          print("No variants left, skipping...")
          next
        }
        variants <- ASE_extract_info(ASE = ASE.somatic.sample.filt)
        somatic.PTCs.transcripts.sample.ASE <- do.call(rbind,variants)
        # 5) Filterings for the NB regression 
        # 5.1) VAF and coverage (done)
        # 5.2) LOEUF score
        filter.LOEUF <- which(somatic.PTCs.transcripts.sample.ASE$transcript_id %in% LOEUF.transcripts.to.remove)
        # 5.3) Non protein-coding
        genes.id <- strsplit(somatic.PTCs.transcripts.sample.ASE$gene_id,";")
        genes.id.df <- as.data.frame(do.call(rbind,genes.id))
        filter.nonCoding <- which( !as.character(genes.id.df$V1)%in%ensembl.v88.coding.transcripts$V2 | !as.character(genes.id.df$V2)%in%ensembl.v88.coding.transcripts$V2 )
        # 5.3) PSG genes
        filter.PSG <- which( as.character(genes.id.df$V1)%in%PSG.ensembl.gene$Ensembl.Gene.ID | as.character(genes.id.df$V2)%in%PSG.ensembl.gene$Ensembl.Gene.ID )
        # 5.4) NMD genes
        NMD.genes <- unique(somatic.PTCs.transcripts.sample.ASE[as.character(somatic.PTCs.transcripts.sample.ASE$transcript_id)%in%NMD.global.ensembl$ensembl_transcript_id,"gene_id"])
        filter.NMD <- which(as.character(somatic.PTCs.transcripts.sample.ASE$gene_id)%in%NMD.genes)
        # 5.5) Transcripts with overlapping somatic mut/CNV        
        filter.somatic.mut.CNV <- somatic_transcripts_filtering(ASE.somatic = somatic.PTCs.transcripts.sample.ASE, transcripts.filtering.matrix = transcripts.to.remove, TCGA.barcode = TCGA.barcode)
        # Add up filters
        all.filters <- unique(c(filter.LOEUF,filter.NMD,filter.nonCoding,filter.somatic.mut.CNV,filter.PSG))
        somatic.PTCs.transcripts.sample.ASE.filt <- somatic.PTCs.transcripts.sample.ASE[-all.filters,]
        # 5.6) Remove redundant transcripts with same ASE counts
        somatic.PTCs.transcripts.sample.ASE.filt <- somatic.PTCs.transcripts.sample.ASE.filt[,c("start_pos","gene_id","refCount","altCount","AF")]
        NMD.final.df <- somatic.PTCs.transcripts.sample.ASE.filt[which(!duplicated(somatic.PTCs.transcripts.sample.ASE.filt)),]
        if (nrow(NMD.final.df) == 0 ){
          if ( PTC.set.name == "stopgain" & type == "NMD-triggering" )  {
            nb.coeff.res[sample,"error"] <- "<1 total genes"
          }
          next
        }
        # 5.7) Sample 100 variants maximum
        if ( nrow(NMD.final.df) > 100) { 
          NMD.final.df <- NMD.final.df[sample(1:nrow(NMD.final.df))[1:100],]
        }
      }
      # Number of PTC genes for the regression
      if ( PTC.set.name == "stopgain" & type == "NMD-triggering" ) {
        num.PTC <- nrow(NMD.final.df)
        nb.coeff.res[sample,"num.PTCs"] <- num.PTC
        print(paste0("Final number of PTCs --> ", num.PTC))
        total.num.PTC.sample[[TCGA.barcode]] <- num.PTC
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

file.path <- gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(nb.coeff.res.path,paths[paths$folder_or_object=="NB_coeff_res","path_or_filename"]))
file.path <- gsub("\\[X2\\]",paste0(as.character(VAF),"-",method), file.path)
file.path <- gsub("ASE_NB_coeff","ASE_somatic_NB_coeff", file.path)
write.table(nb.coeff.res, file = file.path , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

