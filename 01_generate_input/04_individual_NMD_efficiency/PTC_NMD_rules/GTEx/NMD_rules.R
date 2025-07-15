rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

VCF_extract_info <- function(VCF, VCF.type) {

    if (VCF.type == "somatic") {
      VCF.tmp <- VCF
      colnames(VCF.tmp)[colnames(VCF.tmp) %in% c("Alt","Otherinfo12","Otherinfo14")] <- c("ALT","FORMAT","TUMOR")
      VCF.tmp <- GetStrelkaVAF(VCF.tmp, name.of.VCF = NULL)
    }

    variants.list <- list()
    for (i in 1:nrow(VCF)) {
      # Transcripts list from same variant
      transcripts.list <- strsplit(as.character(VCF[i,"AAChange.refGene"]), split=",")
      # Obtain the ID name and other info (exon location, PTC CDS location)
      transcripts.exons.info <- lapply(transcripts.list, function(transcript) {
        variant.CDS.exon.num <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon([0-9]{1,3}).*","\\1",transcript)
        transcript.id <- gsub(".*(ENST[0-9]{11})\\.[0-9]{1,3}:exon.*","\\1",transcript)
        #variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)\\D*:p\\..*","\\1",transcript)
        variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)_?.*:p\\..*","\\1",transcript)
        data.frame(transcript_id=transcript.id,variant_exon_num=variant.CDS.exon.num, variant_CDS_pos = variant.CDS.pos)
      })
      transcripts.exons.info <- transcripts.exons.info[[1]]
      transcripts.exons.info$start_pos <- VCF[i,"Start"]
      transcripts.exons.info$gene_id <- VCF[i,"Gene.refGene"]
      transcripts.exons.info$chr <- VCF[i,"Chr"]
      transcripts.exons.info$genotype <- unlist(strsplit(VCF[i,"Otherinfo13"],":"))[1]
      # Type of stopgain (nonsense or frameshift)
      ref.allele <- as.character(VCF[i,"Ref"])
      mut.allele <- as.character(VCF[i,"Alt"])
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
      if ( (ref.allele.num == 0) & (mut.allele.num >= 1) ) {
        transcripts.exons.info$stopgain <- "frameshift_insertion"
      } else if ( (ref.allele.num >= 1) & (mut.allele.num == 0) ) {
        transcripts.exons.info$stopgain <- "frameshift_deletion"
      } else {
        transcripts.exons.info$stopgain <- "nonsense"
      }
      # Predicted PTC? For FS indels
      transcripts.exons.info$variant_type <- gsub(" ","_",VCF[i,"ExonicFunc.refGene"])
      # Calculate VAF
      if (VCF.type == "germline") {
        # Population MAF
        transcripts.exons.info$VAF <- VCF[i,"AF"]
      } else if (VCF.type == "somatic") {
        # VAF within sample tumor estimation
        transcripts.exons.info$VAF <- VCF.tmp[i,"VAF"]
      }
      variants.list[[i]] <- transcripts.exons.info
  }
  return(variants.list)
}

NMD_rules <- function(VCF, VCF.type) {

  # 1) Obtain ENSEMBL transcript ID(s) and some info for each variant
  # Note: there can be >=1 transcript overlapping the same variant. We will take them all
  PTCs.list <- VCF_extract_info(VCF = VCF, VCF.type = VCF.type)
  # 2) Check NMD rules for each PTC-transcript combination
  # Variants lists from a PTC
  PTC.variants.all <- c()
  #count <- 0
  columns <- c("variant_CDS_exon_num","variant_CDS_exon_length","PTC_CDS_pos","PTC_CDS_exon_num","PTC_CDS_exon_length",
                "PTC_EJC_dist","TSS_PTC_dist","CDS_num_exons_downstream","55_nt", "PTC_stop_codon_type",
                "ATG_downstream","intercistronic_dist","ATG_number","ATG_kozak_score_A","ATG_kozak_score_B","PTC_downstream_4nt",
                "UTR5s_length","introns_length_prePTC", "exons_length_prePTC","introns_length_postPTC","exons_length_postPTC",
                "UTR3s_length","normal_stop_codon_CDS_pos","transcript_GC_content","ATG_CDS_pos")
  #PTC.variants <- PTCs.list[[2]]
  PTCs.list.NMD <- lapply(PTCs.list, function(PTC.variants) {
    #count <<- count + 1
    # Each variant
    for (i in 1:nrow(PTC.variants)) {
        #print(i)
        PTC <- PTC.variants[i,]
        PTC[,columns] <- NA
        PTC[,"comments"] <- ""
        # 2.1) Obtain PTC position using coordinates from the ENSEMBL GTF file
        # Filter ENSEMBL GTF file
        ensembl.v88.gtf.filt <- ensembl.v88.gtf[ensembl.v88.gtf$transcript_id%in%PTC$transcript_id,]
        # Check Chr is the same (some transcript might be in X and Y chr)
        ensembl.v88.gtf.filt <- ensembl.v88.gtf.filt[which(ensembl.v88.gtf.filt$seqnames%in%PTC$chr),]
        # TSS (Translation Start Site) or start codon
        TSS.abs.pos <- as.numeric(ensembl.v88.gtf.filt[ensembl.v88.gtf.filt$type == "start_codon","start"])
        # Strand
        strand <- as.character(unique(ensembl.v88.gtf.filt$strand))
        PTC[,"strand"] <- strand
        # Start / Stop codon
        start.codon.start <- ensembl.v88.gtf.filt[ensembl.v88.gtf.filt$type == "start_codon","start"][1]
        stop.codon <- ensembl.v88.gtf.filt[ensembl.v88.gtf.filt$type == "stop_codon",c("start","end")]
        # UTRs
        UTRs <- ensembl.v88.gtf.filt[ensembl.v88.gtf.filt$type == "UTR",c("start","end")]
        # Obtain 3UTRs and 5UTRs
        if (strand == "+") {
          UTR5 <- UTRs[UTRs$start <= start.codon.start,]
          UTR3 <- UTRs[UTRs$end >= stop.codon$end,]
        } else if (strand == "-") {
          UTR5 <- UTRs[UTRs$start >= start.codon.start,]
          UTR3 <- UTRs[UTRs$end <= stop.codon$end,]
        }
        # Alleles
        ref.allele <- tolower(as.character(PTC[,"Ref"]))
        mut.allele <- tolower(as.character(PTC[,"Alt"]))
        if (mut.allele == "-") {
          mut.allele <- ""
        }
        if (ref.allele == "-"){
          ref.allele <- ""
        }
        ref.allele.num <- nchar(ref.allele)
        mut.allele.num <- nchar(mut.allele)
        # If start or stop codons don't exist, then initialize UTRs to not enter in its conditions
        if (nrow(stop.codon)==0 || is.na(stop.codon)) {
          UTR3 <- data.frame()
          PTC[,"comments"] <- paste0(PTC[,"comments"],"/no normal stop codon")
        } 
        if (length(start.codon.start)==0 || is.na(start.codon.start)) {
          UTR5 <- data.frame()
          PTC[,"comments"] <- paste0(PTC[,"comments"],"/no start codon")
        }
        PTC[,"UTR3"] <- nrow(UTR3)
        PTC[,"UTR5"] <- nrow(UTR5)
        # CDS
        ensembl.v88.gtf.filt.cds  <- ensembl.v88.gtf.filt[ensembl.v88.gtf.filt$type=="CDS",]
        # Order
        ensembl.v88.gtf.filt.cds <- ensembl.v88.gtf.filt.cds[order(ensembl.v88.gtf.filt.cds$start),]
        # Introns, exons size
        introns.size <- (ensembl.v88.gtf.filt.cds$start[-1] -1) - (ensembl.v88.gtf.filt.cds$end[-length(ensembl.v88.gtf.filt.cds$end)])
        exons.size <- ensembl.v88.gtf.filt.cds$end-ensembl.v88.gtf.filt.cds$start + 1
        # Position of the variant relative since first position of first exon or TSS (theoretically it takes into account exon part including UTRs + introns)
        #variant.gene.position <- PTC$start_pos-TSS.abs.pos
        variant.gene.position <- PTC$start_pos-ensembl.v88.gtf.filt.cds$start[1]+1
        # Check in which CDS exon the variant is located (+ strand)
        variant.CDS.exon.num <- which(PTC$start_pos >= ensembl.v88.gtf.filt.cds$start & PTC$start_pos <= ensembl.v88.gtf.filt.cds$end)
        if(length(variant.CDS.exon.num) == 0 ){
          print("Coords are wrong, because exon number is not found")
          next
        }
        # Position of the variant in the CDS (only exons, removing introns size)
        # Variant position of the first exon
        if (variant.CDS.exon.num == 1) {
          variant.cds.position <- variant.gene.position
        } else { # Otherwise
          variant.cds.position <- variant.gene.position-sum(introns.size[1:(variant.CDS.exon.num-1)])
        }
        if (strand == "-") {
          variant.cds.position <- sum(exons.size)-variant.cds.position+1
          variant.CDS.exon.num <- length(exons.size)-variant.CDS.exon.num+1
        }
        # If it is a frameshift stopgain, then we slightly modify the variant CDS position, depending on strand and type of FS
        if ( (PTC[,"stopgain"] == "frameshift_insertion") & (strand == "-") ) {
          variant.cds.position <- variant.cds.position - 1
        } else if ( (PTC[,"stopgain"] == "frameshift_insertion") & (strand == "+") & (nchar(PTC$Alt) == 1) ) {         
          variant.cds.position <- variant.cds.position + 1
        } else if ( (PTC[,"stopgain"] == "frameshift_deletion") & (strand == "-") & (nchar(PTC$Ref) >= 2) ) {
          variant.cds.position <- variant.cds.position - (nchar(PTC$Ref) - 1)
          if (variant.cds.position < 0) {
            variant.cds.position <- 1
          }
        }
        # We have to consider the total number of exons (those that also includes UTR) for the absolute exon number location of the variant
        total.exon.num <- max(as.numeric(na.omit(ensembl.v88.gtf.filt$exon_number)))
        total.cds.exon.num <- length(exons.size)
        PTC[,"transcript_CDS_exon_num"] <- total.cds.exon.num
        PTC[,"transcript_total_exon_num"] <- total.exon.num
        ### DEPRECATED ### This is not correct 
        #PTC.exon.num <- PTC.CDS.exon.num +  (total.exon.num - total.cds.exon.num)
        ##################
        # We need to sum only 3UTRs or 5UTRs exons depending if it is (+) or (-) strand, but some UTRs overlap with first or last exon,
        # so lets take directly the exon number from the GTF file
        tmp.rows <- which(PTC$start_pos >= ensembl.v88.gtf.filt$start & PTC$start_pos <= ensembl.v88.gtf.filt$end)
        variant.exon.num <- as.numeric(unique(na.omit(ensembl.v88.gtf.filt[tmp.rows,"exon_number"])))
        # Check that PTC.cds.pos and PTC.CDS.exon.num from ENSEMBL GTF coincides with VCF file
        if (PTC$variant_exon_num != variant.exon.num ) {
          print("Exon PTC location wrong!")
          PTC[,"comments"] <- paste0(PTC[,"comments"],"/Exon PTC location wrong!")
          print(paste0("------------------ ENSEMBL GTF --> ",variant.exon.num ))
          print(paste0("------------------ VCF ----------> ",PTC$variant_exon_num))
        }
        if (as.character(PTC$variant_CDS_pos) != as.character(variant.cds.position)) {
          print("variant CDS position wrong!")
          PTC[,"comments"] <- paste0(PTC[,"comments"],"/PTC CDS position wrong!")
          print(paste0("------------------ ENSEMBL GTF --> ",variant.cds.position))
          print(paste0("------------------ VCF ----------> ",PTC$variant_CDS_pos))
        }

        # 3) Check NMD_rules for the variant
        # 3.1) Obtain mutated fasta sequence and the PTC CDS position using the variant CDS position
        PTC <- mutate_fasta_sequence(PTC=PTC)
        cols <- c("CDS_mut_length","exon_length_with_UTR3","PTC_EJC_dist","splice_site_3UTR","PTC_CDS_exon_num","PTC_CDS_exon_length","TSS_PTC_dist",
                "CDS_num_exons_downstream","last_exon","55_nt","55_nt_last_exon","PTC_downstream_4nt","ATG_downstream",
                "ATG_number","ATG_kozak_score_A","ATG_kozak_score_B","intercistronic_dist","UTR5s_length","introns_length_prePTC",
                "exons_length_prePTC","introns_length_postPTC","exons_length_postPTC","UTR3s_length","normal_stop_codon_CDS_pos",
                "transcript_GC_content","ATG_CDS_pos")
        if (is.na(PTC[,"PTC_CDS_pos"])) {
          PTC[,cols] <- NA 
          if (is.null(PTC.variants.all)) {
            PTC.variants.all <- PTC
          } else {
            PTC.variants.all <- rbind(PTC.variants.all,PTC)
          }   
          next
        } else {
          PTC.CDS.pos <- PTC[,"PTC_CDS_pos"]
          # 3.2) Obtain some PTC info (exon length, last exon CDS pos, etc)
          find_exon_position <- function(num, vec) {
            # Initialize a variable to store the sum
            total <- 0
            # Iterate through the elements in the vector
            for (i in 1:length(vec)) {
              total <- total + vec[i]
              if (num <= total) {
                # Return the position if the number is lower than the sum
                return(i)
              }
            }
              # Return -1 if the number is not lower than any sum
              return(-1)
          }
          modify_exons_size <- function(pos,exons) {
            variant.CDS.exon.length <- exons[pos]
            if ((PTC[,"stopgain"] == "frameshift_insertion")) {
              exons[pos] <- variant.CDS.exon.length + mut.allele.num
            } else if (PTC[,"stopgain"] == "frameshift_deletion") {
              exons[pos] <- variant.CDS.exon.length - ref.allele.num
            }
            return(exons)
          }

          if (strand == "-") {
            exons.size <- rev(exons.size)
          }
          # Add or substract nucleotides to the exon where variant is located (not PTC) if the variant is INS or DEL
          exons.size.modified <- modify_exons_size(variant.CDS.exon.num,exons.size)
          # variant exon length
          variant.CDS.exon.length <- exons.size.modified[variant.CDS.exon.num]
          # Last exon CDS start position
          last.exon.cds.start.pos <- sum(exons.size.modified)-exons.size.modified[length(exons.size.modified)]
          # PTC exon number 
          # If PTC CDS pos is larger than length of exon there variant is located, PTC is located in at least a downstream exon!
          if (PTC.CDS.pos > variant.CDS.exon.length) {
            PTC.CDS.exon.num <- find_exon_position(PTC.CDS.pos,exons.size.modified)
          } else {
            PTC.CDS.exon.num <- variant.CDS.exon.num
          }
          # Sometimes transcript sequence is longer than sum of exons (errors in obtaining fasta seq...)
          if ( PTC.CDS.exon.num == -1 ) {
            PTC[,cols] <- NA 
            if (is.null(PTC.variants.all)) {
              PTC.variants.all <- PTC
            } else {
              PTC.variants.all <- rbind(PTC.variants.all,PTC)
            }   
            next
          }
          # PTC exon length
          PTC.CDS.exon.length <- exons.size.modified[PTC.CDS.exon.num]
          # CDS end position of the exon where the PTC is located
          PTC.exon.cds.end.pos <- sum(exons.size.modified[1:PTC.CDS.exon.num])

          # Is the PTC on the last exon?
          isLastExon <- (PTC.CDS.exon.num == total.cds.exon.num)
          # If PTC is on last exon then sum +3 nt of the stop codon
          if (isTRUE(isLastExon)) {
            PTC.exon.cds.end.pos <- PTC.exon.cds.end.pos + 3
            PTC.CDS.exon.length <- PTC.CDS.exon.length +3
          }
          # Control
          if ( (length(start.codon.start)==0 || is.na(start.codon.start)) || (nrow(stop.codon)==0 || is.na(stop.codon)) )  {
            next
          }
          # +3 because normal stop codon is not on the fasta_sequence_mut
          PTC[,"CDS_mut_length"] <- nchar(PTC$fasta_sequence_mut)+3
          # 3.2) PTC_EJC_dist
          # Distance between PTC and downstream EJC (i.e CDS end position of the exon where the PTC is located)
          # If #UTR3' >= 2 AND PTC pos is in last exon or <55 in penultime exon
          isPenultimateExonEvading  <- FALSE
          PTC[,"exon_length_with_UTR3"] <- NA
          if (PTC.CDS.exon.num == (total.cds.exon.num-1)) {
            isPenultimateExonEvading <- abs(PTC.CDS.pos - PTC.exon.cds.end.pos) <= 55
          }
          if ( (nrow(UTR3) >= 2) & ( isTRUE(isLastExon) || isTRUE(isPenultimateExonEvading) ) ) {
            PTC[,"comments"] <- paste0(PTC[,"comments"],"/splice site UTR3 and PTC at last exon")
            splice_site_UTR3 <- TRUE
            # downstream EJC position is the splice site found in 3'UTR
            UTR3.length <- UTR3[1,"end"]-UTR3[1,"start"] - 1
            if (strand == "-") {
              #UTR3.length <- abs(UTR3[nrow(UTR3),"start"]-UTR3[nrow(UTR3),"end"])
              UTR3.length <- abs(UTR3[1,"start"]-UTR3[1,"end"]) - 1
            }
            PTC[,"exon_length_with_UTR3"] <- PTC.CDS.exon.length + UTR3.length - 1
            EJC.cds.pos <- PTC.exon.cds.end.pos + UTR3.length
          } else {
            splice_site_UTR3 <- FALSE
            # downstream EJC position is the last CDS base of the exon where the PTC is located
            EJC.cds.pos <- PTC.exon.cds.end.pos + 1
          }
          # NOT +2 because it's from START of PTC till EJC (not the END)
          PTC.EJC.dist <- abs(PTC.CDS.pos - EJC.cds.pos)
          PTC[,"PTC_EJC_dist"] <- PTC.EJC.dist
          # 3.3) Splice site 3'UTR
          if (isTRUE(splice_site_UTR3)) {
            PTC[,"splice_site_3UTR"] <- "yes"
          } else {
            PTC[,"splice_site_3UTR"] <- "no"
          }
          # 3.4) CDS exon length and number (of where PTC and variant are located)
          PTC[,"PTC_CDS_exon_length"] <- PTC.CDS.exon.length
          PTC[,"variant_CDS_exon_length"] <- variant.CDS.exon.length
          PTC[,"PTC_CDS_exon_num"] <- PTC.CDS.exon.num
          PTC[,"variant_CDS_exon_num"] <- variant.CDS.exon.num
          # 3.5) TSS_PTC_dist
          PTC[,"TSS_PTC_dist"] <- PTC.CDS.pos 
          # 3.6) Number of CDS exons downstream of PTC
          PTC[,"CDS_num_exons_downstream"] <- total.cds.exon.num - PTC.CDS.exon.num
          # 3.7) Last_exon rule (without taking into account splice site 3'UTR)
          if (isTRUE(isLastExon)) {
            PTC[,"last_exon"] <- "yes"
          } else {
            PTC[,"last_exon"] <- "no"
          }
          # 3.8) 55_nt rule (PTCs located only in the 55nt boundary of the penultimate exon, without taking into account splice site 3'UTR)
          if(isTRUE(isPenultimateExonEvading)) {
            PTC[,"55_nt"] <- "yes"
          } else {
            PTC[,"55_nt"] <- "no"
          }
          # 3.9) 55_nt + last exon rule
          #  PTC located upstream or downstream of 55nt from last base of last exon (taking into account splice site at 3UTR)
          if (isTRUE(splice_site_UTR3)) {
            if (PTC.EJC.dist > 55) {
              PTC[,"55_nt_last_exon"] <- "NMD-triggering"
            } else if (PTC.EJC.dist <= 55) {
              PTC[,"55_nt_last_exon"] <- "NMD-evading"
            }
          } else if (!isTRUE(splice_site_UTR3)) {
            if (isTRUE(isLastExon)) {
              PTC[,"55_nt_last_exon"] <- "NMD-evading"
            } else if (isTRUE(isPenultimateExonEvading)) {
              PTC[,"55_nt_last_exon"] <- "NMD-evading"
            } else if (!isTRUE(isPenultimateExonEvading)) {
              PTC[,"55_nt_last_exon"] <- "NMD-triggering"
            }
          }
          # If gene has 1 exon
          # if ( total.cds.exon.num == 1) {
          #   PTC[,"55_nt_last_exon"] <-"NMD-evading"
          # } 
          # 3.10) Get Stop Codon type (already in the PTC dataframe)
          # 3.11) ATG_downstream
          # Presence or absence of a inframe ATG downstream of the PTC
          #system(paste0("rm -rf /scratch/2/gpalou/trash/*"))
          ATG.df.res <- ATG_finder(PTC_CDS_pos = PTC.CDS.pos, fasta_sequence = PTC[,"fasta_sequence_mut"])
          # 3.12) Get 4 nt downstream of PTC
          PTC[,"PTC_downstream_4nt"] <- ATG.df.res[,"PTC_downstream_4nt"]
          PTC[,"ATG_downstream"] <- ATG.df.res[,"ATG_downstream"]
          # 3.13) Intercistronic_dist
          PTC[,"intercistronic_dist"] <- ATG.df.res[,"ATG_downstream_CDS_pos"]
          PTC[,"ATG_number"] <- ATG.df.res[,"ATG_number"]
          # 3.14) ATG_nt_context
          PTC[,"ATG_kozak_score_A"] <- ATG.df.res[,"kozak_score_A"]
          PTC[,"ATG_kozak_score_B"] <- ATG.df.res[,"kozak_score_B"]
          # 3.15) New features for Neural Network (Marcell)
          # Introns, exons size
          num.exons <- length(exons.size.modified)
          if (strand == "-") {
              introns.size <- rev(introns.size)
          }
          # Exons and introns size
          num.exons.prePTC <- num.exons - as.numeric(PTC["CDS_num_exons_downstream"]) - 1
          num.introns.prePTC <- num.exons.prePTC
          num.exons.postPTC <- as.numeric(PTC["CDS_num_exons_downstream"])
          num.introns.postPTC <- num.exons.postPTC
          # Check
          if ( !(num.exons.prePTC + 1 + num.exons.postPTC) == as.numeric(PTC["transcript_CDS_exon_num"]) ) {
              print("number of exons is wrong")
          }
          # 3.15.1) Pre-PTC exons/introns length
          if (num.exons.prePTC > 0) {
              exons.prePTC <- paste0(exons.size.modified[1:num.exons.prePTC], collapse =",")
              introns.prePTC <- paste0(introns.size[1:num.introns.prePTC], collapse =",")
              PTC["introns_length_prePTC"] <- paste0(introns.prePTC, collapse =",")
              PTC["exons_length_prePTC"] <- paste0(exons.prePTC, collapse =",")
          }
          # 3.15.2) Post-PTC exons/introns length
          if (num.exons.postPTC > 0) {
              exons.postPTC <- paste0(tail(exons.size.modified, num.exons.postPTC), collapse =",")
              introns.postPTC <- paste0(tail(introns.size, num.introns.postPTC), collapse =",")
              PTC["introns_length_postPTC"] <- paste0(introns.postPTC, collapse =",")
              PTC["exons_length_postPTC"] <- paste0(exons.postPTC, collapse =",")
          }
          # 3.15.3) UTRs length
          if (nrow(UTR3) != 0 ) {
              UTR3.size <- paste0((UTR3[,"end"]-UTR3[,"start"] - 1),collapse=",")
              if (strand == "-") {
                  UTR3.size <- rev(UTR3.size)
              }
              PTC["UTR3s_length"] <- UTR3.size
          }
          if (nrow(UTR5) != 0 ) {
              UTR5.size <- paste0((UTR5[,"end"]-UTR5[,"start"] - 1),collapse=",")
              if (strand == "-") {
                  UTR5.size <- rev(UTR5.size)
              }
              PTC["UTR5s_length"] <- UTR5.size
          }
          # 3.15.5) Normal stop codon CDS position 
          if (!is.na(PTC["CDS_mut_length"])) {
              PTC["normal_stop_codon_CDS_pos"] <- as.numeric(PTC["CDS_mut_length"]) - 3
          }
          # 3.15.6) GC content
          fasta.sequence.mut <- PTC["fasta_sequence_mut"]
          num_g <- str_count(fasta.sequence.mut, "g")
          num_c <- str_count(fasta.sequence.mut, "c")
          gc_content <- (num_g + num_c) / str_length(fasta.sequence.mut) * 100
          PTC["transcript_GC_content"] <- round(gc_content,2)
          # 3.15.7) ATG downstream position
          if (PTC["ATG_downstream"] == "yes") {
              PTC["ATG_CDS_pos"] <- as.numeric(PTC["intercistronic_dist"]) + as.numeric(PTC["PTC_CDS_pos"])
          }
        }
        if (is.null(PTC.variants.all)) {
          PTC.variants.all <- PTC
        } else {
          PTC.variants.all <- rbind(PTC.variants.all,PTC)
        }
    }
    PTC.variants.all
    })
  # Return a df
  return(PTCs.list.NMD)
}

mutate_fasta_sequence <- function(PTC) {

  transcript <- as.character(PTC[,"transcript_id"])
  # Obtain fasta sequence
  index <- grep(transcript,names(transcripts.fasta.sequences))
  transcript.fasta.sequence <- transcripts.fasta.sequences[[index]]
  # Alleles
  ref.allele <- tolower(as.character(PTC[,"Ref"]))
  mut.allele <- tolower(as.character(PTC[,"Alt"]))
  if (mut.allele == "-") {
    mut.allele <- ""
  }
  if (ref.allele == "-"){
    ref.allele <- ""
  }
  ref.allele.num <- nchar(ref.allele)
  mut.allele.num <- nchar(mut.allele)
  # Type stopgain
  stopgain.type <- as.character(PTC[,"stopgain"])
  # Variant CDS position
  variant.cds.position <- as.numeric(as.character(PTC[,"variant_CDS_pos"]))
  if (is.na(variant.cds.position)) {
    print("Variant CDS position weird")
    PTC[,"comments"] <- paste0(PTC[,"comments"],"/Variant CDS position weird")
    PTC[,"PTC_CDS_pos"] <- NA
    return(PTC)
  }
  # Wild-type fasta sequence
  fasta.sequence.wt <- tolower(unlist(strsplit(transcript.fasta.sequence,"")))
  if (PTC$strand == "-") {
    # Reverse complement the fasta sequence
    fasta.sequence.wt <- tolower(as.character(reverseComplement(x = DNAString(paste0(fasta.sequence.wt, collapse = "")))))
    fasta.sequence.wt <- unlist(strsplit(fasta.sequence.wt,""))
    # Reverse complement the Mutated Allele
    if (mut.allele != "") {
      mut.allele <- tolower(as.character(reverseComplement(x = DNAString(mut.allele))))
      mut.allele <- unlist(strsplit(mut.allele,"")) 
    }
    if (ref.allele != "") {
      ref.allele <- tolower(as.character(reverseComplement(x = DNAString(ref.allele))))
      ref.allele <- unlist(strsplit(ref.allele,"")) 
    }
  }
  # SANITY CHECK #1
  # Nucleotide positions 1-3 or 2-4 of wt fasta sequence must be ATG (TSS site)
  TSS.pos.1 <- paste0(fasta.sequence.wt[1:3], collapse = "")
  TSS.pos.2 <- paste0(fasta.sequence.wt[2:4], collapse = "")
  if ( TSS.pos.1 != "atg") {
    if (TSS.pos.2 != "atg") {
      PTC[,"comments"] <- paste0(PTC[,"comments"],"/No ATG found in fasta WT sequence")
    } else if (TSS.pos.2 == "atg") {
      # Remove first nucleotide so that first CDS position starts at ATG
      fasta.sequence.wt <- fasta.sequence.wt[-1]
    }
  }
  # Copy sequence for the mutated
  fasta.sequence.mutant <- fasta.sequence.wt
  # Mutated allele to add
  mut.allele.to.add <- tolower(unlist(strsplit(mut.allele,"")))
  # SANITY CHECK #2
  # Last 3 or 3-1 nucleotide positions of wt fasta sequence must be stop codon
  # It might not in (-) strand, because in some cases the obtained fasta sequence lacks 1-2 nt for unkown reason
  # That's not a problem in our case as we only care about PTC position, and the start CDS position that should start at ATG. The end is not a big deal.
  stop.codon <- paste0(fasta.sequence.mutant[(length(fasta.sequence.mutant)-2):length(fasta.sequence.mutant)],collapse="")
  stop.codon.match <- grep("taa|tag|tga", stop.codon, ignore.case = TRUE)
  if (length(stop.codon.match) != 1) {
    PTC[,"comments"] <- paste0(PTC[,"comments"],"/Last 3 nucleotides are not a stop codon!")
  }
  # SANITY CHECK #3
  # Check REF alleles match with fasta sequences (only for nonsense and FS deletions)
  if ( ( stopgain.type == "frameshift_deletion" | stopgain.type == "nonsense" ) ) {
    cond <- sum( na.omit( unlist(strsplit(ref.allele,"")) != fasta.sequence.wt[(variant.cds.position):(variant.cds.position+ref.allele.num-1)] ) )
    if ( cond != 0 ) {
      print("REF allele does not match with WT fasta sequence")
      PTC[,"comments"] <- paste0(PTC[,"comments"],"/REF allele does not match with WT fasta sequence")
    }
  }


  # Some FS Deletions can occur on the PTC CDS == 1 (but in reality it's on the 5'UTR overlapping with first CDS nt)
  # Even overlapping stopcodon...



  if( stopgain.type == "frameshift_insertion" ) {
    n1 <- 0
    n2 <- 1
  } else if ( stopgain.type == "frameshift_deletion" | stopgain.type == "nonsense" ) {
    n1 <- 1
    n2 <- 0
  }
  fasta.sequence.mutant <- c(fasta.sequence.wt[1:(variant.cds.position-n1)],  mut.allele.to.add, fasta.sequence.wt[length((fasta.sequence.wt[1:(variant.cds.position+ref.allele.num+n2)])):length(fasta.sequence.wt)])
  # Manual Check
  # n <- 15
  # ref.allele
  # mut.allele
  #fasta.sequence.wt[(variant.cds.position-n):(variant.cds.position+n)]
  #fasta.sequence.mutant[(variant.cds.position-n):(variant.cds.position+n)]
  # Save fasta sequences
  PTC[,"fasta_sequence_wt"] <- paste0(fasta.sequence.wt,collapse="")
  PTC[,"fasta_sequence_mut"] <- paste0(fasta.sequence.mutant,collapse="")
  # Check ORFs
  ORFs <- findORFs(
    tolower(paste0(fasta.sequence.mutant,collapse="")),
    startCodon = "atg",
    stopCodon = "taa|tag|tga",
    longestORF = TRUE,
    minimumLength = 0
  )
  ORFs.df <- data.frame(ORFs$`1`)
  # Find location of start codon
  # No ORFs found
  if (nrow(ORFs.df) == 0) {
    print("No ORFs found")
    PTC[,"comments"] <- paste0(PTC[,"comments"],"/No ORFs found")
    PTC[,"PTC_CDS_pos"] <- NA
    return(PTC)
  } else if (!ORFs.df[1,"start"]%in%c(1)) { # The ORF we are looking for starts at first nt position
    print("ORFs do not start at pos 1")
    PTC[,"comments"] <- paste0(PTC[,"comments"],"/ORFs do not start at pos 1")
    PTC[,"PTC_CDS_pos"] <- NA
    return(PTC)
  } else {
    PTC.cds.position <- ORFs.df[1,"end"]-2
    # PTC position should be located after the variant position
    # Not always, if the mutated nt is the second or third position
    if (PTC.cds.position < variant.cds.position) {
      print("PTC position is located before the variant CDS position")
      PTC[,"comments"] <- paste0(PTC[,"comments"],"/PTC position is located before the variant CDS position")
      #return(NA)
    } 
    # Save PTC CDS position
    PTC[,"PTC_CDS_pos"] <- PTC.cds.position
    # Obtain PTC stop codon type
    PTC[,"PTC_stop_codon_type"] <- paste0(fasta.sequence.mutant[PTC.cds.position:(PTC.cds.position+2)],collapse="")
    # PTC position < sum of exons
    # else if (variant.cds.position > sum(exons.size)) {
    #   fs.recurrent.indels.filt[i,"error"] <- "PTC position is higher than the sum of exons"
    #   print(fs.recurrent.indels.filt[i,"error"])
    #   next
    # }
  }
  return(PTC)
}

ATG_finder <- function(PTC_CDS_pos, fasta_sequence) {
        
  res_df <- data.frame(ATG_downstream = NA, ATG_downstream_CDS_pos = NA, ATG_number = NA, 
                          PTC_downstream_4nt = NA, kozak_score_A = NA, kozak_score_B = NA)
  # Split fasta sequence
  fasta_sequence_split <- s2c(fasta_sequence)
  # Obtain fasta sequence from PTC to downstream
  fasta_sequence_downstream_split <- fasta_sequence_split[(PTC_CDS_pos+3):(length(fasta_sequence_split))]
  fasta_sequence_downstream <- paste(fasta_sequence_downstream_split, collapse = "")
  # Obtain 4 nucleotides downstream of PTC
  PTC_downstream_4nt <- paste0(fasta_sequence_downstream_split[1:4],collapse="")
  res_df[,"PTC_downstream_4nt"] <- PTC_downstream_4nt
  # Obtain all ATG positions in all frames
  ATG_positions <- seqinr::words.pos("atg", fasta_sequence_downstream, ignore.case = FALSE, perl = TRUE, fixed = FALSE, useBytes = TRUE)
  if (length(ATG_positions) == 0){
    res_df[,"ATG_downstream"] <- "no"
    return(res_df)
  } else {
    # Check which ones are IN frame ((start position + 2) / 3 should be an integer)
    ATG_positions_inframe <- ATG_positions[which(((ATG_positions+2)/3)%%1 == 0)]
    if (length(ATG_positions_inframe) == 0) {
      res_df[,"ATG_downstream"] <- "no"
      return (res_df)
    } else {
      res_df[,"ATG_downstream"] <- "yes"
      # 1) Obtain the ATG closest to the PTC
      ATG_downstream_CDS_pos <- ATG_positions_inframe[which(ATG_positions_inframe == min(ATG_positions_inframe))]
      # 2) Obtain the total number of ATGs
      ATG_number <- length(ATG_positions_inframe)
      # 3) Check the ATG nucleotide context A and B for the closest one 
      kozak_score_A <- kozak_ATG_context(kozak_pwm = "A", fasta_sequence_downstream_split = fasta_sequence_downstream_split ,ATG_downstream_CDS_pos = ATG_downstream_CDS_pos)
      kozak_score_B <- kozak_ATG_context(kozak_pwm = "B", fasta_sequence_downstream_split = fasta_sequence_downstream_split ,ATG_downstream_CDS_pos = ATG_downstream_CDS_pos)
      # Return output
      res_df[,"ATG_downstream"] <- "yes"
      res_df[,c("ATG_downstream_CDS_pos","ATG_number","kozak_score_A","kozak_score_B")] <- c(ATG_downstream_CDS_pos,ATG_number,kozak_score_A,kozak_score_B)

      return(res_df)
    }
  }
  # This counts the codon usage
  #seqinr::count(fasta_sequence_downstream, wordsize = 3, start = 0, by = 3, freq = FALSE, alphabet = s2c("acgt"))
  # How to specify ambiguous base ? Look for YpR motifs by
  #words.pos("[ct][ag]", myseq) # Should be 1 3
}

kozak_ATG_context <- function(kozak_pwm, fasta_sequence_downstream_split, ATG_downstream_CDS_pos) {
    if (kozak_pwm == "A") {
        nt_context_pre_ATG <- 4
        nt_context_post_ATG <- 1
    } else if (kozak_pwm == "B") {
        nt_context_pre_ATG <- 3
        nt_context_post_ATG <- 2
    }
    if ( ATG_downstream_CDS_pos < nt_context_pre_ATG ) {
        return(NA)
    }
    ATG_sequence_context <- fasta_sequence_downstream_split[(ATG_downstream_CDS_pos-nt_context_pre_ATG):(ATG_downstream_CDS_pos+2+nt_context_post_ATG)]
    ATG_sequence_context <- paste(ATG_sequence_context, collapse="")
    # Create random directory to store outputs 
    # this is to not collapse same files with different samples
    scratch <- sample(1:4)[1]
    folder_name <- stri_rand_strings(n = 1, length = 20, pattern = "[A-Za-z0-9]")
    folder_path <- paste0("/scratch/",scratch,"/gpalou/trash/",folder_name)
    if (file.exists(folder_path)) {
        folder_path <- paste0(folder_path,"_tmp")
    } else {
        dir.create(folder_path, showWarnings = FALSE)
    }
    # Create fasta sequence
    system(paste0("echo -e \">seq1 \n",ATG_sequence_context, "\" > ",folder_path,"/ATG_sequence.txt"))
    # Use MEME fimo to find Kozak sequence context around ATG fasta sequence
    system(paste0("fimo --o ",folder_path,"/fimo_out --thresh 1 /g/strcombio/fsupek_home/gpalou/data/Ignasi/kozak_PWM_",kozak_pwm,".txt ",folder_path,"/ATG_sequence.txt"))
    # Open output
    error <- FALSE
    tryCatch( { 
        kozak_df_res <- read.table(file = paste0(folder_path,"/fimo_out/fimo.tsv"), header = TRUE, sep = "\t")
        kozak_score <- kozak_df_res[kozak_df_res[,"strand"] == "+","score"]
    },error = function(e) {error <<- TRUE})
    # Remove tmp files
    system(paste0("rm -rf ",folder_path))
    if (isTRUE(error)) {
        return(NA)
    } else {
        return(kozak_score)
    }
}

VCF_sample_file <- function(GTEx.sample, type) {
  if (type == "germline") {
    VCF.path.sample <- paste0(VCF.germline.path,paths[paths$folder_or_object=="VCF_germline_file","path_or_filename"])
    VCF.path.sample <- gsub("\\[X1\\]",GTEx.sample,VCF.path.sample)
  } else if (type == "somatic") {
    VCF.path.sample <- paste0(VCF.germline.path,paths[paths$folder_or_object=="VCF_somatic_file","path_or_filename"])
    VCF.path.sample <- gsub("\\[X1\\]",GTEx.sample,VCF.path.sample)
  }
  if (file.exists(VCF.path.sample)) {
    print(VCF.path.sample)
    #VCF.sample <- read.vcfR(VCF.path.sample, verbose = FALSE ) 
    #VCF.sample <- data.frame(VCF.sample@fix)
    VCF.sample <- read.table(file = VCF.path.sample, header = TRUE, sep = "\t")
  } else {
    return(NULL)
  }
  # Fix VAF numeric
  VAFs <- VCF.sample$AF
  VAFs <- ifelse(VAFs == "" | VAFs == ".",0,format(VAFs,scientific = TRUE))
  VAFs <- as.numeric(VAFs)
  VCF.sample$AF <- VAFs
  # Fix ENSEMBL Genes IDs
  # Some transcript swill have different Genes IDs (for whatever reason), but bc it's only used on the regression as gene.id it's not problem
  VCF.sample$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",VCF.sample$Gene.refGene)
  return(VCF.sample)
}

VCF_sample_filtering <- function(VCF.sample, mutation, VAF, filter = NULL) {
    # Filter variants by: nonsense (PTCs) and VAF or synonymous (negative control)
    VCF.sample.filt <- VCF.sample[VCF.sample$ExonicFunc.refGene %in% mutation &  VCF.sample$AF <= VAF,]
    # PASS variants
    if (!is.null(filter)) {
      VCF.sample.filt <- VCF.sample.filt[VCF.sample.filt$Otherinfo10 == filter,]
    }
    return(VCF.sample.filt)
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################
.libPaths( rev( .libPaths() ) )
# Read GTF
library("rtracklayer")
# To read fasta
library("seqinr")
# ORFs
library("ORFik")
# Random strings
library("stringi")
library("stringr")
# somatic VCF strelka VAF calculation
library("ICAMS")

##############################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.sample <- args[2]
VCF.type <- args[3]
GTEx.sample <- gsub("(GTEX\\-\\w{4,5})\\-.*","\\1",GTEx.sample)

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#print(paste0("GTEx tissue --> ",GTEx.tissue))

conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
VCF.germline.path <- paths[paths$folder_or_object=="VCF_germline_file_path","path_or_filename"]
transcripts.fasta.path <- paths[paths$folder_or_object=="transcripts_fasta_path","path_or_filename"]
PTC.transcripts.output.path <- paths[paths$folder_or_object=="PTC_transcripts_output_path","path_or_filename"]

# Output path

output.path <- gsub("\\[X1\\]",GTEx.sample, paste0(PTC.transcripts.output.path,paths[paths$folder_or_object=="germline_PTC_transcripts_output","path_or_filename"]))

if (!file.exists(output.path)) {
  
  # 1.1) ENSEMBL transcripts IDs hg38 GTF
  ensembl.v88.gtf <- rtracklayer::import(paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
  ensembl.v88.gtf <- as.data.frame(ensembl.v88.gtf)
  # Remove versions from IDs
  ensembl.v88.gtf[,c("gene_id","transcript_id")] <- sapply(ensembl.v88.gtf[,c("gene_id","transcript_id")], function(col) { 
    gsub("(.*)\\..*","\\1", col)
  })
  # 1.2) Transcripts fasta sequences
  transcripts.fasta.sequences <- read.fasta(file = paste0(transcripts.fasta.path,paths[paths$folder_or_object=="transcripts_fasta","path_or_filename"]))
  # 2) NMD rules
  print(paste0("Starting NMD rules search for GTEx sample --> ",GTEx.sample))
  print(paste0("PTCs from ",VCF.type," variants"))
  # 2.1) VCF data
  # Obtain VCF germline data file
  # Check if file exists
  VCF.sample <- NULL
  VCF.sample <- VCF_sample_file(GTEx.sample = GTEx.sample, type = VCF.type)
  if(is.null(VCF.sample)){
      print("VCF sample file not found, skipping...")
  } else {
      # 2.2) Filter VCF and check for PTCs (or controls)
      # Filter VCF for variants (stopgain or synonymous)
      mutations <- c("stopgain","frameshift insertion","frameshift deletion")
      VCF.sample.filt <- VCF_sample_filtering(VCF.sample = VCF.sample, 
                                  mutation = mutations, VAF = 1, filter = "PASS")
      if (nrow(VCF.sample.filt) == 0) {
          print("No variants in this sample")
      } else {
          # 2.3) Check NMD rules for each PTC
          variants <- NMD_rules(VCF = VCF.sample.filt, VCF.type = VCF.type)
          transcripts.PTCs <- do.call(rbind,variants)
          write.table(transcripts.PTCs, file = output.path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
      }
  }
} else {
  print("PTCs transcripts file already exist")
}
