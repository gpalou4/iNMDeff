rm(list=ls())
# To read fasta
library("seqinr")
# Reverse complement sequences
library("spgs")
# Split
library("stringr")
library("dplyr")

# 1.1) ENSEMBL transcripts IDs hg38 GTF
ensembl.v88.gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
#ensembl.v88.gtf <- rtracklayer::import("/g/strcombio/fsupek_cancer1/gpalou/human_genome/Homo_sapiens.GRCh38.107.gtf")
ensembl.v88.gtf <- as.data.frame(ensembl.v88.gtf)
# Remove versions from IDs
ensembl.v88.gtf[,c("gene_id","transcript_id")] <- sapply(ensembl.v88.gtf[,c("gene_id","transcript_id")], function(col) { 
  gsub("(.*)\\..*","\\1", col)
}) 
# 1.2) MANE select + protein coding only
#MANE.transcripts <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/ensembl_MANE_v107.txt")

# 2) Build BED12 file with exon blocks to obtain the fasta sequence of the transcripts

transcripts <- na.omit(unique(ensembl.v88.gtf$transcript_id))

# Using ensembl GTF
i <- 0
bed.file.all <- c()
for (transcript in transcripts) {
  i <- i +1
  print(i)
  # Filter GTF for that transcript
  ensembl.v88.gtf.transcript <- ensembl.v88.gtf[which(ensembl.v88.gtf$transcript_id%in%transcript),]
  ensembl.v88.gtf.transcript.exons <- ensembl.v88.gtf.transcript[ensembl.v88.gtf.transcript$type=="CDS",]
  if (nrow(ensembl.v88.gtf.transcript.exons) == 0) {
      next
  }
  # Order
  ensembl.v88.gtf.transcript.exons <- ensembl.v88.gtf.transcript.exons[order(ensembl.v88.gtf.transcript.exons$start),]
  # Obtain information for the BED file
  # chrom
  chrom <- as.character(ensembl.v88.gtf.transcript[1,"seqnames"])
  # strand
  strand <- as.character(ensembl.v88.gtf.transcript[1,"strand"])
  # start
  start.pos <- ensembl.v88.gtf.transcript.exons$start[1]-1
  # end
  end.pos <- ensembl.v88.gtf.transcript.exons$end[length(ensembl.v88.gtf.transcript.exons$end)]+1
  # score
  score <- "1"
  # thickStart, thickEnd, itemRgb
  thick.start <- thick.end <- item.rgb <- "*"
  # blockCount (exons number), blockStarts (exons start positions) and blockSizes (exons size)
  exons.start <- paste0(ensembl.v88.gtf.transcript.exons$start-ensembl.v88.gtf.transcript.exons$start[1],",")
  exons.number <- length(exons.start)
  exons.start <- paste0(exons.start, collapse = "")
  exons.size <- ensembl.v88.gtf.transcript.exons$end-ensembl.v88.gtf.transcript.exons$start + 1
  exons.size.char <- paste0(exons.size, collapse = ",")
  # Create bedfile
  bed.file <- data.frame(chrom=chrom, start=start.pos, end=end.pos, name=transcript, score=score, strand=strand, thickStart=thick.start, thickEnd=thick.end, itemRgb=item.rgb,
                         blockCount=exons.number, blockSizes=exons.size.char, blockStarts=exons.start)
  bed.file.all <- rbind(bed.file.all,bed.file)
}

# Write bed file
write.table(bed.file.all, file = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts.bed", quote=FALSE, sep='\t',row.names=FALSE, col.names = F)
#write.table(bed.file.all, file = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_V107_coding.bed", quote=FALSE, sep='\t',row.names=FALSE, col.names = F)

# 3) Get fasta sequence with exons concatenated (bash)
# For some bedtools version reason, this only works on node1
system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa -name -bed /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts.bed > /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts.fasta"))
#system(paste0("bedtools getfasta -split -fi /g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa -name -bed /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_V107_coding.bed > /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_V107_coding.fasta"))

# 4) Build BED12 file with intron blocks to obtain the fasta sequence of the transcripts

# Using ensembl GTF
i <- 0
bed.file.all <- c()
for (transcript in transcripts) {
  i <- i + 1
  print(i)
  # Filter GTF for that transcript
  ensembl.v88.gtf.transcript <- ensembl.v88.gtf[which(ensembl.v88.gtf$transcript_id%in%transcript),]
  ensembl.v88.gtf.transcript.exons <- ensembl.v88.gtf.transcript[ensembl.v88.gtf.transcript$type=="CDS",]
  if (nrow(ensembl.v88.gtf.transcript.exons) <= 1) {
      next
  }
  # Order
  ensembl.v88.gtf.transcript.exons <- ensembl.v88.gtf.transcript.exons[order(ensembl.v88.gtf.transcript.exons$start),]
  # Obtain information for the BED file
  # chrom (there could be one transcript ID in 2 chromosomes, X and Y)
  chrs <-  as.character(unique(ensembl.v88.gtf.transcript[,"seqnames"]))
  for (chrom in chrs) {
    ensembl.v88.gtf.transcript.exons.filt <- ensembl.v88.gtf.transcript.exons[which(ensembl.v88.gtf.transcript.exons$seqnames %in% chrom),] 
    if (nrow(ensembl.v88.gtf.transcript.exons.filt) <= 1) {
        next
    }
    # strand
    strand <- as.character(ensembl.v88.gtf.transcript.exons.filt[1,"strand"])
    # start
    start.pos <- ensembl.v88.gtf.transcript.exons.filt$end[1]
    # end
    end.pos <- (ensembl.v88.gtf.transcript.exons.filt$end[length(ensembl.v88.gtf.transcript.exons.filt$end)-1])
    # score
    score <- "1"
    # thickStart, thickEnd, itemRgb
    thick.start <- thick.end <- item.rgb <- "*"
    # blockCount (introns number), blockStarts (introns start positions) and blockSizes (introns size)
    introns.start <- paste0( (ensembl.v88.gtf.transcript.exons.filt$end[-length(ensembl.v88.gtf.transcript.exons.filt$end)]) - (ensembl.v88.gtf.transcript.exons.filt$end[1]),",")
    introns.number <- length(introns.start)
    introns.start <- paste0(introns.start, collapse = "")
    introns.size <- (ensembl.v88.gtf.transcript.exons.filt$start[-1] -1) - (ensembl.v88.gtf.transcript.exons.filt$end[-length(ensembl.v88.gtf.transcript.exons.filt$end)])
    introns.size.char <- paste0(introns.size, collapse = ",")
    # Create bedfile
    bed.file <- data.frame(chrom=chrom, start=start.pos, end=end.pos, name=transcript, score=score, strand=strand, thickStart=thick.start, thickEnd=thick.end, itemRgb=item.rgb,
                          blockCount=introns.number, blockSizes=introns.size.char, blockStarts=introns.start)
    bed.file.all <- rbind(bed.file.all,bed.file)
  }
}

# Write bed file
write.table(bed.file.all, file = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_introns.bed", quote=FALSE, sep='\t',row.names=FALSE, col.names = F)

# Transcripts with only 1 intron (== 2 exons) have same start and same end, only 1 feature block --> fasta cannot be retrieved! What do we do?

# 5) Get fasta sequence with introns concatenated (bash)
# For some bedtools version reason, this only works on node1
system(paste0("bedtools getfasta -split -s -fi /g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa -name -bed /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_introns.bed > /g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_introns.fasta"))



bed.file.all[grep("2103365",bed.file.all$start),]

# 6) Obtain intronic sequences separated
transcripts_fasta_sequences <- read.fasta(file = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/bed12_transcripts_introns.fasta")

transcript_and_intron_sequences <- lapply(transcripts_fasta_sequences, function(transcript) {

  ensembl_transcript_id <- gsub("(ENST.*)\\:\\:.*","\\1",attr(transcript,"name"))
  #index <- grep(ensembl_transcript_id,names(transcripts_fasta_sequences))
  #transcript_fasta_sequence <- transcripts_fasta_sequences[[index]]
  fasta_sequence <- tolower(unlist(strsplit(transcript,"")))
  # Coordinates of introns
  bed_file_transcript <- bed.file.all %>% 
                        filter(name %in% ensembl_transcript_id)
  introns_size <- as.numeric(unlist(strsplit(as.character(bed_file_transcript$blockSizes),",")))
  strand <- as.character(bed_file_transcript$strand)
  if (strand == "-") {
    introns_size <- rev(introns_size)
  }

  all_intron_sequences <- c()
  for (i in 1:length(introns_size)) {
    intron_size <- introns_size[i]
    # Obtain intronic sequence
    intron_sequence <- fasta_sequence[(1:intron_size)]
    # Modify fasta sequence
    fasta_sequence <- fasta_sequence[-(1:intron_size)]
    # Save
    intron_sequence <- paste0(intron_sequence,collapse = "")
    # 
    if (length(all_intron_sequences) == 0){
      all_intron_sequences <- intron_sequence
    } else {
      all_intron_sequences <- c(all_intron_sequences,intron_sequence)
    }
  }
  all_intron_sequences_collapsed <- paste0(all_intron_sequences, collapse = ",")
  data.frame(ensembl_transcript_id = ensembl_transcript_id, intron_sequences = all_intron_sequences_collapsed)

})

# Save

df <- do.call(rbind, transcript_and_intron_sequences)
write.table(df, file = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/transcripts_and_intron_sequences.txt", 
            quote = FALSE, sep = '\t',row.names = FALSE, col.names = TRUE)




