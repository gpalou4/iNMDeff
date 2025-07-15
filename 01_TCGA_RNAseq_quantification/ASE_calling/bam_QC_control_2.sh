
#!/usr/bin/bash

STAR_bam_sorted_file_remap=$1
bam_to_remap=$2
bam_to_keep=$3
TCGA_BARCODE=$4
TCGA_PROJECT=$5
THREADS=$6
fastq_var=$7

WASP_DIR="/g/strcombio/fsupek_home/gpalou/software/WASP/mapping"

# Reads Mapping quality
conda activate samtools

samtools view -h -b -q 255 ${STAR_bam_sorted_file_remap} > ${TCGA_BARCODE}_q2.bam
samtools index ${TCGA_BARCODE}_q2.bam

# 3.5) Use filter_remapped_reads.py to create filtered list of reads that correctly remap to same position
conda activate /home/gpalou/anaconda3_envs/general

mkdir filter_remapped_reads

python3  ${WASP_DIR}/filter_remapped_reads.py \
       ${bam_to_remap} \
       ${TCGA_BARCODE}_q2.bam \
       filter_remapped_reads/${TCGA_BARCODE}.remap.keep.bam

# 3.6) # Create a merged BAM containing [1] reads that did not need remapping and [2] filtered remapped reads

mkdir merge_bams

samtools merge merge_bams/${TCGA_BARCODE}.keep.merged.bam \
              ${bam_to_keep} \
              filter_remapped_reads/${TCGA_BARCODE}.remap.keep.bam

# Sort and index the bam file
samtools sort merge_bams/${TCGA_BARCODE}.keep.merged.bam -o merge_bams/${TCGA_BARCODE}.keep.merged.sorted.bam
samtools index merge_bams/${TCGA_BARCODE}.keep.merged.sorted.bam

# 4) Duplicate reads

# With GATK4
# conda activate gatk4
# gatk MarkDuplicates \
#       -I=q1.bam \
#       -O=q2.bam \
#       --REMOVE_DUPLICATES=true \
#       -M=marked_dup_metrics.txt

# WASP in-house script
# Solution to the typical bias from retaining the duplicated read with the highest score, 
# instead of choosing it randomly

if [ ${fastq_var} == "single_end" ]; 
then
        python ${WASP_DIR}/rmdup.py merge_bams/${TCGA_BARCODE}.keep.merged.sorted.bam ${TCGA_BARCODE}.keep.rmdup.merged.sorted.bam
else
        python ${WASP_DIR}/rmdup_pe.py merge_bams/${TCGA_BARCODE}.keep.merged.sorted.bam ${TCGA_BARCODE}.keep.rmdup.merged.sorted.bam
fi
samtools sort ${TCGA_BARCODE}.keep.rmdup.merged.sorted.bam -o ${TCGA_BARCODE}.keep.rmdup.merged.sorted.final.bam
samtools index ${TCGA_BARCODE}.keep.rmdup.merged.sorted.final.bam



