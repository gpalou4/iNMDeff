#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=STAR_align_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=40GB
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

FASTQ_1=$1
FASTQ_2=$2
THREADS=$3
TCGA_BARCODE=$4

# STAR alignment
if [ ${FASTQ_2} == "single_end" ]; 
then
        /bin/bash -c "/src/run_STAR.py \
                /local_scratch/gpalou/tmp/star_index_oh100 \
                ${FASTQ_1} \
                ${TCGA_BARCODE} \
                --threads ${THREADS} \
                --output_dir ./"
else
       /bin/bash -c "/src/run_STAR.py \
                /local_scratch/gpalou/tmp/star_index_oh100 \
                ${FASTQ_1} \
                ${FASTQ_2} \
                ${TCGA_BARCODE} \
                --threads ${THREADS} \
                --output_dir ./" 
fi

# singularity run /g/strcombio/fsupek_home/gpalou/singularity/gtex_rnaseq_V8.sif \
#     /bin/bash -c "/src/run_STAR.py \
#         /home/gpalou/data/TCGA_RNAseq_quantification/star_index_oh100 \
#         ${FASTQ_1_PATH_FIX} \
#         ${FASTQ_2_PATH_FIX} \
#         ${TCGA_BARCODE} \
#         --threads ${THREADS} \
#         --output_dir /home/gpalou/data/TCGA_RNAseq_quantification/star_out_tmp"

# singularity run /g/strcombio/fsupek_home/gpalou/singularity/gtex_rnaseq_V8.sif \
#     /bin/bash -c "/src/run_STAR.py \
#         /home/gpalou/data/TCGA_RNAseq_quantification/star_index_oh100 \
#         /home/gpalou/data/Marina/SRR14301954.fastq \
#         SRR14301954 \
#         --threads 6 \
#         --output_dir /home/gpalou/data/Marina"