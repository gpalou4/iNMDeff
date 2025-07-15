#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=WASP_test
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=45GB
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

FASTQ_1=$1
FASTQ_2=$2
THREADS=$3
TCGA_BARCODE=$4
out_dir=$5

sleep 30

#/home/gpalou/data/TCGA_RNAseq_quantification/star_index_oh100

# STAR alignment
if [ ${FASTQ_2} == "single_end" ]; 
then
        /bin/bash -c "/src/run_STAR.py \
                /local_scratch/gpalou/tmp/star_index_oh100 \
                ${FASTQ_1} \
                ${TCGA_BARCODE} \
                --threads ${THREADS} \
                --output_dir ${out_dir}"
else
       /bin/bash -c "/src/run_STAR.py \
                /local_scratch/gpalou/tmp/star_index_oh100 \
                ${FASTQ_1} \
                ${FASTQ_2} \
                ${TCGA_BARCODE} \
                --threads ${THREADS} \
                --output_dir ${out_dir}" 
fi