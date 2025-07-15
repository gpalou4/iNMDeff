#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=RSEM_quantification_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

THREADS=$1
TCGA_BARCODE=$2
type=$3

# RSEM transcript quantification
if [ ${type} == "single_end" ]; 
then
        /bin/bash -c "/src/run_RSEM.py \
                /home/gpalou/data/TCGA_RNAseq_quantification/rsem_reference \
                ${TCGA_BARCODE}.Aligned.toTranscriptome.out.bam \
                ${TCGA_BARCODE} \
                --paired_end false \
                -o ./ \
                -t ${THREADS}"
else
        /bin/bash -c "/src/run_RSEM.py \
                /home/gpalou/data/TCGA_RNAseq_quantification/rsem_reference \
                ${TCGA_BARCODE}.Aligned.toTranscriptome.out.bam \
                ${TCGA_BARCODE} \
                -o ./ \
                -t ${THREADS}"
fi