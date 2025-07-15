#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=CGC_somatic_mut_CNV_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=4GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

TCGA_cancer=$1
TCGA_sample=$2
DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"

# Check somatic mut and CNV for that sample
Rscript ${DIR}/2_genes_with_som_mut_TCGA_sample.R ${DIR}/2_genes_with_som_mut_TCGA_sample.txt ${TCGA_cancer} ${TCGA_sample}

