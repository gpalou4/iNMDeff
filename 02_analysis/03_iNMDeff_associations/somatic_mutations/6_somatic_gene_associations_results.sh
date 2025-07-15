#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=normal_prio
#SBATCH --job-name=CGC_analysis
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --exclude=fsupeksvr5
#SBATCH --mem=6GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment

conda activate /home/gpalou/anaconda3_envs/general

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"
pancancer=$1 # yes//no
NMD_method=$2 # ASE//endogenous
test=$3 # CGC//subtypes
PCs_optimized=$4 # yes//no	 

# Check somatic mut and CNV for that sample
Rscript ${DIR}/6_somatic_gene_associations_results.R $pancancer $NMD_method $test $PCs_optimized



