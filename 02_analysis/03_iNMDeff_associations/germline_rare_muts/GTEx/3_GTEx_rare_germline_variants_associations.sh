#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=GTEx_RGWAS
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=150GB
##SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate RGWAS

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/GTEx"

tissue=$1 # pantissue // [X]
rare_germline_variants_name=$2
#PTV_Missense_CADD15_0.1perc
#PTV_Missense_CADD25_0.1perc
#PTV_0.1perc
#MRT_25thCentile_0.1perc
#CCR_90orhigher_0.1perc
NMD_method=$3 # ASE // endogenous
VAF=$4 # 0.2 // 0.01
genelist=$5 # yes // no
numberGenesThres=$6
randomGenes=$7 # yes // no
randomization=$8 # yes // no

# Check somatic mut and CNV for that sample
Rscript ${DIR}/3_GTEx_rare_germline_variants_associations.R $tissue $rare_germline_variants_name $NMD_method $VAF $genelist $numberGenesThres $randomGenes $randomization

