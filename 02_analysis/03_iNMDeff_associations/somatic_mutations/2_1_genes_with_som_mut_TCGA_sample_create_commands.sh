#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=CGC_somatic_mut_CNV_cancer
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=1GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

# TCGA cancer tissue
TCGA_cancer=$1
DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"
OUT_DIR="/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/${TCGA_cancer:5}"
# Obtain samples
TCGA_samples=$(ls -lh /g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/${TCGA_cancer:5}/ | grep -v ${TCGA_cancer:5} | awk '{print $9}')
# Create commands_file.txt
> "${OUT_DIR}/${TCGA_cancer}_genes_som_mut_array_commands.txt"

for TCGA_sample in ${TCGA_samples}
    do
    # Check if sample is already done
    TCGA_sample_output_file="/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/${TCGA_cancer:5}/${TCGA_sample}/${TCGA_sample}_CGC_somatic_mut_CNV.txt"
    # Check somatic mut and CNV for that sample
	#if [ -f ${TCGA_sample_output_file} ]; then
    #    echo "${TCGA_sample} sample already done, skipping -->" ${TCGA_sample_output_file}
    #else
    echo ${DIR}/2_genes_with_som_mut_TCGA_sample.sh ${TCGA_cancer} ${TCGA_sample} >> "${OUT_DIR}/${TCGA_cancer}_genes_som_mut_array_commands.txt"
    #fi
    done
