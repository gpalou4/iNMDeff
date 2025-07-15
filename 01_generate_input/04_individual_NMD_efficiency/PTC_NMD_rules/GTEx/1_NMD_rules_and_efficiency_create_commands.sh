#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=NMDeff_rules_create_commands
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

GTEx_tissue=$1
VCF_type=$2
DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx"

# Obtain samples
GTEx_samples=$(head -1 /g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/tissues/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm_${GTEx_tissue}.gct | cut -f3-)
#GTEx_samples=$(cat "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt")

>"/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_tissue}_${VCF_type}_PTC_NMD_rules_samples_commands.sh"

for GTEx_sample in ${GTEx_samples}
    do
    # Check NMD rules for all PTCs from the sample
    echo ${DIR}/NMD_rules_and_efficiency_command.sh NMD_rules_and_efficiency.txt ${GTEx_tissue} ${GTEx_sample} ${VCF_type} >> /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_tissue}_${VCF_type}_PTC_NMD_rules_samples_commands.sh
    done
