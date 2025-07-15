#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMDrules_command
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
#SBATCH --mem=25GB
#SBATCH --time=15:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

#conda activate /home/dnaro/.conda/envs/NMD_regression_3
PATH_FILE=$1
GTEx_tissue=$2
GTEx_sample=$3
VCF_type=$4
DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx"

GTEx_sample_filt=$(echo $GTEx_sample | cut -d '-' -f1-2)

if [ ${VCF_type} == "germline" ]; then
    GTEx_sample_output_file_1="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/${GTEx_sample_filt}_germline_PTC_transcripts.txt"
    GTEx_sample_output_file_2="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_sample_filt}_AAAgermline_PTC_transcripts_metadata.txt"
else
    GTEx_sample_output_file_1="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/${GTEx_sample_filt}_somatic_PTC_transcripts.txt"
    GTEx_sample_output_file_2="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_sample_filt}_somatic_PTC_transcripts_metadata.txt"
fi

if [ -f ${GTEx_sample_output_file_1} ]; then
        echo "${GTEx_sample} sample already done for PTC transcripts file, skipping -->" ${GTEx_sample_output_file_1}
    if [ -f ${GTEx_sample_output_file_2} ]; then
        echo "${GTEx_sample} sample already done for PTC transcripts metadata file, skipping -->" ${GTEx_sample_output_file_2}
    else
        # PTC metadata
        bash ${DIR}/PTCs_NMDeff_metadata.sh ${PATH_FILE} ${GTEx_tissue} ${GTEx_sample} ${VCF_type}
    fi
else
    # Check NMD rules for all PTCs from the sample
    bash ${DIR}/NMD_rules.sh ${PATH_FILE} ${GTEx_sample} ${VCF_type}
    bash ${DIR}/PTCs_NMDeff_metadata.sh ${PATH_FILE} ${GTEx_tissue} ${GTEx_sample} ${VCF_type}
fi

