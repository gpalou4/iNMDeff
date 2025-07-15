#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=ASE_cancer
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --exclude=fsupeksvr5
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

TCGA_CANCER=$1
type=$2 # germline // somatic

conda activate /home/dnaro/.conda/envs/NMD_regression_3

# Cancer samples
TCGA_cancer_path="/g/strcombio/fsupek_cancer1/gpalou/TCGA_${type}_variants/${TCGA_CANCER:5}"
TCGA_cancer_samples=$(ls ${TCGA_cancer_path} | grep -v ${TCGA_CANCER})
DIR="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_CANCER}"

cat > "${DIR}/${TCGA_CANCER}_ASE_samples_commands.sh"

for TCGA_SAMPLE in ${TCGA_cancer_samples}
do
    ASE_calls_file_path="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_CANCER}/${TCGA_SAMPLE}/ASE_${type}/AASE_${type}_strelka_calls.txt"
    if [ -f ${ASE_calls_file_path} ]; then
        echo "${TCGA_SAMPLE} sample already done, skipping -->" ${ASE_calls_file_path}
    else
        echo "/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification/ASE_calling/ASE_calling_sample.sh ${TCGA_CANCER} ${TCGA_SAMPLE} ${type}" >> "${DIR}/${TCGA_CANCER}_ASE_samples_commands.sh"
    fi
done



