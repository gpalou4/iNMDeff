#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=transcripts_filt_somaticMut_CNV
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/data/NMD_project/transcripts_to_remove/ENSEMBL/out/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# Activate base environment
conda activate /home/gpalou/anaconda3_envs/general/

TCGA_project=$1

python3 transcripts_filt_somaticMut_CNV.py ${TCGA_project}
