#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=array_endogenous_TCGA_NMD_efficiency
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
##SBATCH --exclude=fsupeksvr5
#SBATCH --mem=1GB
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA/"
commands_file=$1

bash /home/gpalou/projects/NMD/scripts/slurm_arrays/execute_lines.sh $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_STEP ${commands_file}


