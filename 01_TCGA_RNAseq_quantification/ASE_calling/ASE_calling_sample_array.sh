#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=ASE_array
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
#SBATCH --mem=25GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
##SBATCH --exclude=fsupeksvr5
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

#conda activate /home/dnaro/.conda/envs/NMD_regression_3

commands_file=$1

l=`head -n${SLURM_ARRAY_TASK_ID} ${commands_file} | tail -n1`
eval $l



