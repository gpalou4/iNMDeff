#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMDrules_array
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
##SBATCH --exclude=fsupeksvr5
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

#conda activate /home/dnaro/.conda/envs/NMD_regression_3

commands_file=$1

bash /home/gpalou/projects/NMD/scripts/slurm_arrays/execute_lines.sh $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_STEP ${commands_file}


