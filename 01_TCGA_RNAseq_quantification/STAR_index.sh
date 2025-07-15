#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=STAR_index_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=60GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# Activate base environment

# build the STAR index:
## sjdbOverhang is 75 in GTEx but 100 in TCGA

# singularity run /g/strcombio/fsupek_home/gpalou/singularity/gtex_rnaseq_V8.sif \
#     /bin/bash -c "STAR \
#         --runMode genomeGenerate \
#         --genomeDir /home/gpalou/data/TCGA_RNAseq_quantification/star_index_oh100 \
#         --genomeFastaFiles /home/gpalou/data/human_genome/TCGA/GRCh38.d1.vd1.fa \
#         --sjdbGTFfile /home/gpalou/data/NMD_project/conversor_tables/gencode.v26.annotation.gtf \
#         --sjdbOverhang 100 \
#         --runThreadN 8"

scratch=$1

singularity run -B /scratch/${scratch}/gpalou/ /g/strcombio/fsupek_home/gpalou/singularity/gtex_rnaseq_V8.sif \
    /bin/bash -c "STAR \
        --runMode genomeGenerate \
        --genomeDir /scratch/${scratch}/gpalou/tmp/star_index_oh100 \
        --genomeFastaFiles /scratch/${scratch}/gpalou/tmp/TCGA/GRCh38.d1.vd1.fa \
        --sjdbGTFfile  /scratch/${scratch}/gpalou/tmp/gencode_annotation/gencode.v26.annotation.gtf \
        --sjdbOverhang 100 \
        --runThreadN 8"
