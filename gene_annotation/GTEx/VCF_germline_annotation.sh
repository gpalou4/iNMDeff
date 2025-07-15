#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_annovar
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

while read VCF_germline_sample
do
        VCF_germline_annotation_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs"

	# if [ -f "${VCF_germline_annotation_path}/${VCF_germline_sample}_liftOver_annotated.HG_multianno.txt.gz" ]; then
        #         echo "Already done, skipping --> ${VCF_germline_annotation_path}/${VCF_germline_sample}_liftOver_annotated.HG_multianno.txt.gz"
        # else
		sbatch VCF_germline_annotation_sample.sh $VCF_germline_sample $VCF_germline_annotation_path
                sleep 2
        # fi
done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt

