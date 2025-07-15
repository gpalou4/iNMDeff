#!/usr/bin/bash

WD="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/GTEx"

> ${WD}/sbatch_commands.txt

# With gene list
# while read tissue
# do
#     echo $tissue
#     for randomization in {no,yes}
#     do
#         for randomGenes in {no,yes}
#         do
#             for NMD_method in {ASE,endogenous}
#             do
#                 for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc}
#                 do
#                     echo "${WD}/3_GTEx_rare_germline_variants_associations.sh $tissue $dataset $NMD_method 0.2 yes 2 $randomGenes $randomization" >> /home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/GTEx/sbatch_commands.txt
#                 done
#             done
#         done
#     done
# done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt

# Without gene list

while read tissue
do
    echo $tissue
    for randomization in {no,yes}
    do
        for NMD_method in {ASE,endogenous}
        do
            for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc}
            do
                echo "${WD}/3_GTEx_rare_germline_variants_associations.sh $tissue $dataset $NMD_method 0.2 no 2 no $randomization" >> /home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/GTEx/sbatch_commands.txt
            done
        done
    done
done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt