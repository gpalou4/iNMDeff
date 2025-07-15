

echo $cancer
ll /home/gpalou/analysis_results/NMD_project/associations/TSG_OG_associations/lambda/${cancer}_ASE_CGC_somatic_mut_CNV_PCs_*_lambda.txt | grep "oct" | awk '{print $9}' > ${cancer}_tmp.txt

for file in $(cat ${cancer}_tmp.txt);
do
    newfile=$(echo $file | sed 's/ASE_CGC/ASE_0.01_CGC/g')
    mv $file $newfile
done



