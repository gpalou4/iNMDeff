
for alpha in {1e-04,0.00015,2e-04,3e-04,5e-04,7e-04} 
do
    for PC in {50,150,200}
    do 
        sbatch 1_TCGA_PCA_CNV.sh pancancer $alpha no $PC
    done
done


