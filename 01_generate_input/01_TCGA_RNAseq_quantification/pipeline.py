import pandas as pd
import subprocess
import os
import sys

gdc_token = "gdc-user-token.2022-03-04T09_27_50.252Z.txt
nf_path = "/home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification/nextflow"
batches_path = "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/sample_batches"
wd = os.getcwd()

########### DO THIS ONLY ONCE TO GENERATE THE INPUT BATCH FILES #############
# TCGA_metadata_path = "/g/strcombio/fsupek_home/gpalou/data/TCGA_fastq_RNAseq/TCGA_fastq_RNAseq_metadata.txt"
# TCGA_metadata = pd.read_csv(filepath_or_buffer = TCGA_metadata_path, sep = "\t")
# # Create sample_batches_files
# # Each batch containing 5 samples
# # Each nextflow will run for each batch
# # 9195/5=1839 batches
# batch_n = 0
# for i in range(0,int(TCGA_metadata.shape[0]/5)):
#      batch_n = batch_n + 5
#      TCGA_metadata_filt = TCGA_metadata.iloc[batch_n:batch_n+5,:]
#      TCGA_metadata_filt = TCGA_metadata_filt.astype(str)
#      batch_path = "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/sample_batches_all/" + "batch_" + str(i) + ".csv"
#       #Save the file
#      TCGA_metadata_filt.to_csv(batch_path, index = False)
###########################################
#TCGA_metadata.iloc[[232,522,4532,4923,6534,8121,6123,6321],].to_csv("/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/sample_batches_ASE/prueba.csv", index = False)


batch_numbers = sys.argv[1]
batch_numbers = list(batch_numbers.split(","))
node = sys.argv[2]

# 1) Go to /scratch working directory from the specified node
cmd = ['cd','/scratch/'+str(node)+'/gpalou/']
subprocess.run(cmd, 
     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
     check=True)

# Create tmp_dir for the WorkDir nextflow folder 
# Each tmp_dir must be different each time we execute the script to not overlap rm -rf WorkDir between samples)
tmp_dir = "tmp_"+"_".join(batch_numbers)
workDirPath = "/scratch/"+str(node)+"/gpalou/tmp/"+tmp_dir
os.makedirs(workDirPath, exist_ok = True)

for i in range(int(batch_numbers[0]),int(batch_numbers[1])+1):

    # 2) Obtain Sample Batch metadata file path
    print("BATCH number--> "+str(i))
    TCGA_metadata_batch_path = batches_path+'/batch_'+str(i)+'.csv'
    if os.path.isfile(TCGA_metadata_batch_path):
        TCGA_metadata_batch = pd.read_csv(filepath_or_buffer = TCGA_metadata_batch_path, sep = "\t")    
        print(TCGA_metadata_batch)
    else:
        print("Metadata file is not found, probably the sample batch is already done")
        continue

    # 3) Run Nextflow
    cmd = ['nextflow','run',nf_path+'/main.nf','-c',nf_path+'/nextflow_slurm.config', '--input_csv', TCGA_metadata_batch_path, 
            '--scratchDir', '-resume','--cpus', '2','--gdc_token',gdc_token,'--gdc_download_tries', '10']
    
    print("Nextflow command --> \n"+" ".join(cmd))

    # Create work directory for the batch
    batch_workDirPath = workDirPath+'/batch_'+str(i)
    os.makedirs(batch_workDirPath, exist_ok = True)

    results = subprocess.run(cmd, cwd = batch_workDirPath,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)

    # 4) If successful, just move the metadata file to the "samples done" directory
    if results.returncode == 0:
        subprocess.run(['mv',TCGA_metadata_batch_path, batches_path+'_downloaded_done'], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)

    # 5) Remove nextflow working directory for the specific batch
    subprocess.run(['rm','-rf',batch_workDirPath], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)

# 6) Remove tmp_dir for all the batches
subprocess.run(['rm','-rf',workDirPath], 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    check=True)
    

