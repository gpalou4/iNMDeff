# Individual NMD efficiency (iNMDeff)

Repository to reproduce analysis, figures and tables of the research article titled "_Variable efficiency of nonsense-mediated mRNA decay across human tissues, tumors and individuals_"

DOI: X



<p align="center">
  <img
    src="./Fig1.png"       
    alt="Variability and determinants of NMD efficiency across human tissues, tumors, and individuals"
    width="630">
</p>



# Repository structur

1. [01_generate_input](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input) --> Contains code to generate all input data required for the analysis.

    - [01_TCGA_RNAseq_quantification](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/01_TCGA_RNAseq_quantification) --> Nextflow pipelines for 1) RNA-seq gene/transcript-level quantification and 2) allele-specific expression (ASE) variant calling. It includes _fastq_ alignment to human genome with _STAR_, transcriptomic quantification using _RSEM_ or allele-specific expression (ASE) variant calling for germline variants using _Strelka2_ across all ~10k samples from the TCGA cohort.

2. 


