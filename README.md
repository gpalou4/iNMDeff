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
    - [02_VCF_variant_annotation](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/02_VCF_variant_annotation) --> Variant annotation using _ANNOVAR_, for GTex and TCGA germline _VCF_ files.
    - [03_NMD_gene_sets](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/03_NMD_gene_sets) -->
        1) Checks for NMD-triggering/evading features across the v88 ENSEMBL _gtf_ annotation file: uORFs at 5'UTR and GC content at 3'UTR.
        2) Obtains NMD gene sets, from different sourced articles
        3) Updates gene symbol to the same current version
        4) For each NMD geneset from 2), generates different new features for the ENSEMBL transcripts (#ORFs and length, transcript length, RNA-seq TPM expression, MANE transcript) and chooses a tentative pair of transcripts (NMD-target and control) for each gene

2. 


