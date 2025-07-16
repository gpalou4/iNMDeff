# Individual NMD efficiency (iNMDeff)

Repository to reproduce analysis, figures and tables of the research article titled "_Variable efficiency of nonsense-mediated mRNA decay across human tissues, tumors and individuals_"

<p align="center">
  <img
    src="./Fig1.png"       
    alt="Variability and determinants of NMD efficiency across human tissues, tumors, and individuals"
    width="630">
</p>

# Repository structure

1. [01_generate_input](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input) --> Generates all input data required for the **individual NMD efficiency (iNMDeff)** estimation and downstream analysis.

    - [01_TCGA_RNAseq_quantification](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/01_TCGA_RNAseq_quantification) --> Nextflow pipelines for 1) RNA-seq gene/transcript-level quantification and 2) allele-specific expression (ASE) variant calling. It includes _fastq_ alignment to human genome with _STAR_, transcriptomic quantification using _RSEM_ or allele-specific expression (ASE) variant calling for germline variants using _Strelka2_ across all ~10k samples from the TCGA cohort.
    - [02_VCF_variant_annotation](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/02_VCF_variant_annotation) --> Variant annotation using _ANNOVAR_, for GTex and TCGA germline _VCF_ files.
    - [03_NMD_gene_sets](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/03_NMD_gene_sets) -->
        1) Checks for NMD-triggering/evading features across the v88 ENSEMBL _gtf_ annotation file: uORFs at 5'UTR and GC content at 3'UTR.
        2) Obtains NMD gene sets, from different sourced articles.
        3) Updates gene symbols to the same current version.
        4) For each NMD gene set generated above, creates different new features for the ENSEMBL transcripts (#ORFs and length, transcript length, RNA-seq TPM expression, MANE transcript) and classifies transcripts as NMD-target or controls (non NMD-target) for each gene.
        5) Chooses the best NMD-target and control transcripts for every gene across all our NMD gene sets.
    - [04_individual_NMD_efficiency](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/04_individual_NMD_efficiency)
        1) [PTC_NMD_rules](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/04_individual_NMD_efficiency/PTC_NMD_rules). Classifies premature termination codons (PTCs, from nonsense and indel mutations) as NMD-triggering or NMD-evading based on canonical known genomic NMD rules. Done for GTex and TCGA individuals separately. For indels, the position of the downstream generated PTC is firstly predicted.
        2) Estimates **individual NMD efficiency (iNMDeff)** separately for [**ASE**](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/04_individual_NMD_efficiency/ASE) and [**endogenous target gene (ETG)**](https://github.com/gpalou4/iNMDeff/tree/main/01_generate_input/04_individual_NMD_efficiency/endogenous_target_gene) methods (and separately for every NMD gene set).
           
2. [02_analysis](https://github.com/gpalou4/iNMDeff/tree/main/02_analysis) --> Main analysis are iNMDeff variability and associations with genetic mutations, CNAs, survival and immunotherapy response.
    - [01_iNMDeff_robustness](https://github.com/gpalou4/iNMDeff/tree/main/02_analysis/01_iNMDeff_robustness) --> Assess robustness of ETG and ASE iNMDeff methods by different orthogonal approaches: cell lines with UPF1 KD, PRO-seq data, ASE-ETG correlations, comparison with tissue rankings from Teran et. al. (GTex only), etc.
    - [02_iNMDeff_variability](https://github.com/gpalou4/iNMDeff/tree/main/02_analysis/02_iNMDeff_variability) --> iNMDeff variability analysis
       1) Randomization tests for iNMDeff variability: Inter-Tissue iNMDeff variability deviation (ITNVD) test and Tissue iNMDeff Deviation (TND)
       2) Plot for ITNVD test
       3) General linear model (glm) to predict iNMDeff and assess variability explained by the different explanatory variables such as CNA, cancer or tissue type, sex, tumor purity, etc.
       4) Inter-individual iNMDeff variability tests
    - [03_iNMDeff_associations](https://github.com/gpalou4/iNMDeff/tree/main/02_analysis/03_iNMDeff_associations) --> iNMDeff associations
       1) Cell type deconvolution in brain (proportion of glia vs neuron cells)
       2) Rare variant association analysis (RVAS) of putative loss-of-function (pLoF) germline variants across all protein coding genes. This was done using SKAT-O, separately for ETG and ASE iNMDeff, for 18 matched GTex tissues and TCGA cancer types, and 3 pLoF datasets: i) One threshold involved only NMD-triggering PTC variants, and the other two additionally included the predicted deleterious missense variants using CADD scores at two different cutoffs: >= 25 (more stringent) and >= 15 (permissive).
       3) Immune infiltration analysis from TCGA and Hartwig cohorts for CD8+ cytotoxic T-cells
       4) Immunotherapy response associations for external cohorts from different cancer types. Includes survival analysis for OS and PFS, and immunotherapy response prediction by logistic regression.
       5) Somatic CNA associations. Firstly, [sparse-PCA](https://github.com/gpalou4/iNMDeff/tree/main/02_analysis/03_iNMDeff_associations/somatic_CNAs/01_sparse_PCA_CNA) is used to obtain CNA principal component signatures (CNA-PCs)


