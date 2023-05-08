# S. cerevisiae response to lactic acid enanthiomers

This repository contains the data and code for the manuscript on transcriptional responses of *S. cerevisiae* to D-lactic acid and L-lactic acid (to be submitted).

Steps and tools: 
 1. QC: FastQC v0.11.9 and MultiQC v1.13
 2. Reference: R64 release of S. cerevisiae strain S288C genome
 3. Alignment: hisat2 v2.2.1
 4. Sorted with samtools v1.9 and quantified with the featureCounts tool from subread v2.0.4.
 5. Differentially expressed genes identified with DESeq2 v1.34.0 for R v4.1.2.



Other useful pieces of information (not exactly used in the final version):
  1. https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
  2. on linear models: https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_2021_June/Markdowns/09_Linear_Models.html
  3. edgeR lab: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
  4. comparison of DE methods: https://www.jove.com/t/62528/three-differential-expression-analysis-methods-for-rna-sequencing
  
