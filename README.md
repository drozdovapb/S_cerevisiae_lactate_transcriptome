# S. cerevisiae response to lactic acid enanthiomers

This repository contains the data and code for the manuscript on transcriptional responses of *S. cerevisiae* to D-lactic acid and L-lactic acid (to be submitted).

Steps and tools: 
 1. QC: FastQC v0.11.9 and MultiQC v1.13
 2. Reference: R64 release of S. cerevisiae strain S288C genome
 3. Alignment: hisat2 v2.2.1
 4. Sorted with samtools v1.9 and quantified with the featureCounts tool from subread v2.0.4.
 5. Differentially expressed genes identified with DESeq2 v1.34.0 for R v4.1.2.
