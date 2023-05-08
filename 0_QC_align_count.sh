## quality control
fastqc * 
multiqc .

## get reference
cd ~/SC_DL/1_refs
wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz 
gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz 

## navigate to the alignment directory and set paht
cd ~/SC_DL/2_aln_2
export dir="../0_raw_data/S6614_1/01_fastq/"

## hisat: build index and prepare splice file
hisat2-build ../1_refs/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index
hisat2_extract_splice_sites.py ../1_refs/Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt

## hisat: align & samtools: make sorted bam
for sample in `ls ../0_raw_data/S6614_1/01_fastq/RNA_S6614Nr*1.fastq.gz`; \
  do base=$(basename $sample ".1.fastq.gz"); hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 8 \
  -1 $dir/$base.1.fastq.gz -2 $dir/$base.2.fastq.gz | samtools view --threads 6 -bS | samtools sort --threads 6 -o $base.bam; \
  done

## featureCounts: quantify
/media/secondary/apps/subread-2.0.4-Linux-x86_64/bin/featureCounts -a ../1_refs/Saccharomyces_cerevisiae.R64-1-1.108.gtf \
  -s 2 -T 6 -p -o allSamples.featureCounts.txt $(ls *.bam)
