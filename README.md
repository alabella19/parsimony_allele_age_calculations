# parsimony_allele_age_calculations
This repository contains the scripts and instructions for calculating parsimony based age calculations of SNP alleles from mammalian alignments

---
## Step 1 - obtaining hg18, hg19 and hg38 positions for all loci based on RS number

# Requirements
R 

# Input 
text file with 1 rs value per line
example: input_rs.txt

# Output
input_rs.chrN.txt - hg38 info for chromosome N
input_rs.rs_hg18_chrN.txt - hg18 info for chromosomes N*
input_rs.rs_hg19_chrN.txt - hg19 info for chromosome N*
input_rs.rs_hg19.txt - hg19 info for all chromosomes

NOTE: Biomart returns hits from other chromosomes when retreiving hg18 and hg19. These should be deleted before moving on.

# Execution
Rscript biomart.R input_rs.txt

