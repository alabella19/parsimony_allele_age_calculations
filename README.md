# parsimony_allele_age_calculations
This repository contains the scripts and instructions for calculating parsimony based age calculations of SNP alleles from mammalian alignments

---
# Step 1 - obtaining hg18, hg19 and hg38 positions for all loci based on RS number

### Requirements
R 

R - biomaRt package (https://bioconductor.org/packages/release/bioc/html/biomaRt.html) 

R - readr

### Input 
text file with 1 rs value per line

example: input_rs.txt

### Output
input_rs.chrN.txt - hg38 info for chromosome N

input_rs.rs_hg18_chrN.txt - hg18 info for chromosomes N* (no rs values)

input_rs.rs_hg18_chrN.ref - hg18 info for chromosomes N* (with rs values)

input_rs.rs_hg19_chrN.txt - hg19 info for chromosome N* (no rs values)

input_rs.rs_hg19_chrN.ref - hg19 info for chromosome N* (with rs values)

input_rs.rs_hg19.txt - hg19 info for all chromosomes

NOTE: Biomart returns hits from other chromosomes when retreiving hg18 and hg19. These should be deleted before moving on.

### Execution
`Rscript biomart.R input_rs.txt`

# Step 2 - obtain multiz data from hg38. Done for each chromosome

### Requirements
PHAST package (https://github.com/CshlSiepelLab/phast) 

### Input
MAF alignment obtained from UCSC human genome database

hg38 information - example: input_rs.chr21.txt

### Output
Multiz alignment for each chromosome at each site - input_rs.chrN.multiz

NOTE: The multiz alignment can be mapped to the SNP based on the position

### Execution (done for each chromosome)
`wget -nd -q "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/maf/chr21.maf.gz"

gunzip chr21.maf.gz

maf_parse --features input_rs.chr21.txt chr21.maf > input_rs.chr21.multiz

rm chr21.maf.gz

rm chr21.maf`
