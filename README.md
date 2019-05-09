# parsimony_allele_age_calculations
This repository contains the scripts and instructions for calculating parsimony based age calculations of SNP alleles from mammalian alignments

---
# Step 1 - obtaining hg18, hg19 and hg38 positions for all loci based on RS number

NOTE: This step will remove any loci that are indels longer than 1 base pair. 

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
`wget -nd -q "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/maf/chr21.maf.gz"`

`gunzip chr21.maf.gz`

`maf_parse --features input_rs.chr21.txt chr21.maf > input_rs.chr21.multiz`

`rm chr21.maf.gz`

`rm chr21.maf`

# Step 3 - Get Neanderthal data. This is optional but better to do it now than want it later! 

### Requirements
VCFtools https://vcftools.github.io/index.html

### Input
hg19 information - example: input_rs.rs_hg19_chr21.txt

### Output
VCF file of the neanderthal information for each neanderthal genome

### Execution 
Get the neanderthal data and rename it

`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_altai.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_denisova.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/LBK/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_lbk.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Loschbour/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_loschbour.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Mez1/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_mez1.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Ust_Ishim/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_ishim.vcf.gz`
`wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz"`
`mv chr21_mq25_mapab100.vcf.gz chr21_vindija.vcf.gz`

Extract just the relevent positions

`vcftools --gzvcf chr21_altai.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.altai.out`
`vcftools --gzvcf chr21_denisova.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.denisova.out`
`vcftools --gzvcf chr21_lbk.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.lbk.out`
`vcftools --gzvcf chr21_loschbour.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.loschbour.out`
`vcftools --gzvcf chr21_mez1.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.mez1.out`
`vcftools --gzvcf chr21_ishim.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.ishim.out`
`vcftools --gzvcf chr21_vindija.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.vindija.out`


