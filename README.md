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

input_rs.rs.rev_data.txt - information needed for reverse data in step 5

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
```bash
wget -nd -q "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/maf/chr21.maf.gz"
gunzip chr21.maf.gz
maf_parse --features input_rs.chr21.txt chr21.maf > input_rs.chr21.multiz
rm chr21.maf.gz
rm chr21.maf
```

# Step 3 - Get Neanderthal data. This is optional but better to do it now than want it later! 

### Requirements
VCFtools https://vcftools.github.io/index.html

### Input
hg19 information - example: input_rs.rs_hg19_chr21.txt

### Output
VCF file of the neanderthal information for each neanderthal genome - example: input_rs.chr21.vindija.out.recode.vcf

Log for each file - example: input_rs.chr21.vindija.out.log

### Execution 
Get the neanderthal data and rename it

```bash
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_altai.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_denisova.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/LBK/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_lbk.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Loschbour/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_loschbour.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Mez1/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_mez1.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Ust_Ishim/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_ishim.vcf.gz
wget -nd -q "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz"
mv chr21_mq25_mapab100.vcf.gz chr21_vindija.vcf.gz
```

Extract just the relevent positions

```bash
vcftools --gzvcf chr21_altai.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.altai.out
vcftools --gzvcf chr21_denisova.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.denisova.out
vcftools --gzvcf chr21_lbk.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.lbk.out
vcftools --gzvcf chr21_loschbour.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.loschbour.out
vcftools --gzvcf chr21_mez1.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.mez1.out
vcftools --gzvcf chr21_ishim.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.ishim.out
vcftools --gzvcf chr21_vindija.vcf.gz --positions input_rs.rs_hg19_chr21.txt --recode --out input_rs.chr21.vindija.out
```
# Step 4 - obtain greatape data
### Requirements
VCFtools https://vcftools.github.io/index.html

### Input
hg18 information - example: input_rs.rs_hg18_21.txt

### Output
VCF for each each species - example: rs_input.chr21.Gorilla.vcf.recode.vdf

Log for each species - example: rs_input.chr21.Gorilla.vcf.log

### Execution 
Get great ape data - these files are LARGE and do not include indels

```bash
wget -nd -q "https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Gorilla.vcf.gz"
wget -nd -q "https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_paniscus.vcf.gz"
wget -nd -q "https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_troglodytes.vcf.gz"
wget -nd -q "https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_abelii.vcf.gz"
wget -nd -q "https://eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_pygmaeus.vcf.gz"
```

Extract positions
```bash
vcftools --gzvcf Gorilla.vcf.gz --positions input_rs.rs_hg18_chr21.txt --recode --out input_rs.chr21.Gorilla.vcf
vcftools --gzvcf Pan_paniscus.vcf.gz --positions input_rs.rs_hg18_chr21.txt --recode --out input_rs.chr21.Pan_paniscus.vcf
vcftools --gzvcf Pan_troglodytes.vcf.gz --positions input_rs.rs_hg18_chr21.txt --recode --out input_rs.chr21.Pan_troglodytes.vcf
vcftools --gzvcf Pongo_abelii.vcf.gz --positions input_rs.rs_hg18_chr21.txt --recode --out input_rs.chr21.Pongo_abelii.vcf
vcftools --gzvcf Pongo_pygmaeus.vcf.gz --positions input_rs.rs_hg18_chr21.txt --recode --out input_rs.chr21.Pongo_pygmaeus.vcf
```



# Step 5 - Get SNPs that may be listed in the reverse orientation
Due to differences between the databases some SNPs are listed on the - strand, while others are listed on the + strand.

To account for these issues we used the VCF for hg38

### Requirements
VCFtools https://vcftools.github.io/index.html

### Input
positions formatted for VCFtools - example: input_rs.rev_data.txt 

VCF of all the SNPs with orientation information from hg38 - LARGE FILE

### Output
VCF of the reverse information - example: input_rs.snpdb.vcf.recode.vcf

Log of reverse information - example: input_rs.snpdb.vcf.log

### Execution
Get the reference database (LARGE)
```bash
wget -nd -q "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz"
```

Get the reverse data
```bash
vcftools --gzvcf All_20180418.vcf.gz --positions input_rs.rev_data.txt --recode --recode-INFO "RV" --out input_rs.snpdb.vcf
```

# Step 6 - Get the allele data for indels 

Due to the differences in databases indels are sometimes listed in different ways. For example A- and AT vs -G and TG for the same SNP 

### Requirements
VCFtools https://vcftools.github.io/index.html

### Input
1000 genomes VCF file http://www.internationalgenome.org/faq/are-1000-genomes-variant-calls-phased/

hg19 snp positions - example: input_rs.rs_hg19_chr21.txt

### Output
VCF with only the positions of interest - example: input_rs.chr21.recode.vcf

Files with the 3 different types of SNPs

No indel - example: input_rs.chr21.alleles

Indels - examples: input_rs.chr21.gap1.alleles, input_rs.chr21.gap2.alleles

### Execution

```bash
vcftools --gzvcf ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions input_rs.rs_hg19_chr21.txt --indv HG00103 --recode --out input_rs.chr21
```

NOTE: filtered out on a random individual

Separate the indels

```bash
grep -o -P '^\d+\t\d+\trs\d+\t.\t.\t' input_rs.chr21.recode.vcf > input_rs.chr21.alleles
grep -o -P '^\d+\t\d+\trs\d+\t..\t.\t' input_rs.chr21.recode.vcf > input_rs.chr21.gap1.alleles
grep -o -P '^\d+\t\d+\trs\d+\t.\t..\t' input_rs.chr21.recode.vcf > input_rs.chr21.gap2.alleles
```

NOTE: sequences without alleles in hg19 WILL NOT BE ANALYZED 



# Step 7 - Combine all the data! 

### Requirements
Perl - Data::Dumper & POSIX

### Input
Order of the input files

1 - multiz alignment produced in Step 2

2 - hg38 position information (.txt)

3 - allelic data generated in Step 6 (will automatically find the gap1 and gap2 files based on files generated in step 6)

4 - hg18 reference data (.ref) produced in Step 1

5 - hg19 reference data (.ref) produced in Step 1

6 - header for neanderthal data generated in Step 3 

7 - header for the great ape data generated in Step 4

8 - Reverse data generated in Step 5 (VCF)


### Output
The script will print the SNPs not included in the analysis and the reason why

The script will create an xmfa for EACH ALLELE for the variant positions provided. The xmfa files will be split into a manageable size


### Execution 
```bash
split_multiz_v2.pl input_rs.chr21.multiz input_rs.chr21.txt input_rs.chr21.alleles input_rs.rs_hg18_chr21.ref input_rs.rs_hg19_chr21.ref input_rs.chr21 input_rs.chr21 input_rs.snpdb.vcf.recode.vcf >multiz.out
```

