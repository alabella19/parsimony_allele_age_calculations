.libPaths("~/R/rlib-3.4.3")
library(biomaRt)
library(readr)

##Get the file with the RS values only
args = commandArgs(trailingOnly=TRUE)
rs_file = args[1]
rs_out_file<-gsub("txt","",rs_file)


suppressMessages(rs_input<-read_delim(rs_file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names=FALSE))

#set up the biomart
variation=useEnsembl(biomart="snp", dataset="hsapiens_snp", host="useast.ensembl.org")


rs_hg38<-getBM(filters='snp_filter',values=unlist(rs_input),mart=variation, attributes=c('refsnp_id','chr_name','chrom_start','chrom_end','ensembl_transcript_chrom_strand','allele'))


#FOR the VCF file we need the second/lower number for ALL SNPs with indels
#For the multiz alignment we want the higher number! 
#For example see SNP rs201650033
#In hg38 it is listed at 8837543 - this is the position we want from the multiz alignment
#in the VCF file it is listed at 8837542 - this is the position we want for the RV data


rv_data<-rs_hg38
to_subtract<-grep('-',rv_data$allele)
indels<-rv_data[to_subtract,]
no_indels<-rv_data[-(to_subtract),]
indels$chrom_end<-indels$chrom_start-1

rv_data<-rbind(indels,no_indels)

write.table(rv_data[,c(2,4)],file=paste(rs_out_file,"rev_data.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)

#print(rs_hg38)

#get removed due to missing from biomart
missing<-rs_input[which(rs_input$X1 %in% rs_hg38$refsnp_id == FALSE),]

#removed due to having designation on an alternative chromosome
alternative_chromosome<-rs_hg38[grep("CHR", rs_hg38$chr_name),]

#new subset
if(nrow(alternative_chromosome)>0){
	rs_hg38<-rs_hg38[-grep("CHR",rs_hg38$chr_name),]
	}


#remove those with length >0
rs_hg38$diff=rs_hg38$chrom_end - rs_hg38$chrom_start

#If the difference is -1 then the reference allele in the databse is a gap
#Use first number

rs_saved<-rs_hg38

rs_hg38<-subset(rs_hg38, diff > -2)



gff<-(rs_hg38[,1:2])
gff<-cbind(gff,rep("Ensembl",nrow(gff)))
gff<-gff[,-1]
gff<-cbind(gff, rs_hg38[,1])
gff<-cbind(gff, rs_hg38[,4])
gff<-cbind(gff, rs_hg38[,4])
gff<-cbind(gff,rep(".",nrow(gff)))
gff<-cbind(gff,rep(".",nrow(gff)))
gff<-cbind(gff,rep(".",nrow(gff)))
gff<-cbind(gff,rep(".",nrow(gff)))



colnames(gff)<-c("chr","source","rs","start","end","d1","d2","d3","d4")
#print(gff)
#split into chromosomes 1-22 & X

write.table(subset(gff, chr==1),file=paste(rs_out_file,"chr1.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==2),file=paste(rs_out_file,"chr2.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==3),file=paste(rs_out_file,"chr3.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==4),file=paste(rs_out_file,"chr4.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==5),file=paste(rs_out_file,"chr5.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==6),file=paste(rs_out_file,"chr6.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==7),file=paste(rs_out_file,"chr7.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==8),file=paste(rs_out_file,"chr8.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==9),file=paste(rs_out_file,"chr9.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==10),file=paste(rs_out_file,"chr10.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==11),file=paste(rs_out_file,"chr11.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==12),file=paste(rs_out_file,"chr12.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==13),file=paste(rs_out_file,"chr13.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==14),file=paste(rs_out_file,"chr14.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==15),file=paste(rs_out_file,"chr15.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==16),file=paste(rs_out_file,"chr16.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==17),file=paste(rs_out_file,"chr17.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==18),file=paste(rs_out_file,"chr18.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==19),file=paste(rs_out_file,"chr19.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==20),file=paste(rs_out_file,"chr20.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==21),file=paste(rs_out_file,"chr21.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr==22),file=paste(rs_out_file,"chr22.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr=="X"),file=paste(rs_out_file,"chrX.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(gff, chr=="Y"),file=paste(rs_out_file,"chrY.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)


#need to get hg19 for Ancient human
#print("past hg38")
#
mart19<-useMart(host = "grch37.ensembl.org",biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
rs_hg19<-getBM(filters='snp_filter',values=unlist(rs_input),mart=mart19, attributes=c('refsnp_id','chr_name','chrom_start','allele'))
missing_hg19<-rs_input[which(rs_input$X1 %in% rs_hg19$refsnp_id == FALSE),]


#Need subtract 1 from SNPs with an insertion or deletion 
to_subtract<-grep('-',rs_hg19$allele)
indels<-rs_hg19[to_subtract,]
no_indels<-rs_hg19[-(to_subtract),]
indels$chrom_start<-indels$chrom_start-1

rs_hg19<-rbind(indels,no_indels)

rs_hg19<-rs_hg19[,-4]

write.table(subset(rs_hg19, chr_name=="1"),file=paste(rs_out_file,"rs_hg19_chr1.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="2"),file=paste(rs_out_file,"rs_hg19_chr2.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="3"),file=paste(rs_out_file,"rs_hg19_chr3.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="4"),file=paste(rs_out_file,"rs_hg19_chr4.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="5"),file=paste(rs_out_file,"rs_hg19_chr5.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="6"),file=paste(rs_out_file,"rs_hg19_chr6.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="7"),file=paste(rs_out_file,"rs_hg19_chr7.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="8"),file=paste(rs_out_file,"rs_hg19_chr8.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="9"),file=paste(rs_out_file,"rs_hg19_chr9.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="10"),file=paste(rs_out_file,"rs_hg19_chr10.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="11"),file=paste(rs_out_file,"rs_hg19_chr11.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="12"),file=paste(rs_out_file,"rs_hg19_chr12.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="13"),file=paste(rs_out_file,"rs_hg19_chr13.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="14"),file=paste(rs_out_file,"rs_hg19_chr14.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="15"),file=paste(rs_out_file,"rs_hg19_chr15.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="16"),file=paste(rs_out_file,"rs_hg19_chr16.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="17"),file=paste(rs_out_file,"rs_hg19_chr17.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="18"),file=paste(rs_out_file,"rs_hg19_chr18.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="19"),file=paste(rs_out_file,"rs_hg19_chr19.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="20"),file=paste(rs_out_file,"rs_hg19_chr20.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="21"),file=paste(rs_out_file,"rs_hg19_chr21.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="22"),file=paste(rs_out_file,"rs_hg19_chr22.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="X"),file=paste(rs_out_file,"rs_hg19_chrX.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="Y"),file=paste(rs_out_file,"rs_hg19_chrY.ref", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)


##print(rs_hg19)
#
rs_hg19<-cbind(rep("chr",nrow(rs_hg19)),rs_hg19)
rs_hg19<-rs_hg19[,-2]
#rs_hg19$chr_name<-paste(rs_hg19$`rep("chr", nrow(rs_hg19))`,rs_hg19$chr_name, sep = "")
rs_hg19<-rs_hg19[,-1]

#print(rs_hg19)




write.table(subset(rs_hg19, chr_name=="1"),file=paste(rs_out_file,"rs_hg19_chr1.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="2"),file=paste(rs_out_file,"rs_hg19_chr2.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="3"),file=paste(rs_out_file,"rs_hg19_chr3.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="4"),file=paste(rs_out_file,"rs_hg19_chr4.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="5"),file=paste(rs_out_file,"rs_hg19_chr5.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="6"),file=paste(rs_out_file,"rs_hg19_chr6.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="7"),file=paste(rs_out_file,"rs_hg19_chr7.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="8"),file=paste(rs_out_file,"rs_hg19_chr8.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="9"),file=paste(rs_out_file,"rs_hg19_chr9.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="10"),file=paste(rs_out_file,"rs_hg19_chr10.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="11"),file=paste(rs_out_file,"rs_hg19_chr11.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="12"),file=paste(rs_out_file,"rs_hg19_chr12.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="13"),file=paste(rs_out_file,"rs_hg19_chr13.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="14"),file=paste(rs_out_file,"rs_hg19_chr14.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="15"),file=paste(rs_out_file,"rs_hg19_chr15.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="16"),file=paste(rs_out_file,"rs_hg19_chr16.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="17"),file=paste(rs_out_file,"rs_hg19_chr17.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="18"),file=paste(rs_out_file,"rs_hg19_chr18.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="19"),file=paste(rs_out_file,"rs_hg19_chr19.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="20"),file=paste(rs_out_file,"rs_hg19_chr20.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="21"),file=paste(rs_out_file,"rs_hg19_chr21.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="22"),file=paste(rs_out_file,"rs_hg19_chr22.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="X"),file=paste(rs_out_file,"rs_hg19_chrX.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg19, chr_name=="Y"),file=paste(rs_out_file,"rs_hg19_chrY.txt", sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)


#
#get hg18 for greatape

mart18<-useMart(host = "may2009.archive.ensembl.org",biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
rs_hg18<-getBM(filters='refsnp',values=unlist(rs_input),mart=mart18, attributes=c('refsnp_id','chr_name','chrom_start','allele'))

to_subtract<-grep('-',rs_hg18$allele)

indels<-rs_hg18[to_subtract,]
no_indels<-rs_hg18[-(to_subtract),]
indels$chrom_start<-indels$chrom_start-1

rs_hg18<-rbind(indels,no_indels)


#print(rs_hg18)
missing_hg18<-rs_input[which(rs_input$X1 %in% rs_hg18$refsnp_id == FALSE),]
alternative_chromosome_hg18<-rs_hg18[grep("NT", rs_hg18$chr_name),]
alternative_chromosome2_hg18<-rs_hg18[grep("c", rs_hg18$chr_name),]

if(nrow(alternative_chromosome_hg18)>0){
	rs_hg18<-rs_hg18[-grep("NT",rs_hg18$chr_name),]
}
if(nrow(alternative_chromosome2_hg18)>0){
	rs_hg18<-rs_hg18[-grep("c",rs_hg18$chr_name),]
}

rs_hg18<-rs_hg18[,-4]

write.table(subset(rs_hg18, chr_name==1),file=paste(rs_out_file,"rs_hg18_chr1.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==2),file=paste(rs_out_file,"rs_hg18_chr2.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==3),file=paste(rs_out_file,"rs_hg18_chr3.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==4),file=paste(rs_out_file,"rs_hg18_chr4.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==5),file=paste(rs_out_file,"rs_hg18_chr5.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==6),file=paste(rs_out_file,"rs_hg18_chr6.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==7),file=paste(rs_out_file,"rs_hg18_chr7.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==8),file=paste(rs_out_file,"rs_hg18_chr8.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==9),file=paste(rs_out_file,"rs_hg18_chr9.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==10),file=paste(rs_out_file,"rs_hg18_chr10.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==11),file=paste(rs_out_file,"rs_hg18_chr11.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==12),file=paste(rs_out_file,"rs_hg18_chr12.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==13),file=paste(rs_out_file,"rs_hg18_chr13.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==14),file=paste(rs_out_file,"rs_hg18_chr14.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==15),file=paste(rs_out_file,"rs_hg18_chr15.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==16),file=paste(rs_out_file,"rs_hg18_chr16.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==17),file=paste(rs_out_file,"rs_hg18_chr17.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==18),file=paste(rs_out_file,"rs_hg18_chr18.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==19),file=paste(rs_out_file,"rs_hg18_chr19.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==20),file=paste(rs_out_file,"rs_hg18_chr20.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==21),file=paste(rs_out_file,"rs_hg18_chr21.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name==22),file=paste(rs_out_file,"rs_hg18_chr22.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="X"),file=paste(rs_out_file,"rs_hg18_chrX.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="Y"),file=paste(rs_out_file,"rs_hg18_chrY.ref",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)

rs_hg18<-cbind(rep("chr",nrow(rs_hg18)),rs_hg18)
rs_hg18<-rs_hg18[,-2]
rs_hg18$chr_name<-paste(rs_hg18$`rep("chr", nrow(rs_hg18))`,rs_hg18$chr_name, sep = "")
rs_hg18<-rs_hg18[,-1]

write.table(subset(rs_hg18, chr_name=="chr1"),file=paste(rs_out_file,"rs_hg18_chr1.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr2"),file=paste(rs_out_file,"rs_hg18_chr2.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr3"),file=paste(rs_out_file,"rs_hg18_chr3.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr4"),file=paste(rs_out_file,"rs_hg18_chr4.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr5"),file=paste(rs_out_file,"rs_hg18_chr5.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr6"),file=paste(rs_out_file,"rs_hg18_chr6.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr7"),file=paste(rs_out_file,"rs_hg18_chr7.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr8"),file=paste(rs_out_file,"rs_hg18_chr8.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr9"),file=paste(rs_out_file,"rs_hg18_chr9.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr10"),file=paste(rs_out_file,"rs_hg18_chr10.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr11"),file=paste(rs_out_file,"rs_hg18_chr11.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr12"),file=paste(rs_out_file,"rs_hg18_chr12.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr13"),file=paste(rs_out_file,"rs_hg18_chr13.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr14"),file=paste(rs_out_file,"rs_hg18_chr14.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr15"),file=paste(rs_out_file,"rs_hg18_chr15.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr16"),file=paste(rs_out_file,"rs_hg18_chr16.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr17"),file=paste(rs_out_file,"rs_hg18_chr17.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr18"),file=paste(rs_out_file,"rs_hg18_chr18.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr19"),file=paste(rs_out_file,"rs_hg18_chr19.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr20"),file=paste(rs_out_file,"rs_hg18_chr20.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr21"),file=paste(rs_out_file,"rs_hg18_chr21.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chr22"),file=paste(rs_out_file,"rs_hg18_chr22.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chrX"),file=paste(rs_out_file,"rs_hg18_chrX.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)
write.table(subset(rs_hg18, chr_name=="chrY"),file=paste(rs_out_file,"rs_hg18_chrY.txt",sep=""),quote=FALSE, sep="\t", col.names=FALSE, row.names = FALSE)



