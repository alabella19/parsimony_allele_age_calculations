#!/usr/bin/perl
use lib '/home/labellal/lib';
#split multiz file into the corresponding SNPS
use Data::Dumper;
use POSIX;
#use Data::Dumper::Sortkeys = 1;

##This script takes an xmfa file and splits it so that parsimony test can be run on it. 

if($#ARGV<0){
    print "*******************************************\nSyntax: xmfa_parsimony.pl xmfa_file\n\n";
    exit;
}

$xmfa_file = $ARGV[0];

@xmfa_num = split('\.',$xmfa_file);
$xmfa_num = $xmfa_num[0];


open (INPUT, $xmfa_file);
@xmfa = <INPUT>;
close(INPUT);
@this_fasta=();
$first = 0;
foreach $line (@xmfa){
	if($line =~ /\=/){
		if($first != 0){
			open FILE , ">", "$this_file";
			print FILE @this_fasta;
			close FILE;
			
			`Rscript parsimony_date_one.R $this_file $xmfa_num`;
			`rm $this_file`;
		}
		@this_fasta=();
		$line =~ s/\n//g;
		$this_file=$line;
		$this_file =~ s/\= //g;
		$first = 1;
	}else{
		push(@this_fasta, $line);
	}
	

	
	
	

	
	#open FILE, ">>", "$file_num.all_fasta.xmfa";
	#print FILE "= $out_file2\n";
	#print FILE @out_array2;
	#close FILE;
	#
	#$rs_total=$rs_total+1;

}
open FILE , ">", "$this_file";
print FILE @this_fasta;
close FILE;

`Rscript parsimony_date_one.R $this_file $xmfa_num`;
`rm $this_file`;

#all the species in the dataset
