#!/usr/bin/perl
use lib '/home/labellal/lib';
#split multiz file into the corresponding SNPS
use Data::Dumper;
use POSIX;
#use Data::Dumper::Sortkeys = 1;

#V2 now uses the POSITION to get the allele data

print "running split multiz\n";

if($#ARGV<7){
    print "*******************************************\nSyntax: split_multiz.pl multiz_from_ucsc\n input_ucsc\n snpdata.txt\n hg18_snp_data.txt\n hg19_snp_data.txt\n greatapebase\n neanderthalbase\n reverse_file.vcf\n\n\n";
    exit;
}

$multiz_file = $ARGV[0];
$input_file = $ARGV[1];
$snp_data_file = $ARGV[2];
#$ancient_file = $ARGV[3];
$hg18_file = $ARGV[3];
$hg19_file=$ARGV[4];
$greatapebase=$ARGV[5];
$nbase=$ARGV[6];
$rev_file=$ARGV[7];


open (INPUT, $multiz_file) or die "Can't open < $multiz_file $!";;
@multiz = <INPUT>;
close(INPUT);

open (INPUT, $input_file) or die "Can't open < $input_file $!";;
@input = <INPUT>;
close(INPUT);

open(INPUT, $snp_data_file) or die "Can't open < $snp_data_file $!";;
@snp_data = <INPUT>;
close(INPUT);

#open (INPUT, $ancient_file);
#@ancient_data = <INPUT>;
#close(INPUT);

open(INPUT, $hg18_file) or die "Can't open < hg18 $hg18_file $!";;
@hg18 = <INPUT>;
close(INPUT);

open(INPUT, $hg19_file) or die "Can't open < hg19 $hg19_file $!";;
@hg19=<INPUT>;
close(INPUT);

##Get the info for the SNPs containing a GAP
$allele_base=$snp_data_file;
$allele_base=~ s/alleles//;
$allele_gap1_file=$allele_base."gap1.alleles";
$allele_gap2_file=$allele_base."gap2.alleles";

open(INPUT, $allele_gap1_file) or die "Can't open < $allele_gap1_file $!";;
@allele_gap1=<INPUT>;
close(INPUT);

open(INPUT, $allele_gap2_file) or die "Can't open < $allel_gap2_file $!";;
#two bases in second pos
@allele_gap2=<INPUT>;
close(INPUT);


$amb{"C"}{"T"}="Y";
$amb{"C"}{"G"}="S";
$amb{"C"}{"A"}="M";
$amb{"A"}{"T"}="W";
$amb{"A"}{"C"}="M";
$amb{"A"}{"G"}="R";
$amb{"T"}{"A"}="W";
$amb{"T"}{"C"}="Y";
$amb{"T"}{"G"}="K";
$amb{"G"}{"A"}="R";
$amb{"G"}{"T"}="K";
$amb{"G"}{"C"}="S";



foreach $line (@snp_data){
	$line=~ s/\n//g;
	@line= split('\t',$line);
	$pos = $line[1];
	$a1 = $line[3];
	$a2 = $line[4];
	if(exists $allele1{$pos}){
		print "more than one SNP at position $pos... removing\n";
		$allele1{$pos}="";
		$allele2{$pos}="";
	}else{
		$allele1{$pos}=$a1;
		$allele2{$pos}=$a2;
	}
}
#print"checkpoint1\n";
foreach $line (@allele_gap1){
	$line=~ s/\n//g;
	@line= split('\t',$line);
	$pos = $line[1];
	$a1 = $line[3];
	$a2 = $line[4];
	#two bases in first allele
	#print"a1:$a1\ta2:$a2\n";
	if(exists $allele1{$pos}){
		print "more than one SNP at position $pos... removing\n";
		$allele1{$pos}="";
		$allele2{$pos}="";
	}else{
		@a1=split('',$a1);
		
		#print Dumper\@a1;
		if($a1[0] eq $a2){
			$a1=$a1[1];	
			$allele1{$pos}=$a1;
			$allele2{$pos}='-';
		}elsif($a1[1] eq $a2){
			$a1=$a1[0];
			$allele1{$pos}=$a1;
			$allele2{$pos}='-';
		}else{
			print "rs $pos has alleles NN -: SKIPPED\n";	
		}
		
		#print"a1:$a1\ta2:$a2\n";
	}

}

#print"checkpoint2\n";
%all_gap2;
foreach $line (@allele_gap2){
	#GAP 2 NEED SPECIAL ATTENTION!!! 
	#Need to get only the position that has the gap
	$line=~ s/\n//g;
	@line= split('\t',$line);
	$pos = $line[1];
	$a1 = $line[3];
	$a2 = $line[4];
	#two bases in first allele
	#print"a1:$a1\ta2:$a2\n";
	
	if(exists $allele1{$pos}){
		print "more than one SNP at position $pos... removing\n";
		$allele1{$pos}="";
		$allele2{$pos}="";
	}else{
		@a2=split('',$a2);
		
		#print Dumper\@a1;
		if($a2[0] eq $a1){
			$a2=$a2[1];	
			$allele2{$pos}=$a2;
			$allele1{$pos}='-';
			$all_gap2{$pos}=1;
		}elsif($a2[1] eq $a1){
			$a2=$a2[0];
			$allele2{$pos}=$a2;
			$allele1{$pos}='-';
			$all_gap2{$pos}=1;
		}else{
			print "rs $pos has alleles - NN: SKIPPED\n";	
		}
		
		#print"a1:$a1\ta2:$a2\n";
	}

}

#print Dumper \%allele1;
#print Dumper \%allele2;
#
#exit;

#print"checkpoint2\n";
##NEED TO FIX THESE!!!!!!
#shift(@hg18);
foreach $line (@hg18){
	$line =~ s/\n//g;
	@pos = split("\t",$line);
	$name = $pos[0];
	$pos = $pos[1]."\t".$pos[2];
	$pos = "chr".$pos;
	$snp_pos_info{$pos}=$name;
}
#print"checkpoint3\n";
#shift(@hg19);
foreach $line (@hg19){
	#print $line;
	$line =~ s/\n//g;
	@pos = split("\t",$line);
	$name = $pos[0];
	
	$snp_rs_info19{$name}=$pos[2];
	
	$pos = $pos[1]."\t".$pos[2];
	#$pos = "chr".$pos;
	$snp_pos_info19{$pos}=$name;
	
	
}
#print "hg19 info\n";
#print Dumper \%snp_pos_info19;
#exit;
#print Dumper \%snp_info;


@all_rs = ();


#print"checkpoint4\n";
foreach $line (@input){
	#print $line;
	$line =~ s/\n//g;
	@line = split('\t', $line);
	$rs = $line[2];
	$pos = $line[3];
	#print "rs is $rs and pos is $pos\n";
	#When looking in the multiz file the numbers are 0 offset SO we need to look for one less than what
	#is in this file
	$pos = $pos - 1;
	$multiz_pos{$pos}=$rs;
	push(@all_rs, $rs);
}
#print "multiz pos\n";
#print Dumper \%multiz_pos;
#print Dumper \@all_rs;
#exit;
#print"checkpoint5\n";
#print Dumper \@all_rs;
$n = 0;

foreach $line (@multiz){
	if($line =~ /score/){
		#start of a new alignment
	}else{
		if($line =~ /^s/){
			if($line =~ /hg38/){
				#Get RS for this 
				@line = split('\s+',$line);
				$this_pos = $line[2];
				$this_rs=$multiz_pos{$this_pos};
				#if($this_rs =~ /rs4965031/){
				#	print $this_rs;
				#	print Dumper \@line;
				#	exit;
				#}
				$seq = $line[6];
				$spec = $line[1];
				@spec = split('\.',$spec);
				$spec = $spec[0];
				$spec =~ s/\t//g;
				$fasta = ">$spec\n$seq\n";
				$len_hash{$this_rs} = length($seq);
				$all_fasta{$this_rs}{$spec}=$fasta;
					
			}else{
				@line = split('\s+',$line);
				$seq = $line[6];
				$spec = $line[1];
				@spec = split('\.',$spec);
				$spec = $spec[0];
				$spec =~ s/\t//g;
				$fasta = ">$spec\n$seq\n";
				$len_hash{$this_rs} = length($seq);
				$all_fasta{$this_rs}{$spec}=$fasta;
			}
			#exit;
		}elsif($line=~/^e/){
			#this species is alignable BUT has a gap 
			#if value is C = gap
			#if value is I = N
			#if value is M = N
			#if value is n = n
			$line =~ s/\n//g;
			@line=split('\s+',$line);
			$spec=$line[1];
			@spec = split('\.',$spec);
			$spec = $spec[0];
			$spec =~ s/\t//g;
			$seq=$line[6];
			if($seq eq "C"){
				$fasta=">$spec\n-\n";
				$all_fasta{$this_rs}{$spec}=$fasta;
			}
		}
		
	}
}
#print "All fastas\n";
#print Dumper \%all_fasta;
#print"checkpoint6\n";
$reverse_hash{"G"}="C";
$reverse_hash{"g"}="c";
$reverse_hash{"C"}="G";
$reverse_hash{"c"}="g";
$reverse_hash{"T"}="A";
$reverse_hash{"t"}="a";
$reverse_hash{"A"}="T";
$reverse_hash{"a"}="t";
$reverse_hash{"N"}="N";
$reverse_hash{"-"}="-";


#SKIPPING NEANDERTHAL DAT
@all_neand=("altai","denisova","lbk","loschbour","mez1","ishim","vindija");
#change this to use VCF files
foreach $species (@all_neand){
	#print $species;
	$neand_file = $nbase.".".$species.".out.recode.vcf";
	open (INPUT, $neand_file) or die "Can't open < $neand_file $!";
	@nd_input = <INPUT>;
	close(INPUT);
	
	foreach $line (@nd_input){
		if ($line !~ /^#/){
			#print $line;
			@line = split('\t',$line);
			$pos = $line[0]."\t".$line[1];
			$this_rs = $snp_pos_info19{$pos};
			$ref = $line[3];
			$alt = $line[4];
			if($alt eq "."){
				$this_base = $ref;	
			}else{
				$this_base=$amb{$ref}{$alt};
			}
			#$new_name = $great_ape_names{$species};
			$nd_data{$this_rs}{$species}=$this_base;
		}
	}
}

#print Dumper \%nd_data;
#exit;


#print Dumper \%len_hash;
@all_spec=("panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1","cerAty1","papAnu3","chlSab2","manLeu1","nasLar1","colAng1","rhiRox1","rhiBie1","calJac3","saiBol1","cebCap1","aotNan1","tarSyr2","micMur3","proCoq1","eulMac1","eulFla1","otoGar3","mm10","canFam3","dasNov3","Pongo_pygmaeus");

#print"checkpoint7\n";

$great_ape_names{"Gorilla"}="gorGor5";
$great_ape_names{"Pan_paniscus"}="panPan2";
$great_ape_names{"Pan_troglodytes"}="panTro5";
$great_ape_names{"Pongo_abelii"}="ponAbe2";
$great_ape_names{"Pongo_pygmaeus"}="Pongo_pygmaeus"; #there is no pongo_pygmaeus in the greatape dataset


#print Dumper \%snp_pos_info;
#get all of the greatape data
foreach $species (keys %great_ape_names){
	#print $species;
	open (INPUT, "$greatapebase.$species.vcf.recode.vcf") or die "Can't open < $greatapebase.$species.vcf.recode.vcf: $!";
	@ga_input = <INPUT>;
	close(INPUT);
	
	foreach $line (@ga_input){
		if ($line !~ /^#/){
			#print $line;
			@line = split('\t',$line);
			$pos = $line[0]."\t".$line[1];
			$this_rs = $snp_pos_info{$pos};
			$ref = $line[3];
			$alt = $line[4];
			$this_base=$amb{$ref}{$alt};
			$new_name = $great_ape_names{$species};
			$ga_data{$this_rs}{$new_name}=$this_base;
		}
	}
	

}

	#print Dumper \%ga_data;
	#exit;


##print"checkpoint8\n";
#@all_rs = keys %all_fasta;
#print "all rs is ";
#print Dumper \@all_rs;
#exit;
	$rs_total=0;
foreach $rs (keys %all_fasta){
	$RV_this=0;
	@out_fasta=();
	##NEED to make an out hash
	#print "processing rs $rs\n";
	#print"checkpoint8a\n";
	#need to check to see if there is population level information for the greatape data
	if ($rs eq ''){
		next;	
	}
	
	
	foreach $spec (@all_spec){
		$current_val=$all_fasta{$rs}{$spec};
		#print "current: $rs $spec $current_val:\n";
		if(exists $ga_data{$rs}{$spec}){
			#USE the greatape data instead of the mutliz data
			#print"using ga data for $rs\n";
			$base = $ga_data{$rs}{$spec};
			$out_hash{$spec}=$base;
			#push(@out_fasta, ">$spec\n$base\n");
		}
		elsif(exists $all_fasta{$rs}{$spec}){
			$seq = $all_fasta{$rs}{$spec};
			@seq = split('\n',$seq);
			$seq = $seq[1];
			$out_hash{$spec}=$seq;
			#push(@out_fasta, $all_fasta{$rs}{$spec});
		}else{
			#print "replacing with N for $rs $spec\n";
			$len = $len_hash{$rs};
			$seq = "N"x$len;
			$out_hash{$spec}=$seq;
			#push(@out_fasta, ">$spec\n$seq\n");
		}
	}
	
	#NEW out hash;
	#print Dumper \%out_hash;
	#add in neanderthal data
	#print"checkpoint8b\n";
	foreach $spec (@all_neand){
		if(exists $nd_data{$rs}{$spec}){
			$base = $nd_data{$rs}{$spec};
			$out_hash{$spec}=$base;
			#push(@out_fasta, ">$spec\n$base\n");
		}else{
			$len = $len_hash{$rs};
			$seq = "N"x$len;
			$out_hash{$spec}=$seq;
			#push(@out_fasta, ">$spec\n$seq\n");
		}
	}
	#print "after neander\n";
	#print Dumper \%out_hash;
	#check for Pon Pygm
	#if(exists $ga_data{$rs}{"Pongo_pygmaeus"}){
	#	$base = $ga_data{$rs}{$spec};
	#	push(@out_fasta, ">Pongo_pygmaeus\n$base\n");
	#}
	
	##add in the ancestor
	#print"checkpoint8c\n";
	
	#ADD HUMAN DATA HERE
	#ONE for each allele
	#check to make sure at least ONE allele matches what's in the allele file
	$this_pos=$snp_rs_info19{$rs};
	
	$this_allele1=$allele1{$this_pos};
	$this_allele2=$allele2{$this_pos};
	
	$human = $all_fasta{$rs}{"hg38"};
	@human = split('\n',$human);
	$human = $human[1];
	
	$human_test=$human;
	$human_test =~ s/\-//g;
	
	#print"checkpoint8d\n";
	if(exists $all_gap2{$this_pos}){
		$human_test = "-";
		#print "is an all gap2\n";
	}
	#print "human is $human\n";
	#print "this $rs comparing $human_test to $this_allele1 and $this_allele2\n";
	if(lc $this_allele1 eq lc $human_test){
		#allele matches
		$match=1;
	}elsif(lc $this_allele2 eq lc $human_test){
		$match=2	
	}else{
		$match=0;
		$rev = `grep -w $rs $rev_file`;
		#print $rev;
		@rev=split('\t',$rev);
		$rev = $rev[7];
		$rev =~ s/\n//g;
		#print "rev is:$rev\n";
		if($rev eq "RV\=1"){
			print "need to reverse $rs\n";
			print "for $rs changing $this_allele1 and $this_allele2\n";
			$RV_this=1;
			$new_allele2 = $reverse_hash{$this_allele2};
			$new_allele1 = $reverse_hash{$this_allele1};
			
			#print Dumper \%out_hash;
			#print "New alleles are $new_allele1 and $new_allele2\n";
			$this_allele1=$new_allele1;
			$this_allele2=$new_allele2;
			
			print "for $rs replacing $this_allele1 and $this_allele2\n";
			if(lc $new_allele1 eq lc $human_test || lc $new_allele2 eq lc $human_test){
				#print "reverse is a match!\n";	
			}else{
				print "Skipping rs $rs because reverse is NOT a match\n";
				print "$rs allele 1: $this_allele1\tallele 2: $this_allele2\thuman_test is $human_test\n";
				next;
			}
			
			
		}else{
			#this snp has an unusual format skipping
			print "Skipping rs $rs because SNP bases do not match\n";
			print "$rs allele 1: $this_allele1\tallele 2: $this_allele2\thuman_test is $human_test\n";
			next;
		}#check the REVERSE file 
	}
	#HANDLE THE GAP2 SNPS DIFFERENTLY!!!!
	
	$len = length($human);
	if($len>2){
		print "Skipping $rs because the multiz alignment is greater than 2\n";
		#Can't handle these
	}elsif($len == 2){
		#length is 2
		#is the first allele listed in alleles always the one we are looking for?!
		#print "length of $rs is exactly 2\n";
		#print "human is $human\n";
		#print Dumper \%out_hash;
		
		#now we need to do something different if it's a gap2 allele
		if(exists $all_gap2{$this_pos}){
			#print "gap 2 example\n";
			#this is a gap2 allele
			#take the second position!
			#print Dumper \%out_hash;
			foreach $f_spec (keys %out_hash){
				$bases=$out_hash{$f_spec};
				@bases=split('',$bases);
				$base=$bases[1];
				$new_out_hash{$f_spec}=$base;
			}
			#print "removing FIRST base!\n";
			%out_hash=%new_out_hash;
			#print Dumper \%new_out_hash;
			
			
		}else{
			#this is not a gap 2 allele, just take the first position	
			foreach $f_spec (keys %out_hash){
				$bases=$out_hash{$f_spec};
				@bases=split('',$bases);
				$base=$bases[0];
				$new_out_hash{$f_spec}=$base;
			}
			#print "removing second base!\n";
			%out_hash<-%new_out_hash;
			#print Dumper \%new_out_hash;
		}
		
	}elsif($len == 1){
		
		#CHECK TO SEE IF IT'S A GAP2 ALLELE
		if(exists $all_gap2{$rs}){
			#This means that humans are the ONLY species with this base
			foreach $f_spec (keys %out_hash){
				$new_out_hash{$f_spec}="-";
			}			
			#print "Human Specific Isert\n";
			$out_hash=$new_out_hash;
			#print Dumper \%new_out_hash;
		}else{
			#just one base! don't have to do anything	
		}
		
	}
	
	#Now print the out_hash!! one for each 
	#print "$rs out_hash\n";
	#print Dumper \%out_hash;
	if($RV_this == 1){
		$out_file1="$rs.$this_allele1.rev.fasta";	
	}else{
		$out_file1="$rs.$this_allele1.fasta";	
	}
	
	#ALLELE1
	
	$out_hash{"hg38"}=$this_allele1;
	@out_array1 = ();
	foreach $spec (keys %out_hash){
		$seq = $out_hash{$spec};
		$line = ">$spec\n$seq\n";
		push(@out_array1, $line);
	}
	
	$file_num=$rs_total/5000;
	$file_num=ceil($file_num);
	#@file_num=split('',$file_num);
	#$file_num=$file_num[0];	

	#open FILE , ">", "$out_file1";
	#print FILE @out_array1;
	#close FILE;
	
	open FILE, ">>", "$file_num.all_fasta.xmfa";
	print FILE "= $out_file1\n";
	print FILE @out_array1;
	close FILE;
	
	#CALCULATE PARSIMONY DATE
	
	#`Rscript parsimony_date_one.R $out_file1`;
	#`rm $out_file1`;

	
	if($RV_this == 1){
		$out_file2="$rs.$this_allele2.rev.fasta";	
	}else{
		$out_file2="$rs.$this_allele2.fasta";	
	}
	
	$out_hash{"hg38"}=$this_allele2;
	@out_array2 = ();
	foreach $spec (keys %out_hash){
		$seq = $out_hash{$spec};
		$line = ">$spec\n$seq\n";
		push(@out_array2, $line);
	}
	
	#open FILE , ">", "$out_file2";
	#print FILE @out_array2;
	#close FILE;
	
	#`Rscript parsimony_date_one.R $out_file2`;
	#`rm $out_file2`;
	

	
	open FILE, ">>", "$file_num.all_fasta.xmfa";
	print FILE "= $out_file2\n";
	print FILE @out_array2;
	close FILE;
	
	$rs_total=$rs_total+1;

}

#all the species in the dataset
