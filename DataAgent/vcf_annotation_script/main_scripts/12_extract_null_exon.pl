#!/usr/bin/perl -w
use strict;
#extrcat exon_rank eq null
my $in_vcf_file=$ARGV[0];
my $out_exon_null_file=$ARGV[1];
my $out_snpeff_input_file=$ARGV[2];

my @headers;
my @recs;
my $header;
my %header_rec=();
my $exon;
my $header_chr;
my $chr;
my $start;
my $ref;
my $alt;
my $varfreq;
my $cov;
my $var_cov;
open(IN,"<$in_vcf_file") or die"can not open $in_vcf_file\n";
open(OUT,">$out_exon_null_file") or die"can not open $out_exon_null_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~/CHROM/){
		@headers=split('\t', $_);
		printf OUT "$_\n";
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
#		$exon=$header_rec{"Exon_Rank"};
#print "$exon\n";
		if($header_rec{"Exon_Rank"} eq ""){
			#print "$_\n";
			printf OUT "$_\n";
		}
	}
}
close IN;
close OUT;

open(IN,"<$out_exon_null_file") or die"can not open $out_exon_null_file\n";
open(OUT,">$out_snpeff_input_file") or die"can not open $out_snpeff_input_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~/CHROM/){
		printf OUT "#CHROM	POS	ID	REF	ALT	VarFreq_precent	Cov	Var_Cov\n";
		@headers=split('\t', $_);
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
			if($header=~/CHROM/){
				$header_chr=$header;
			}
		}
		$chr=$header_rec{"$header_chr"};
		$start=$header_rec{"ID"};
		$varfreq=$header_rec{"VarFreq_precent"};
		$cov=$header_rec{"Cov"};
		$var_cov=$header_rec{"Var_Cov"};
		$ref=$header_rec{"REF"};
		$ref=substr($ref,0,1);
		if($ref eq 'A'){
			$alt='T';
		}
		elsif($ref eq 'T'){
			$alt='A';
		}
		elsif($ref eq 'C'){
			$alt='G';
		}
		elsif($ref eq 'G'){
			$alt='C';
		}
		printf OUT "$chr\t$start\t\.\t$ref\t$alt\t$varfreq\t$cov\t$var_cov\n";
	}
}
close IN;
close OUT;
