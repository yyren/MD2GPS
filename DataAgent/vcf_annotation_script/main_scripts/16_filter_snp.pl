#!/usr/bin/perl -w
use strict;
my $in_vcf=$ARGV[0];
my $out_vcf=$ARGV[1];

my @headers;
my @recs;
my $header;
my $snp;
my %header_rec=();


open(IN, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT, ">", $out_vcf) or die "cannot open < $out_vcf: $!";
while (<IN>) {
    chomp $_;
    $_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		$_=~s/#//g;
		printf OUT "$_\n";
		@headers=split("\t", $_);
	}
	else{
		@recs=split("\t", $_);
		for(my $i=0; $i<scalar(@headers); $i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$snp=$header_rec{"avsnp142"};
		if($snp=~m/rs/){
			#it is common snp
		}else{
			printf OUT "$_\n";
		}
	}
}
close IN;
close OUT;
