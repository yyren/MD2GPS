#!/usr/bin/perl -w
use strict;

my $in_target_region_bed=$ARGV[0];
my $in_clean_snpeff=$ARGV[1];
my $out_target_variants=$ARGV[2];

my @recs;
my $chr;
my $start;
my $end;
my $loci;
my $id;
my %id_target=();
my $header;
my @headers;
my %header_rec;

open(IN,"<$in_target_region_bed") or die "can not open $in_target_region_bed\n";

while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	@recs=split('\t', $_);
	$chr=$recs[0];
	#$chr=~s/chr//g;
	$start=$recs[1];
	$end=$recs[2];
	my $i;
	for($i=$start;$i<=$end;$i++){
		$loci=$i;
		$id=join("_", ($chr,$loci));
		$id_target{$id}=$_;
	}
}

close IN;

open(IN,"<$in_clean_snpeff") or die "can not open $in_clean_snpeff\n";
open(OUT,">$out_target_variants") or die "can not open $out_target_variants\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/){
		#printf OUT "$_\n";
		@headers=split('\t', $_);
		
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@recs);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"CHROM"};
		#$chr=~s/chr//g;
		$start=$header_rec{"ID"};
		$id=join("_", ($chr,$start));
		if(exists $id_target{$id}){
			printf OUT "$_\n";
		}
	}

}
close IN;
close OUT;
