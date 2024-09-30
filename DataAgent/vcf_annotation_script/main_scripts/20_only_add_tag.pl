#!/usr/bin/perl
use strict;

my $in_file=$ARGV[0];
my $out_file=$ARGV[1];
my (@recs,@headers,$header,$header_cosmic,$cosmic_id,@cosmics,%header_rec);
open(IN, "<$in_file") or die "can not open $in_file\n";
open(OUT, ">$out_file") or die "can not open $out_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		@headers=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			if($header=~m/cosmic/g){
				$header_cosmic=$header;
			}
		}
		printf OUT "Tag\tCOSMIC\t$_\n";
	}else{
		@recs=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$cosmic_id=$header_rec{$header_cosmic};
		if($cosmic_id eq "NA" ){
			printf OUT "target\tNA\t$_\n";
		}
		else{
			if($cosmic_id=~/\;/g){
				@cosmics=split("\;", $cosmic_id);
				$cosmic_id=$cosmics[0];
			}
            $cosmic_id=~s/ID=//g;
            $cosmic_id=~s/COSM//g;
			$cosmic_id=~s/,/;/g;
			printf OUT "target\t$cosmic_id\t$_\n";
		}
	}

}