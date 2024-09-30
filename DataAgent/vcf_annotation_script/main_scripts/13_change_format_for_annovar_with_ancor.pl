#!/usr/bin/perl -w
use strict;
my $in_file=$ARGV[0];
my $out_file=$ARGV[1];

my @header;
open(IN,"<$in_file") or die "can not open $in_file:$!";
open(OUT,">$out_file") or die "can not open $out_file:$!";
while(<IN>){
	chomp $_;
	my @recs=split('\t',$_);
	if($_=~m/CHROM/){
		for(my $i=0;$i<scalar(@recs);$i++){
			$header[$i]=$recs[$i];
		}
	}else{
		my %header_rec=();
		for(my $i=0;$i<scalar(@recs);$i++){
			$header_rec{$header[$i]}=$recs[$i];
		}
		my $chr=$header_rec{"CHROM"};
		my $start=$header_rec{"ID"};
		my $ref=$header_rec{"REF"};
		my $alt=$header_rec{"ALT"};
		my $end=length($ref)+$start-1;
		my $id_rec=join("\t",($chr,$start,$end,$ref,$alt));
		printf OUT "$id_rec\n";
	}
}
close IN;
close OUT;