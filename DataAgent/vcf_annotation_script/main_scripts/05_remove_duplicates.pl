#!/usr/bin/perl
use strict;

my $in_file=$ARGV[0];
my $out_file=$ARGV[1];

my %rec_hash;
open(IN,"<$in_file") or die "can not open $in_file:$!";
open(OUT,">$out_file") or die "can not open $out_file:$!";
while(<IN>){
	chomp $_;
	if($_=~m/CHROM/){
		printf OUT "$_\n";
	}else{
		my @recs=split('\t',$_);
		my $id=join('_',@recs[0,1,3,4]);
		$rec_hash{$id}=$_;
	}
}
close IN;
foreach my $id(keys %rec_hash){
	printf OUT "$rec_hash{$id}\n";
}
close OU;