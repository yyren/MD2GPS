#!/usr/bin/perl
use strict;

my $in_file=$ARGV[0];
open(IN,"<$in_file") or die "can not open $in_file:$!";
while(<IN>){
	chomp $_;
	my $format_idx=0;
	if($_=~m/#CHROM/){
		my @recs=split('\t',$_);
		for(my $i=0;$i<scalar(@recs);$i++){
			if($recs[$i] eq 'FORMAT'){
				$format_idx=$i;
			}
		}
		my $sample_numb=$#recs-$format_idx;
		if($sample_numb==1){
			print "one";
		}elsif($sample_numb==2){
			print "two";
		}else{
			print "multi";
		}
		last;
	}
}
close IN;