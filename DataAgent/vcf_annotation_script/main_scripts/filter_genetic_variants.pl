#!/usr/bin/perl
use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];

open(IN,"<$infile") or die "can not open $infile";
open(OUT,">$outfile") or die "can not open $outfile";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^Tag/){
		print OUT "$line\n";
	}else{
		my @recs=split(/\t/,$line);
		my $impact_func='pass';
		my $pop_tag='pass';
		my $func_tag='filter';
		if($recs[15] =~m/LOW/){
			$impact_func='filter';
		}
		for(my $i=25;$i<=30;$i++){
			if(($recs[$i] ne 'NA') and ($recs[$i] ne '.')){
				if($recs[$i] >= 0.05){
					$pop_tag='filter';
				}
			}
		}
		if(($recs[31] eq 'NA') or ($recs[31] eq '.')){#sift score
			$func_tag='pass';
		}else{
			if($recs[31]<=0.05){
				$func_tag='pass';
			}
		}
		if(($recs[33] eq 'NA') or ($recs[33] eq 'D') or ($recs[33] eq 'P') or ($recs[33] eq '.')){
			$func_tag='pass';
		}
		if(($recs[35] eq 'NA') or ($recs[35] eq 'D') or ($recs[35] eq 'P') or ($recs[35] eq '.')){
			$func_tag='pass';
		}
		if(($impact_func eq 'pass') and ($pop_tag eq 'pass') and ($func_tag eq 'pass')){
			print OUT "$line\n";
		}
	}
}
close IN;
close OUT;
