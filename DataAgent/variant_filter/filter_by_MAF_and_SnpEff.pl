use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];

open(IN,"<$infile") or die "cannot open $infile";
open(OUT,">$outfile") or die "cannot open $outfile";

while(my $line=<IN>){
	chomp $line;
	if($line=~m/^Tag/){
		print OUT "$line\n";
	}else{
		my @recs=split(/\t/,$line);
		
		
		my $impact_func='filter';
		my $pop_tag='pass';
		
		if($recs[15] !~m/LOW/){
			$impact_func='pass';
		}
		for(my $i=25;$i<=30;$i++){
			if(($recs[$i] ne 'NA') and ($recs[$i] ne '.')){
				if($recs[$i] >= 0.01){
					$pop_tag='filter';
				}
			}
		}
		
		if(($impact_func eq 'pass') and ($pop_tag eq 'pass')){
			print OUT "$line\n";
		}
	}
}
close IN;
close OUT;
