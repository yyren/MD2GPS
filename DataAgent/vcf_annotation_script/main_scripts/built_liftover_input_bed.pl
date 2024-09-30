use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];

open(IN,"<$infile") or die "can not open $infile\n";
open(OUT,">$outfile") or die "can not open $outfile\n";
while(my $line=<IN>){
    $line=~s/\r|\n//g;
    if($line=~m/CHROM/){
        #print OUT "$line\n";
    }else{
        my @recs=split(/\t/,$line);
        if($recs[0]!~m/chr/){
            $recs[0]='chr'.$recs[0];
        }
		my $start=$recs[1]-1;
        my $new_line=join("\t",($recs[0],$start,$recs[1],$recs[1]));
        print OUT "$new_line\n";
    }
}
close IN;
close OUT;
