use strict;
use warnings;

my $in_vcf=$ARGV[0];
my $in_hg38_bed=$ARGV[1];
my $out_hg19_vcf=$ARGV[2];
my $out_hg38_vcf=$ARGV[3];

my %hg19_to_hg38=();
open(IN,"<$in_hg38_bed") or die "cannot open $in_hg38_bed\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/^chr/){
        $line=~s/^chr//;
    }
    my @recs=split(/\t/,$line);
    my $hg19_pos=join("_",@recs[0,3]);
    $hg19_to_hg38{$hg19_pos}=$recs[2];
}
close IN;

open(IN,"<$in_vcf") or die "cannot open $in_vcf\n";
open(OUT1,">$out_hg19_vcf") or die "cannot open $out_hg19_vcf\n";
open(OUT2, ">$out_hg38_vcf") or die "cannot open $out_hg38_vcf\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/^#/){
        print OUT1 "$line\n";
        print OUT2 "$line\n";
    }else{
        my @recs=split(/\t/,$line);
        my $chr=$recs[0];
        if($chr=~m/^chr/){
            $chr=~s/chr//;
        }
        my $pos=join("_",($chr,$recs[1]));
        if(not exists $hg19_to_hg38{$pos}){
            print OUT1 "$line\n";
        }else{
            $recs[2]=$recs[1];
            $recs[1]=$hg19_to_hg38{$pos};
            my $new_line=join("\t",@recs);
            print OUT2 "$new_line\n";
        }
    }
}
close IN;
close OUT1;
close OUT2;