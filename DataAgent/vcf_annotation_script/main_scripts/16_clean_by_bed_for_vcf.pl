#!/usr/bin/perl -w
use strict;
my $in_target=$ARGV[0];
my $in_vcf=$ARGV[1];
my $out_vcf=$ARGV[2];
my $chr;
my $loci;
my $id;
my %id_target=();
my $start;
my $end;
my @headers;
my $header;
my %header_rec=();
my @recs;
open(IN_target, "<", $in_target) or die "can not open < $in_target: $!";
open(IN, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT, ">", $out_vcf) or die "cannot open < $out_vcf: $!";
while (<IN_target>) {
    chomp $_;
    $_=~s/\r|\n//g;
    @_=split('\t', $_);
    if ($_[0] eq "TargetA") {
        #code
    }else{
        $chr=$_[0];
        $chr=~s/chr//g;
        $start=$_[1];
        $end=$_[2];
        my $i=0;
        for($i=$start; $i<$end; $i++){
            $loci=$i;
            $id=join("\_", ($chr, $loci));
            if (exists $id_target{$id}) {
                $id_target{$id}=$id_target{$id}+1;
            }
            else{
                $id_target{$id}=1;
            }
			#print "1:$id_target{$id}\n";
        }
    }    
}
close IN_target;

while (<IN>) {
    chomp $_;
    $_=~s/\r|\n//;
	if($_=~/^##/){
		print OUT "$_\n";
	}elsif($_=~m/^#CHROM/){
		printf OUT "$_\n";
		$_=~s/#//;
		@headers=split("\t", $_);
	}
	else{
		@recs=split("\t", $_);
		for(my $i=0; $i<scalar(@headers); $i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"CHROM"};
		$chr=~s/chr//;
		if(exists $header_rec{"POS"}){
			$start=$header_rec{"POS"};
		}elsif(exists $header_rec{"ID"}){
			$start=$header_rec{"ID"};
		}
		$id=join("\_", ($chr, $start));
		#print "2:$id\n";
		if (exists $id_target{$id}) {
			printf OUT "$_\n";
		}
	}
}
close IN;
close OUT;
