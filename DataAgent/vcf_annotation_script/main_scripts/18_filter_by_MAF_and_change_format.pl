#!/usr/bin/perl -w
use strict;
my $in_vcf=$ARGV[0];
my $out_MAF_less_1_precent=$ARGV[1];
my $out_MAF_more_than_1_precent=$ARGV[2];
my $out_all_vcf=$ARGV[3];
my $race=$ARGV[4];

my $population_tag='gnomAD_exome_'.$race;
my @rec_numb;
my @headers;
my $header;
my $header_cosmic;
my $header_1000g;
my @recs;
my %header_rec=();
my $chr;
my $id;
my $mutation;
my %mutations=();
my $gene;
my $dna_change;
my $p_change;
my $cosmic_id;
my $loci;
my $one_k;
my @cosmics;

open(IN, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT, ">", $out_MAF_less_1_precent) or die "cannot open < $out_MAF_less_1_precent: $!";
open(OUT_MAF_more_than_1_precent, ">", $out_MAF_more_than_1_precent) or die "cannot open < $out_MAF_more_than_1_precent: $!";
open(OUT_all,">", $out_all_vcf) or die "cannot open < $out_all_vcf: $!";
my $allel;
my $ref;

while (<IN>) {
    chomp $_;
    $_=~s/\r|\n|\#//g;
    @rec_numb=split("\t", $_);
    if (scalar(@rec_numb) > 15){
		if($_=~m/CHROM/g){
			printf OUT "COSMIC\t$_\n";
            printf OUT_MAF_more_than_1_precent "COSMIC\t$_\n";
			printf OUT_all "COSMIC\t$_\n";
			@headers=split("\t", $_);
			for(my $i=0;$i<scalar(@headers);$i++){
				$header=$headers[$i];
				if($header=~m/cosmic/g){
					$header_cosmic=$header;
				}
				# if($header=~m/1000g/g){
					# $header_1000g=$header;
				# }
				if($header eq $population_tag){
					$header_1000g=$header;
				}
			}
			
		}
		else{
			@recs=split("\t", $_);
			for(my $i=0;$i<scalar(@headers);$i++){
				$header=$headers[$i];
				$header_rec{$header}=$recs[$i];
			}
			$chr=$header_rec{"CHROM"};
			$loci=$header_rec{"ID"};
			$ref=$header_rec{"REF"};
			$allel=$header_rec{"ALT"};
			$one_k=$header_rec{$header_1000g};
			$cosmic_id=$header_rec{$header_cosmic};
					
			if($cosmic_id eq "NA" ){
				printf OUT_all "NA\t$_\n";
				if ($one_k eq "NA") {
                    printf OUT "NA\t$_\n";
               	 }
				else{
                    if ($one_k <= 0.01) {
						printf OUT "NA\t$_\n";
                  	}
					else{
                      printf OUT_MAF_more_than_1_precent "NA\t$_\n";
					}
				}
			}
			else{
				if($cosmic_id=~/\;/g){
					@cosmics=split("\;", $cosmic_id);
					$cosmic_id=$cosmics[0];
				}
                $cosmic_id=~s/ID=//g;
                $cosmic_id=~s/COSM//g;
				$cosmic_id=~s/,/;/g;
				printf OUT_all "$cosmic_id\t$_\n";
                if ($one_k eq "NA") {
					printf OUT "$cosmic_id\t$_\n";
                }
				else{
                    if ($one_k <= 0.01) {
						printf OUT "$cosmic_id\t$_\n";
                    }else{
                        printf OUT_MAF_more_than_1_precent "$cosmic_id\t$_\n";
                    }
                }
            }  
        } 
    }
}
close IN;
close OUT;
close OUT_MAF_more_than_1_precent;
close OUT_all;
