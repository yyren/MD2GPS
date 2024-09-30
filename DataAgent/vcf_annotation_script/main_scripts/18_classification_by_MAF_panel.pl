#!/usr/bin/perl -w
use strict;
my $in_vcf=$ARGV[0];
my $out_MAF_less_1_precent=$ARGV[1];
my $race=$ARGV[2];

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
my @covs;
my @varcovs;
my $tumor_varfreq;
my $normal_varfreq;
my $single_or_compare;

open(IN, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT, ">", $out_MAF_less_1_precent) or die "cannot open $out_MAF_less_1_precent: $!";

my $allel;
my $ref;

while (<IN>) {
    chomp $_;
    $_=~s/\r|\n|\#//g;
    @rec_numb=split("\t", $_);
    if (scalar(@rec_numb) > 15){
		if($_=~m/CHROM/){
			printf OUT "Tag\tCOSMIC\t$_\n";
			@headers=split("\t", $_);
			for(my $i=0;$i<scalar(@headers);$i++){
				$header=$headers[$i];
				if($header=~m/cosmic/){
					$header_cosmic=$header;
				}
				# if($header=~m/1000g/i){
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
			$normal_varfreq=0;
			$tumor_varfreq=$header_rec{"VarFreq_precent_max"};
			$single_or_compare="single";
			if($header_rec{"Cov_max"}=~m/\_/g){
				$single_or_compare="compare";
				@covs=split('_',$header_rec{"Cov_max"});
				@varcovs=split('_',$header_rec{"Var_Cov_max"});
				$tumor_varfreq=$varcovs[0]/$covs[0];
				if($covs[1]==0){
					$normal_varfreq=0;
				}else{
					$normal_varfreq=$varcovs[1]/$covs[1];
				}
			}
			if(($tumor_varfreq!=-1)&&($tumor_varfreq < $normal_varfreq)){
				if($cosmic_id=~/\;/g){
					@cosmics=split("\;", $cosmic_id);
					$cosmic_id=$cosmics[0];
				}
				$cosmic_id=~s/ID=//g;
				$cosmic_id=~s/COSM//g;
				$cosmic_id=~s/,/;/g;
				printf OUT "abnormal\t$cosmic_id\t$_\n";
			}else{
				if($cosmic_id eq "NA" ){
					if ($one_k eq "NA") {
						printf OUT "target\tNA\t$_\n";
					 }
					else{
						if ($one_k <= 0.01) {
							printf OUT "target\tNA\t$_\n";
						}
						else{
							if($single_or_compare eq "single"){
								printf OUT "chemotherapy\tNA\t$_\n";
							}else{
								printf OUT "common-somatic\tNA\t$_\n";
							}
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
					if ($one_k eq "NA") {
						printf OUT "target\t$cosmic_id\t$_\n";
					}
					else{
						if ($one_k <= 0.01) {
							printf OUT "target\t$cosmic_id\t$_\n";
						}else{
							if($single_or_compare eq "single"){
								printf OUT "chemotherapy\t$cosmic_id\t$_\n";
							}else{
								printf OUT "common-somatic\t$cosmic_id\t$_\n";
							}
						}
					}
				}
			}	
        } 
    }
}
close IN;
close OUT;

