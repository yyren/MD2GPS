#!/usr/bin/perl -w
use strict;
my $in_target=$ARGV[0];
my $in_vcf=$ARGV[1];
my $out_vcf=$ARGV[2];
my $race=$ARGV[3];

my $population_tag='gnomAD_exome_'.$race;
my @headers;
my @recs;
my $header;
my %header_rec=();
my $chr;
my $loci;
my $ref;
my $allel;
my $cosmic_id;
my $dna_change;
my $gene;
my $dna_change_ori;
my $mutation;
my %mutations=();
my %mutations_ori=();
my @alleles;
my $allel_frac;
my $amino_change_idx;
my $cosmic;
my @cosmics;
my @dna_changes;
my $one_k_header;
my $one_k;
my @covs;
my @varcovs;
my $tumor_varfreq;
my $normal_varfreq;
my $single_or_compare;

open(IN_target, "<", $in_target) or die "can not open < $in_target: $!";
open(IN_vcf, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT_target, ">", $out_vcf) or die "cannot open > $out_vcf: $!";

while (<IN_target>) {
    chomp $_;
    $_=~s/\r|\n//g;
    $_=~s/ //g;
	if($_=~m/Gene/g){
		@headers=split("\t", $_);
	}
	else{
		@recs=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"chr"};
		$chr=~s/chr//g;
		$loci=$header_rec{"start"};
		$cosmic_id=$header_rec{"COSMIC_id"};
		$dna_change=$header_rec{"CDS_mut_syntax"};
		$gene=$header_rec{"Gene"};
		if ($dna_change=~m/del/) {
			@_=split('del', $dna_change);
			$dna_change=$_[0];
			#print "$dna_change\n";
		}else{
			if ($dna_change=~m/_/g) {
				$dna_change_ori=$dna_change;
				$dna_change=~s/_\d+//g;
				$mutation=join("\t", ($chr, $gene, $dna_change, $loci));
				$mutations_ori{$mutation}=$dna_change_ori;
			}
		}
		#$mutation=join("\t", ($chr, $loci, $dna_change));
		$mutation=join("\t", ($chr, $gene, $dna_change, $loci));
		$mutations{$mutation}=$cosmic_id;
		#print "$mutation\n";
	}
}
close IN_target;


while (<IN_vcf>) {
    chomp $_;
    $_=~s/\r|\n//g;
    if ($_=~m/CHROM/) {
        $_=~s/#//g;
        printf OUT_target "Tag\tCOSMIC\t$_\n";
		@headers=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			if($header=~m/Amino\_Change/g){
				$amino_change_idx=$i;
			}
			if($header=~m/cosmic/g){
				$cosmic=$header;
			}
			# if($header=~m/1000g/g){
				# $one_k_header=$header;
			# }
			if($header eq $population_tag){
				$one_k_header=$header;
			}
		}
    }
    else{
		@recs=split("\t", $_);
		if (scalar(@recs) > 15){
			for(my $i=0;$i<scalar(@recs);$i++){
				$header=$headers[$i];
				$header_rec{$header}=$recs[$i];
			}
			
			$allel_frac=$header_rec{"VarFreq_precent_max"};
			$chr=$header_rec{"CHROM"};
			if($chr=~m/chr/g){
				$chr=~s/chr//g;
			}
			$loci=$header_rec{"ID"};
			$ref=$header_rec{"REF"};
			$allel=$header_rec{"ALT"};
			$dna_change=$header_rec{"Amino_Change"};
			$gene=$header_rec{"Gene_Name"};
			$one_k=$header_rec{$one_k_header};
			$cosmic_id=$header_rec{$cosmic};
			#print "$loci\n$dna_change\n";
			if ($dna_change=~m/p\./) {
				@dna_changes=split("\/", $dna_change);
				$dna_change=$dna_changes[1];
			}
			if ($dna_change=~m/del/) {
				@dna_changes=split('del', $dna_change);
				$dna_change=$dna_changes[0];
			}else{
				if ($dna_change=~m/_/g) {
					$dna_change=~s/_\d+//g;
				}
			}

			#print "$loci\t$dna_change\n";
			#$mutation=join("\t", ($chr, $loci, $dna_change));
			$normal_varfreq=0;
			$tumor_varfreq=$header_rec{"VarFreq_precent_max"};
			$mutation=join("\t", ($chr, $gene, $dna_change, $loci));
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
				printf OUT_target "abnormal\t$cosmic_id\t$_\n";
			}elsif ((exists $mutations{$mutation}) &&(($one_k eq "NA") ||($one_k <=0.01))) {
				$cosmic_id=$mutations{$mutation};
				if(exists $mutations_ori{$mutation}) {
					$recs[$amino_change_idx]=~s/$dna_change/$mutations_ori{$mutation}/g;
					$_=join("\t", @recs);
					printf OUT_target "target\t$cosmic_id\t$_\n";
				}else{
					printf OUT_target "target\t$cosmic_id\t$_\n";
				}           
			}
			else{
				if ($cosmic_id eq "NA" ) {
					#code
				}
				else{
					if($cosmic_id=~m/\;/g){
						@cosmics=split("\;", $cosmic_id);
						$cosmic_id=$cosmics[0];
					}
					$cosmic_id=~s/ID=//g;
					$cosmic_id=~s/COSM//g;
					$cosmic_id=~s/,/;/g;
				}
				if(($one_k eq "NA") ||($one_k <=0.01)){
					printf OUT_target "non-target\t$cosmic_id\t$_\n";
				}else{
					if($single_or_compare eq "single"){
						printf OUT_target "chemotherapy\t$cosmic_id\t$_\n";
					}else{
						printf OUT_target "common-somatic\t$cosmic_id\t$_\n";
					}
				}
			}
		}
	} 
}
close IN_vcf;
close OUT_target;
