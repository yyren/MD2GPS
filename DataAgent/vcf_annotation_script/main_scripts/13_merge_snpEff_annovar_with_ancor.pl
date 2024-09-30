#!/usr/bin/perl -w
use strict;

my $in_snpEff_file=$ARGV[0];
my $in_annovar_file=$ARGV[1];
my $out_merged_file=$ARGV[2];

my (@header_names,$header_cosmic,$header_clinvar,$header_snp,$annovar_header,%annovar_record,@snpeff_header);
open(IN,"<$in_annovar_file") or die "can not open $in_annovar_file:$!";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//;
	if($_=~m/Start/){
		my @headers=split('\t',$_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header_names[$i]=$headers[$i];
			if($header_names[$i]=~m/cosmic/){
				$header_cosmic=$header_names[$i];
			}
			elsif($header_names[$i]=~m/clinvar/){
				$header_clinvar=$header_names[$i];
			}
			elsif($header_names[$i]=~m/snp/){
				$header_snp=$header_names[$i];
			}
		}
	}else{
		my @recs=split('\t',$_);
		my %hash_rec=();
		for(my $i=0;$i<scalar(@recs);$i++){
			$hash_rec{$header_names[$i]}=$recs[$i];
		}
		my $chr=$hash_rec{"Chr"};
		my $start=$hash_rec{"Start"};
		my $ref=$hash_rec{"Ref"};
		my $alt=$hash_rec{"Alt"};
		my $id=join("_",($chr,$start,$ref,$alt));
		my $gnomAD_all=$hash_rec{"gnomAD_exome_ALL"};# ALL
		my $gnomAD_eas=$hash_rec{"gnomAD_exome_EAS"};# east Asian
		my $gnomAD_sas=$hash_rec{"gnomAD_exome_SAS"};# Sourth Asian
		my $gnomAD_eur=$hash_rec{"gnomAD_exome_NFE"};#european without finnish
		my $gnomAD_amr=$hash_rec{"gnomAD_exome_AMR"};#Latino
		my $gnomAD_afr=$hash_rec{"gnomAD_exome_AFR"};#African
		my $sift_scrore=$hash_rec{"SIFT_score"};
		my $poly_hdiv_score=$hash_rec{"Polyphen2_HDIV_score"};
		my $poly_hdiv_pred=$hash_rec{"Polyphen2_HDIV_pred"};
		my $poly_hvar_score=$hash_rec{"Polyphen2_HVAR_score"};
		my $poly_hvar_pred=$hash_rec{"Polyphen2_HVAR_pred"};
		my $cadd_raw=$hash_rec{"CADD_raw"};
		my $cadd_phred=$hash_rec{"CADD_phred"};
		my $cosmic=$hash_rec{"$header_cosmic"};
		my $clinvar=$hash_rec{"$header_clinvar"};
		my $snp=$hash_rec{"$header_snp"};
		my $id_records=join("\t",($gnomAD_all,$gnomAD_eas,$gnomAD_sas,$gnomAD_eur,$gnomAD_amr,$gnomAD_afr,$sift_scrore,$poly_hdiv_score,$poly_hdiv_pred,$poly_hvar_score,$poly_hvar_pred,$cadd_raw,$cadd_phred,$cosmic,$clinvar,$snp));
		$annovar_header=join("\t",('gnomAD_exome_all','gnomAD_exome_eas','gnomAD_exome_sas','gnomAD_exome_eur','gnomAD_exome_amr','gnomAD_afr','SIFT_score','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','CADD_raw','CADD_phred',$header_cosmic,$header_clinvar,$header_snp));
		$annovar_record{$id}=$id_records;
	}
}
close IN;
open(IN,"<$in_snpEff_file") or die "can not open $in_snpEff_file:$!";
open(OUT,">$out_merged_file") or die "can not open $out_merged_file:$!";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//;
	if($_=~m/CHROM/){
		my @recs=split('\t',$_);
		printf OUT "$_\t$annovar_header\n";
		for(my $i=0;$i<scalar(@recs);$i++){
			$snpeff_header[$i]=$recs[$i];
		}
	}else{
		my @recs=split('\t',$_);
		my %snpeff_rec=();
		for(my $i=0;$i<scalar(@recs);$i++){
			$snpeff_rec{$snpeff_header[$i]}=$recs[$i];
		}
		my $id=join("_",($snpeff_rec{"CHROM"},$snpeff_rec{"ID"},$snpeff_rec{"REF"},$snpeff_rec{"ALT"}));
		my $merged_rec=join("\t",($_,$annovar_record{$id}));
		printf OUT "$merged_rec\n";
	}
}
close IN;
close OUT;
