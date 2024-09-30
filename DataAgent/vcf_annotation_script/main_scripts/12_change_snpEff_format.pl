#!/usr/bin/perl -w
use strict;

my $in_vcf=$ARGV[0];
my $out_vcf=$ARGV[1];
my $single_or_compare=$ARGV[2];
my $rec_basic;
my $eff;
my @effs=();
my @annotation=();
my $annotation;
my $Effect;
my $Effect_Impact;
my $Functional_Class;
my $Codon_Change;
my $Amino_Acid_Change;
my $Amino_Acid_Length;
my $Gene_Name;
my $Transcript_BioType;
my $Gene_Coding;
my $Transcript_ID;
my $Exon_Rank;
my $Genotype_Number;
my @alts;
my @covs;
my @var_covs;
my $tumor_freq;
my ($normal_freq,$tumor_freq_temp);


open(IN, "<", $in_vcf) or die "can not open < $in_vcf: $!";
open(OUT, ">", $out_vcf) or die "cannot open < $out_vcf: $!";
printf OUT "Source	VarFreq_precent_max	Cov_max	Var_Cov_max	CHROM	ID	REF	ALT	VarFreq_precent	Cov	Var_Cov	Genome_Type	Effect	Effect_Impact	Functional_Class	Codon_Change	Amino_Change	Amino_Acid_Length	Gene_Name	Transcript_BioType	Gene_Coding	Transcript_ID	Exon_Rank\n";
while (<IN>) {
    chomp $_;
    $_=~s/\r|\n//g;
    @_=split('\t', $_);
    if (!($_[0] =~m "#")) {
        
        $eff=$_[7];
        @effs=split('\;EFF=', $eff);
		if($_[4]=~m/\,/g){
			@alts=split("\,", $_[4]);
			if($alts[0] eq "$_[3]"){
				$_[4]=$alts[1];
			}
			else{
				$_[4]=$alts[0];
			}
		}
		if($_[6]=~m/\_/g){
			@covs=split('_', $_[6]);
			@var_covs=split('_', $effs[0]);
			if($single_or_compare eq 'multi'){
				for(my $i=0;$i<scalar(@covs);$i++){
					if($covs[$i]==0){
						$tumor_freq=0.00;
					}elsif(($covs[$i]==-1)||($var_covs[$i]==-1)){
						$tumor_freq=-1;
					}else{
						$tumor_freq=sprintf("%0.2f",($var_covs[$i]/$covs[$i]*100));
					}
					if($i==0){
						$tumor_freq_temp=$tumor_freq;
					}else{
						$tumor_freq=join("_",($tumor_freq_temp,$tumor_freq));
						$tumor_freq_temp=$tumor_freq;
					}
				}
			}elsif($single_or_compare eq 'compare'){
				if($covs[0]==0){
					$tumor_freq=0;
				}elsif(($covs[0]==-1)||($var_covs[0]==-1)){
					$tumor_freq=-1;
				}else{
					$tumor_freq=sprintf("%0.2f",($var_covs[0]/$covs[0]*100));
				}
			}
			$_[5]=$tumor_freq;
		}
		$rec_basic=join("\t", ($_[9],$_[5],$_[6],$effs[0],$_[0],$_[1],@_[3..6],$effs[0],$_[8]));
        $eff=$effs[1];
        if($eff=~m/\;/g){
			@effs=split("\;", $eff);
			$eff=$effs[0];
		}
		if($eff=~m/\,/g){
			@effs=split('\,', $eff);
			for(my $i=0; $i < scalar(@effs); $i++){
				@annotation=split("\\(", $effs[$i]);
				$Effect=$annotation[0];
				$annotation=$annotation[1];
				$annotation=~s/\)//g;
				($Effect_Impact, $Functional_Class, $Codon_Change, $Amino_Acid_Change, $Amino_Acid_Length, $Gene_Name, $Transcript_BioType, $Gene_Coding, $Transcript_ID, $Exon_Rank, $Genotype_Number)=split('\|', $annotation);
				#if ((!($Amino_Acid_Change eq "c\." ))&&(!($Effect=~m/UTR/))&&(!($Effect=~m/TF_binding_site/))&&(!($Effect eq "downstream_gene_variant"))&&(!($Effect eq "intron_variant"))&&(!($Effect eq "intragenic_variant"))&&(!($Effect eq "intergenic_region"))&&(!($Effect eq "upstream_gene_variant"))&&(!($Effect eq "downstream_gene_variant"))&&(!($Effect eq "synonymous_variant")))  {
					if($Effect ne "transcript"){
						$annotation=join("\t", ($Effect_Impact, $Functional_Class, $Codon_Change, $Amino_Acid_Change, $Amino_Acid_Length, $Gene_Name, $Transcript_BioType, $Gene_Coding, $Transcript_ID, $Exon_Rank));
						printf OUT "$rec_basic\t"."$Effect\t$annotation\n";
					}
				#}
            }
		}
		else{
			@annotation=split("\\(", $eff);
			$Effect=$annotation[0];
			$annotation=$annotation[1];
			$annotation=~s/\)//g;
			($Effect_Impact, $Functional_Class, $Codon_Change, $Amino_Acid_Change, $Amino_Acid_Length, $Gene_Name, $Transcript_BioType, $Gene_Coding, $Transcript_ID, $Exon_Rank, $Genotype_Number)=split('\|', $annotation);
			#if ((!($Amino_Acid_Change eq "c\." ))&&(!($Effect=~m/UTR/))&&(!($Effect=~m/TF_binding_site/))&&(!($Effect eq "downstream_gene_variant"))&&(!($Effect eq "intron_variant"))&&(!($Effect eq "intragenic_variant"))&&(!($Effect eq "intergenic_region"))&&(!($Effect eq "upstream_gene_variant"))&&(!($Effect eq "downstream_gene_variant"))&&(!($Effect eq "synonymous_variant")))  {
				if($Effect ne "transcript"){
					$annotation=join("\t", ($Effect_Impact, $Functional_Class, $Codon_Change, $Amino_Acid_Change, $Amino_Acid_Length, $Gene_Name, $Transcript_BioType, $Gene_Coding, $Transcript_ID, $Exon_Rank));
					printf OUT "$rec_basic\t"."$Effect\t$annotation\n";
				}
			#}
        }
    }
}
close IN;
close OUT;
