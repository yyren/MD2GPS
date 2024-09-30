#!/usr/bin/perl -w
use strict;
die "USAGE:\n\tperl $0 <IN_file><OUT_file><indel_length>\n" if(@ARGV!=3);
##e.g perl  nn_indel_code130.pl ./results_meaningful_ci_and_p_delin.txt ./out.filter.txt 130


my $in_file=$ARGV[0];
my $out_file=$ARGV[1];
my $indel_length=$ARGV[2];

my @headers;
my $header;
my %header_idx;
my @recs;
my $ref;
my $ref_length;
my $tumor_vf;
my $amino_Change;
my $tumor_cov;
my $tumor_var_cov;
my $del_seq;
my $del_length;
my $line;
my $ins_seq;
open(IN,"<$in_file") or die "can not open $in_file\n";
open(OUT1,">$out_file") or die "can not open $out_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		printf OUT1 "$_\n";
		my @headers=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header_idx{$headers[$i]}=$i;
		}
	}else{
		@recs=split("\t", $_);
		$ref=$recs[$header_idx{"REF"}];
		$ref_length=length($ref);
		$amino_Change=$recs[$header_idx{"Amino_Change"}];
		$tumor_vf=$recs[$header_idx{"VarFreq_precent_max"}];
		$tumor_cov=$recs[$header_idx{"Cov_max"}];
		$tumor_var_cov=$recs[$header_idx{"Var_Cov_max"}];
		#if(($tumor_vf<=100)&&($tumor_vf>=5)){#filter variants with variant frequency<=5 (date:2019.2.22)
		if($ref_length >= 100){
			$ref=$ref_length;
		}
		####if the length of insertion is bigger than $indel_length, remove it
		if(($amino_Change=~m/ins([ATCG]+)$/g) && (length($1) > $indel_length)){
		
		}elsif(($amino_Change=~m/dup([ATCG]+)$/g) && (length($1) > $indel_length)){  
			#delete this site
		}elsif(($amino_Change=~/del([ATCG]+)/) && (length($1) > $indel_length)){
			$del_seq=$1;
			$del_length=length($del_seq);
			if($del_length >= 100){
				$amino_Change=~s/$del_seq/$del_length/;
			}
			$line=join("\t",(@recs[0..($header_idx{"REF"}-1)],$ref,@recs[($header_idx{"REF"}+1)..($header_idx{"Amino_Change"}-1)],$amino_Change,@recs[($header_idx{"Amino_Change"}+1)..$#recs]));
			print OUT1 "$line\n";
		}else{
			$line=join("\t",(@recs[0..($header_idx{"REF"}-1)],$ref,@recs[($header_idx{"REF"}+1)..$#recs]));
			print OUT1 "$line\n";
		}
		#} 
	}
}


close IN;
close OUT1;
