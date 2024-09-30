#!/usr/bin/perl -w
use strict;
my $in_rank_null_vcf_file=$ARGV[0];
my $in_snpeff_annotated_vcf=$ARGV[1];
my $in_rank_null_all_file=$ARGV[2];
my $out_combain_vcf=$ARGV[3];

my @headers;
my @recs;
my $header;
my %header_rec=();
my $exon;
my $chr;
my $start;
my $id;
my %exon_num=();
my $ref;
my $alt;
my $id2;
my $id3;
my $col_num;
my $no_exon_rec;
my %new_recs=();
my $id4;
my $eff;
my $eff_idx;
my @effs;
my @annotations;
my $annotation;
my $transcript_id;

open(IN1,"<$in_rank_null_vcf_file") or die"can not open $in_rank_null_vcf_file\n";
open(IN2,"<$in_snpeff_annotated_vcf") or die"can not open $in_snpeff_annotated_vcf\n";
open(IN3,"<$in_rank_null_all_file") or die"can not open $in_rank_null_all_file\n";
open(OUT,">$out_combain_vcf") or die"can not open $out_combain_vcf\n";
while(<IN2>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/##/g){
		#nothing
	}
	elsif($_=~m/CHROM/g){
		$_=~s/#//g;
		@headers=split('\t', $_);
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
			if($recs[$i]=~m/EFF\=/g){
				$eff_idx=$i;
			}
		}
		$chr=$header_rec{"CHROM"};
		$start=$header_rec{"POS"};
		$eff=$recs[$eff_idx];
		@effs=split("EFF=", $eff);
		$eff=$effs[1];
		if($eff=~m/\;/g){
			@effs=split('\;', $eff);
			$eff=$effs[0];
			if($eff=~m/\,/g){
				@effs=split('\,', $eff);
				for(my $i=0; $i < scalar(@effs); $i++){
					@annotations=split('\(', $effs[$i]);
					$annotation=$annotations[1];
					$annotation=~s/\)//g;
					@annotations=split('\|', $annotation);
					$transcript_id=$annotations[8];
					$exon=$annotations[9];
					$id=join("_", ($chr, $start, $transcript_id));
					$exon_num{$id}=$exon;
				}
			}
			else{
				@annotations=split('\(', $eff);
				$annotation=$annotations[1];
				$annotation=~s/\)//g;
				@annotations=split('\|', $annotation);
				$transcript_id=$annotations[8];
				$exon=$annotations[9];
				$id=join("_", ($chr, $start, $transcript_id));
				$exon_num{$id}=$exon;
			}
		}
		else{
			if($eff=~m/\,/g){
				@effs=split(',', $eff);
				for(my $i=0; $i < scalar(@effs); $i++){
					@annotations=split('\(', $effs[$i]);
					$annotation=$annotations[1];
					$annotation=~s/\)//g;
					@annotations=split('\|', $annotation);
					$transcript_id=$annotations[8];
					$exon=$annotations[9];
					$id=join("_", ($chr, $start, $transcript_id));
					$exon_num{$id}=$exon;
				}
			}
			else{
				@annotations=split('\(', $eff);
				$annotation=$annotations[1];
				$annotation=~s/\)//g;
				@annotations=split('\|', $annotation);
				$transcript_id=$annotations[8];
				$exon=$annotations[9];
				$id=join("_", ($chr, $start, $transcript_id));
				$exon_num{$id}=$exon;
			}
		}
	}
}

while(<IN1>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		$_=~s/\#//g;
		@headers=split('\t', $_);
		$col_num=scalar(@headers)-2;
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"CHROM"};
		$start=$header_rec{"ID"};
		$ref=$header_rec{"REF"};
		$alt=$header_rec{"ALT"};
		$transcript_id=$header_rec{"Transcript_ID"};
		$id2=join("_", ($chr,$start,$transcript_id));
		$id3=join("_", ($chr,$start,$ref,$alt,$transcript_id));
		$no_exon_rec=join("\t", @recs[0..$col_num]);
		if(exists $exon_num{$id2}){
			$new_recs{$id3}=join("\t", ($no_exon_rec, $exon_num{$id2}));
		}
		else{
			$new_recs{$id3}=$_;
		}
	}
}

while(<IN3>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		printf OUT "$_\n";
		$_=~s/#//g;
		@headers=split('\t', $_);
	}
	else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"CHROM"};
		$start=$header_rec{"ID"};
		$ref=$header_rec{"REF"};
		$alt=$header_rec{"ALT"};
		$transcript_id=$header_rec{"Transcript_ID"};
		$id4=join("_", ($chr,$start,$ref,$alt,$transcript_id));
		if(exists $new_recs{$id4}){
			printf OUT "$new_recs{$id4}\n";
		}
		else{
			printf OUT "$_\n";
		}
	}
}

