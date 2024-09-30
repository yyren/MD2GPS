#!/usr/bin/perl -w
use strict;

my $in_snpeff=$ARGV[0];
my $in_annovar=$ARGV[1];
my $out_merge=$ARGV[2];
my @headers;
my $annovar_header;
my @recs;
my $header;
my $anno_rec;
my %annova_header_rec=();
my $chr;
my $start;
my $end;
my $ref;
my $alt;
my $anno_id;
my %anno_recs=();
my @snpeff_headers;
my %snpeff_header_rec=();
my $numRef;
my $numAlt;
my $newstart;
my $newref;
my $newalt;
my @arryAlt;
my $snpeff_id;
my $header_1000g;
my $header_cosmic;
my $header_clinvar;
my $header_snp;
my $one_k_score;
my $sift_scrore;
my $poly_hdiv_score;
my $poly_hdiv_pred;
my $poly_hvar_score;
my $poly_hvar_pred;
my $cadd_raw;
my $cadd_phred;
my $cosmic;
my $clinvar;
my $snp;
my $header_chr;
my $header_1000g_all;
my $one_k_all_score;
my $one_k_header_same_numb=0;
###use hash to store the column to be merged to snpEff
open(IN, "<$in_annovar") or die "can not open < $in_annovar\n";
while (<IN>) {
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/Chr/) {
		@headers=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			if($headers[$i]=~m/1000g2014oct\_all/){
				$one_k_header_same_numb=$one_k_header_same_numb+1;
			}
		}
		#print "$one_k_header_same_numb\n";
		@headers=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			if($one_k_header_same_numb==1){
				if($header eq '1000g2014oct_all'){
					$header_1000g_all=$header;
					print "all:$header_1000g_all\n";
				}elsif($header=~m/1000g/){
					$header_1000g=$header;
					print "$header_1000g\n";
				}
			}else{
				if($header=~m/1000g/){
					$header_1000g_all=$header;
					$header_1000g=$header;
				}
			}
			if($header=~m/cosmic/){
				$header_cosmic=$header;
			}
			elsif($header=~m/clinvar/){
				$header_clinvar=$header;
			}
			elsif($header=~m/snp/){
				$header_snp=$header;
			}
		}
		$annovar_header=join("\t",($header_1000g,'1k_genome_all','SIFT_score','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','CADD_raw','CADD_phred',$header_cosmic,$header_clinvar,$header_snp));
	}else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@recs);$i++){
			$header=$headers[$i];
			$annova_header_rec{"$header"}=$recs[$i];
		}
		$chr=$annova_header_rec{"Chr"};
		if($chr=~m/chr/g){
			$chr=~s/chr//g;
		}
		$start=$annova_header_rec{"Start"};
		$end=$annova_header_rec{"End"};
		$ref=$annova_header_rec{"Ref"};
		$alt=$annova_header_rec{"Alt"};
		$one_k_score=$annova_header_rec{"$header_1000g"};
		$one_k_all_score=$annova_header_rec{"$header_1000g_all"};
		$sift_scrore=$annova_header_rec{"SIFT_score"};
		$poly_hdiv_score=$annova_header_rec{"Polyphen2_HDIV_score"};
		$poly_hdiv_pred=$annova_header_rec{"Polyphen2_HDIV_pred"};
		$poly_hvar_score=$annova_header_rec{"Polyphen2_HVAR_score"};
		$poly_hvar_pred=$annova_header_rec{"Polyphen2_HVAR_pred"};
		$cadd_raw=$annova_header_rec{"CADD_raw"};
		$cadd_phred=$annova_header_rec{"CADD_phred"};
		$cosmic=$annova_header_rec{"$header_cosmic"};
		$clinvar=$annova_header_rec{"$header_clinvar"};
		$snp=$annova_header_rec{"$header_snp"};
		$anno_id=join("_",($chr,$start,$end,$ref,$alt));
		$anno_rec=join("\t",($one_k_score,$one_k_all_score,$sift_scrore,$poly_hdiv_score,$poly_hdiv_pred,$poly_hvar_score,$poly_hvar_pred,$cadd_raw,$cadd_phred,$cosmic,$clinvar,$snp));
		$anno_recs{$anno_id}=$anno_rec;
		#if($anno_id eq'2_212652764_212652764_T_-'){
		#print "1:	$anno_id\n";
		#}
	}
}
close IN;
###change snpeff variants format(chrom start end ref alt) to annovar format and then map and merge
open(IN,"<$in_snpeff") or die "can not open $in_snpeff\n";
open(OUT,">$out_merge") or die "can not open $out_merge\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	if($_=~m/CHROM/g){
		@snpeff_headers=split('\t', $_);
		printf OUT "$_\t$annovar_header\n";
	}else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@recs);$i++){
			$header=$snpeff_headers[$i];
			if($header=~m/CHROM/g){
				$header_chr=$header;
			}
			$snpeff_header_rec{"$header"}=$recs[$i];
		}
		$chr=$snpeff_header_rec{"$header_chr"};
		$chr=~s/chr//g;
		#print "2:	$chr\n";
		$start=$snpeff_header_rec{"ID"};
		#$end=$annova_header_rec{"End"};
		$ref=$snpeff_header_rec{"REF"};
		$alt=$snpeff_header_rec{"ALT"};
		$numRef=length($ref);
		$numAlt=length($alt);
		if(($numRef==$numAlt)&&($numRef==1)){#if it is snp
			$snpeff_id = join("_",$chr,$start,$start,$ref,$alt);
		}
		else{
			if($alt=~m/\,/){
				@arryAlt=split("\,",$alt);
				if($ref eq $arryAlt[0]){
					$alt = $arryAlt[1];
				}
				else{
					$alt = $arryAlt[0];
				}
				my $snpeff_id = join("\t",($chr,$start,$start,$ref,$alt));
			}elsif($numRef==$numAlt){#if it is mnp
				$end=$start+$numRef-1;
				$snpeff_id=join("_",($chr,$start,$end,$ref,$alt));
			}elsif((($numRef>$numAlt)&&($numAlt>1))||(($numRef<$numAlt)&&($numRef>1))){#it is complex indel, annovar can not annotated
				# $newref=$ref;
				# $newref=~ s/^[A-Z]{1}//;
				# $newalt=$alt;
				# $newalt=~ s/^[A-Z]{1}//;
				# $newstart=$start+1;
				$end=$start+$numRef-1;
				$snpeff_id=join("_",($chr,$start,$end,$ref,$alt));
				#print "2:$snpeff_id\n";
			}elsif($numRef>$numAlt){#if it is deletion
				$newref=$ref;
				$newref=~ s/^[A-Z]{1}//;
				$newstart=$start+1;
				$end=$newstart+$numRef-$numAlt-1;
				$snpeff_id=join("_",$chr,$newstart,$end,$newref,'-');
			}elsif($numRef<$numAlt){#如果是插入
				$newalt=$alt;
				$newalt=~ s/^[A-Z]{1}//;
				$snpeff_id = join("_",$chr,$start,$start,'-',$newalt);
			}
		}
		#print "3:	$snpeff_id\n";
		if(exists $anno_recs{$snpeff_id}){#mapping and merge the same variants from the snpEff and annovar by using hash
			printf OUT "$_\t$anno_recs{$snpeff_id}\n";
		}
	}
}

close IN;
close OUT;
