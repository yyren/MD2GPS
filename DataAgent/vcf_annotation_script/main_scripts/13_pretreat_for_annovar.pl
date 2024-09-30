#!/usr/bin/perl -w
use strict;

my $in_file_1=$ARGV[0];
my $out_file_1=$ARGV[1];

my @headers;
my @recs;
my $header;
my %header_rec=();
my $chr;
my $start;
my $end;
my $ref;
my $alt;
my $numRef;
my $numAlt;
my $newstart;
my $newref;
my $newalt;
my $anno_input1;
my $anno_input2;
my ($anno_input3,$id_rec);
my %snp=();
my %del=();
my %ins=();
my @arryAlt;
open(IN, "<", $in_file_1) or die "can not open < $in_file_1: $!";
open(OUT, ">", $out_file_1) or die "cannot open < $out_file_1: $!";
 while (<IN>) {
	chomp $_;
    $_=~s/\r|\n|\#//g;   
    if ($_=~m/CHROM/) {
		@headers=split("\t", $_);
    }
    else{
		@recs=split("\t", $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$header=$headers[$i];
			$header_rec{$header}=$recs[$i];
		}
		$chr=$header_rec{"CHROM"};
		$chr=~s/chr//g;
		$start=$header_rec{"ID"};
		$ref=$header_rec{"REF"};
		$alt=$header_rec{"ALT"};
		$numRef=length($ref);
		$numAlt=length($alt);
		if(($numRef==$numAlt)&&($numRef==1)){#if it is snp
			my $anno_input1 = join("\t",$chr,$start,$start,$ref,$alt);
			$snp{$anno_input1}=$_;
			printf OUT "$anno_input1\n";
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
				my $newline = join("\t",($chr,$start,$start,$ref,$alt));
				print OUT $newline."\n";
			}elsif($numRef==$numAlt){#if it is mnp
				$end=$start+$numRef-1;
				my $newline=join("\t",($chr,$start,$end,$ref,$alt));
				printf OUT "$newline\n";
			}elsif((($numRef>$numAlt)&&($numAlt>1))||(($numRef<$numAlt)&&($numRef>1))){#it is complex indel, annovar can not annotated
				# $newref=$ref;
				# $newref=~ s/^[A-Z]{1}//;
				# $newalt=$alt;
				# $newalt=~ s/^[A-Z]{1}//;
				# $newstart=$start+1;
				$end=$start+$numRef-1;
				$id_rec=join("\t",($chr,$start,$end,$ref,$alt));
				printf OUT "$id_rec\n";
			}elsif($numRef>$numAlt){#if it is deletion
				$newref=$ref;
				$newref=~ s/^[A-Z]{1}//;
				$newstart=$start+1;
				$end=$newstart+$numRef-$numAlt-1;
				$anno_input2=join("\t",$chr,$newstart,$end,$newref,'-');
				print OUT "$anno_input2\n";
			}elsif($numRef<$numAlt){#如果是插入
				$newalt=$alt;
				$newalt=~ s/^[A-Z]{1}//;
				$anno_input3 = join("\t",$chr,$start,$start,'-',$newalt);
				print OUT "$anno_input3\n";
			}
	
		}
	}
}
close IN;   
close OUT;
