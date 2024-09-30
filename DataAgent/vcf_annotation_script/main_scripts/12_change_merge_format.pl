use strict;
use warnings;
use Bio::DB::Fasta;

my $in_file=$ARGV[0];
my $out_file=$ARGV[1];
my $reference=$ARGV[2];
my $db = Bio::DB::Fasta->new($reference);
#my ($source_tag,@recs,$start_pos,$end_pos,$ref_length,$new_start,$new_ref,$new_alt,$ref_end_pos,$alt_end_pos,$ref_start_pos,$alt_start_pos);
#my (@genotypes,$gt,$gt_ref,$gt_alt);
open(IN,"<$in_file")or die "can not open $in_file\n";
open(OUT,">$out_file") or die "can not open $out_file\n";
while(<IN>){
	chomp $_;
	if($_=~m/CHROM/){
		printf OUT "$_\n";
	}else{
		my @recs=split('\t', $_);
		my $ref_length=length($recs[3]);
		my $alt_length=length($recs[4]);
		my $rec_id='';
		my $genotype='';
		my $gt_type='';
		if($recs[8]=~m/\//){
			my @gts=split('\/',$recs[8]);
			if($gts[0] eq $gts[1]){
				if($gts[0] eq $recs[3]){
					$gt_type='0/0';
				}else{
					$gt_type='1/1';
				}
			}else{
				$gt_type='0/1';
			}
		}else{
			$gt_type='0/1';
		}
		#print "$ref_length\n";
		if($ref_length==$alt_length){#snp or mnp
			if($ref_length==1){#snp
				$rec_id=join("\t",@recs);
			}else{#mnp
				for(my $i=0;$i<scalar(@recs);$i++){#remove the same base from left
					my $ref_one_str=substr($recs[3],0,1);
					my $alt_one_str=substr($recs[4],0,1);
					if($ref_one_str eq $alt_one_str){
						$recs[1]=$recs[1]+1;
						substr($recs[3],0,1)="";
						substr($recs[4],0,1)="";
					}else{
						last;
					}
				}
				$ref_length=length($recs[3]);
				for(my $i=$ref_length;$i>0;$i--){#remove the same base from right
					my $ref_one_str=substr($recs[3],-1);
					my $alt_one_str=substr($recs[4],-1);
					if($ref_one_str eq $alt_one_str){
						substr($recs[3],-1,1)="";
						substr($recs[4],-1,1)="";
					}else{
						last;
					}
				}
				$recs[8]=Genotype_result($recs[3],$recs[4],$gt_type);
				$rec_id=join("\t",@recs);
			}
		}elsif($ref_length>$alt_length){#deletion
			if($alt_length==1){#ATT->A
				#do nothing
			}else{#ACCT->ACT:AC->A #ACCT->ACG:CCT->CG #ACCG -> CG
				for(my $i=$alt_length;$i>0;$i--){# remove same base from right end position firstly
					my $ref_one_str=substr($recs[3],-1);
					my $alt_one_str=substr($recs[4],-1);
					if($ref_one_str eq $alt_one_str){
						substr($recs[3],-1,1)="";
						substr($recs[4],-1,1)="";
						#print "1:$recs[3]|$recs[4]\n";
					}else{
						last;
					}
				}
				if($recs[4] eq ''){
					$alt_length=0;
				}else{
					$alt_length=length($recs[4]);
				}
				#print "$alt_length\n";
				if($alt_length==1){#ATT->A
					#print "2:$recs[3]|$recs[4]|$alt_length\n";
					#do nothing
				}elsif($alt_length==0){#ACCG -> CG: TAC->T
					#print "alt\n";
					my $chr_temp=$recs[0];
					$chr_temp=~s/chr//;
					my $start_temp=$recs[1]-1;
					my $end_temp=$recs[1]-1;
					my $ancor_base=$db->seq($chr_temp, $start_temp => $end_temp);
					chomp($ancor_base);
					$recs[3]=$ancor_base.$recs[3];
					$recs[4]=$ancor_base;
					$recs[1]=$recs[1]-1;
				}else{#ACCG ->ACT; TAC->AA
					$ref_length=length($recs[3]);
					$alt_length=length($recs[4]);
					my ($ref_temp,$alt_temp);
					my $trim_str='no';
					my $pos_raw=$recs[1];
					my $ref_raw=$recs[3];
					my $alt_raw=$recs[4]; 
					for(my $i=0;$i<$ref_length;$i++){# remove same base from left position
						my $ref_one_str=substr($recs[3],0,1);
						my $alt_one_str=substr($recs[4],0,1);
						if($ref_one_str eq $alt_one_str){
							$ref_temp=$recs[3];
							$alt_temp=$recs[4];
							substr($recs[3],0,1)="";
							substr($recs[4],0,1)="";
							$trim_str='yes'
						}else{
							if($trim_str eq 'yes'){
								$recs[3]=$ref_temp;
								$recs[4]=$alt_temp;
								my $adjust_numb=$ref_length-(length($recs[3]));
								$recs[1]=$recs[1]+$adjust_numb;
								last;
							}
						}
					}
					if($trim_str eq 'no'){#TAC->AA: GTAC->GAA
						$recs[1]=$pos_raw-1;
						my $chr_temp=$recs[0];
						$chr_temp=~s/chr//;
						my $ancor_base=$db->seq($chr_temp, $recs[1] => $recs[1]);
						chomp($ancor_base);
						$recs[3]=$ancor_base.$ref_raw;
						$recs[4]=$ancor_base.$alt_raw;
					}
				}
			}
			$recs[8]=Genotype_result($recs[3],$recs[4],$gt_type);
			$rec_id=join("\t",@recs);
		}else{#insertion
			if($ref_length==1){#A->ATT
				#do nothing
			}else{#ACT->ACCT:A->AC #ACG->ACCT:CG->CCT #CG->ACCG 
				for(my $i=$ref_length;$i>0;$i--){# remove same base from right end position firstly
					my $ref_one_str=substr($recs[3],-1);
					my $alt_one_str=substr($recs[4],-1);
					if($ref_one_str eq $alt_one_str){
						substr($recs[3],-1,1)="";
						substr($recs[4],-1,1)="";
					}else{
						last;
					}
				}
				if($recs[3] eq ''){
					$ref_length=0;
				}else{
					$ref_length=length($recs[3]);
				}
				if($ref_length==1){#A->ATT
					#do nothing
				}elsif($ref_length==0){#CG->ACCG: T->TAC
					my $chr_temp=$recs[0];
					$chr_temp=~s/chr//;
					my $start_temp=$recs[1]-1;
					my $end_temp=$recs[1]-1;
					my $ancor_base=$db->seq($chr_temp, $start_temp => $end_temp);
					chomp($ancor_base);
					$recs[3]=$ancor_base;
					$recs[4]=$ancor_base.$recs[4];
					$recs[1]=$recs[1]-1;
				}else{#ACT->ACCG, ACC ->ACCG
					$ref_length=length($recs[3]);
					$alt_length=length($recs[4]);
					my ($ref_temp,$alt_temp);
					for(my $i=0;$i<$alt_length;$i++){# remove same base from left position
						my $ref_one_str=substr($recs[3],0,1);
						my $alt_one_str=substr($recs[4],0,1);
						if($ref_one_str eq $alt_one_str){
							$ref_temp=$recs[3];
							$alt_temp=$recs[4];
							substr($recs[3],0,1)="";
							substr($recs[4],0,1)="";
						}else{
							$recs[3]=$ref_temp;
							$recs[4]=$alt_temp;
							my $adjust_numb=$ref_length-(length($recs[3]));
							$recs[1]=$recs[1]+$adjust_numb;
							last;
						}
					}
				}
			}
			$recs[8]=Genotype_result($recs[3],$recs[4],$gt_type);
			$rec_id=join("\t",@recs);
		}
		printf OUT "$rec_id\n";
	}
}

sub Genotype_result{
	my ($ref_temp,$alt_temp,$gt_types)=@_;
	my $gt_result='';
	if($gt_types eq '1/1'){
		$gt_result=$alt_temp.'/'.$alt_temp;
	}elsif($gt_types eq '0/0'){
		$gt_result=$ref_temp.'/'.$ref_temp;
	}elsif($gt_types eq '0/1'){
		$gt_result=$ref_temp.'/'.$alt_temp;
	}
	return($gt_result);
}

close IN;
close OUT;