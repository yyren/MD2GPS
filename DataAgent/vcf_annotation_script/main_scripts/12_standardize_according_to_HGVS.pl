#!/usr/bin/perl -w
#query_file format:
#CHROM	ID	Database_infor	REF	ALT	VarFreq_precent	Cov	Var_Cov	Genome_Type	Source

#bed format:
#CHROM	start	end	strand	enst_id	offical_genename

#reference fasta file:ucsc.fasta
use strict;
use Bio::DB::Fasta;

my $in_query_file=$ARGV[0];
my $in_bed_file=$ARGV[1];
my $in_reference=$ARGV[2];
my $out_file=$ARGV[3];

#my (%strand_of_pos,$strand,$id_rec_adjust,$chr);
my (%strand_of_pos,$id_rec_adjust,$strand,$bed_rec,$query_rec,$bed_chr,%bed_hash);
&storeBedInfor();
my $db = Bio::DB::Fasta->new($in_reference);
#standardize the coordinate, ref and alt based the standard of HGVS: the variant accour at the most 3' of the transcript
open(IN,"<$in_query_file") or die "can not open $in_query_file:$!";
open(OUT,">$out_file") or die "can not open $out_file:$!";
while($query_rec=<IN>){
	$strand='NA';
	chomp $query_rec;
	$query_rec=~s/\r|\n//g;
	my @recs=split('\t',$query_rec);
	if($query_rec=~m/CHROM/){
		printf OUT "$query_rec\tStrand\tAdjusted\n";
	}else{
		my $ref_length=length($recs[3]);
		my $alt_length=length($recs[4]);
		my $chr=$recs[0];
		my $start=$recs[1];
		my $ref=$recs[3];
		$ref=~tr/atcg/ATCG/;
		my $alt=$recs[4];
		$alt=~tr/atcg/ATCG/;
		$bed_chr=$recs[0];
		$bed_chr=~s/chr//;
		my $ref_final_base=substr($recs[3],-1);
		my $alt_final_base=substr($recs[4],-1);
		my $id_pos=join("_",@recs[0..1]);
		my $adjust_or_not='no';
		foreach my $bed_start(sort {$a<=>$b} keys %{$bed_hash{$bed_chr}{'+'}}){
			next if ($start < $bed_start);
			if(($start>=$bed_start)&&($start<=$bed_hash{$bed_chr}{'+'}{$bed_start})){
				$strand='+';
				last;
			}
			if($bed_start > $start){
				last;
			}
		}
		if($strand eq 'NA'){
			foreach my $bed_start(sort {$a<=>$b} keys %{$bed_hash{$bed_chr}{'-'}}){
				next if ($start < $bed_start);
				if(($start>$bed_start)&&($start<=$bed_hash{$bed_chr}{'-'}{$bed_start})){
					$strand='-';
					last;
				}
				if($bed_start > $start){
					last;
				}
			}
		}
		#print "strand:$strand\n";
		if($strand ne 'NA'){#adjust the coordinate only if the raw seq and adjusted seq has the same sequence.
			if(length($recs[3])==length($recs[4])){#snp or mnp
				$id_rec_adjust=join("\t",(@recs,'NA','no'));
				#printf OUT "$id_rec_adjust\n";
			}elsif((length($recs[3]) > length($recs[4])) && ((length($recs[4])==1))){#deletion	
				my $sequence= sequence_from_bed($bed_chr,$start,$ref,$strand,'del',$in_reference);#extract the sequence with 200 bases
				my $ref_length=length($ref);
				if($strand eq '-'){
					my $raw_ref_seq=substr($sequence,-100,(100-$ref_length)).substr($sequence,-1,1);
					for(my $i=1;$i<100;$i++){
						my $new_ref_seq=substr($sequence,-100,(100-$ref_length-$i)).substr($sequence,-$i);
						my $new_ref=substr($sequence,-($ref_length+$i),$ref_length);
						my $new_alt=substr($new_ref,0,1);
						if($new_ref_seq eq $raw_ref_seq){
							$start=$start-1;
							$ref=$new_ref;
							$alt=$new_alt;
							$adjust_or_not='yes';
						}else{
							#do not move to the left
							last;
						}
					}
					my $gt="$ref".'/'."$alt";
					$id_rec_adjust=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,$recs[9],$strand,$adjust_or_not));
				}elsif($strand eq '+'){
					my $raw_ref_seq=substr($sequence,0,1).substr($sequence,$ref_length,(100-$ref_length));
					for(my $i=1;$i<100;$i++){
						my $new_ref_seq=substr($sequence,0,($i+1)).substr($sequence,($i+$ref_length),(100-$ref_length-$i));
						my $new_ref=substr($sequence,$i,$ref_length);
						my $new_alt=substr($new_ref,0,1);
						if($new_ref_seq eq $raw_ref_seq){
							$start=$start+1;
							$ref=$new_ref;
							$alt=$new_alt;
							$adjust_or_not='yes';
						}else{
							#do not move to the right
							last;
						}
					}
					my $gt="$ref".'/'."$alt";
					$id_rec_adjust=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,$recs[9],$strand,$adjust_or_not));
				}else{
					$id_rec_adjust=join("\t",(@recs,'NA','no'));
				}
			}elsif((length($recs[3]) < length($recs[4])) && (length($recs[3])==1)){#insertion
				my $sequence= sequence_from_bed($bed_chr,$start,$ref,$strand,'ins',$in_reference);#extract the sequence with 200 bases
				#print "$bed_chr,$start,$ref,$strand,$in_reference\n$sequence\n";
				my $ins_alt_seq=substr($alt,1,($alt_length-1));
				if($strand eq '-'){
					my $raw_alt_seq=substr($sequence,-100,100)."$ins_alt_seq";
					for(my $i=1;$i<100;$i++){
						my $new_alt_seq=substr($sequence,-100,(100-$i))."$ins_alt_seq".substr($sequence,(-$i));
						my $new_ref=substr($sequence,(-$i-1),1);
						my $new_alt="$new_ref"."$ins_alt_seq";
						if($new_alt_seq eq $raw_alt_seq){
							$start=$start-1;
							$ref=$new_ref;
							$alt=$new_alt;
							$adjust_or_not='yes';
						}else{
							#do not move to the left
							last;
						}
					}
					if($adjust_or_not ne 'yes'){
						my $ins_length=length($ins_alt_seq);
						my $alt_seq_with_dup=substr($sequence,-100,(100-$ins_length))."$ins_alt_seq".substr($sequence,-($ins_length));#CCA(GCGTGGACA->ins in this pos)[GCGTGGACA]ACCC
						my $ref_dup_seq=substr($sequence,-($ins_length));
						if(($alt_seq_with_dup eq $raw_alt_seq)&&($ref_dup_seq eq $ins_alt_seq)){
							$start=$start-$ins_length;
							$ref=substr($sequence,-($ins_length+1),1);
							$alt=$ref.$ins_alt_seq;
							#print "$start|$ref|$alt\n";
							for(my $i=1;$i<=$ins_length;$i++){
								my $ref_dup_most_left_base=substr($ref_dup_seq,-($i),1);
								my $ref_seq_detect_base=substr($sequence,-($ins_length+$i),1);
								if($ref_dup_most_left_base eq $ref_seq_detect_base){
									$ref=substr($sequence,-($ins_length+$i+1),1);
									$alt=$ref.substr($sequence,-($ins_length+$i),$ins_length);
									$start=$start-1;
									$adjust_or_not='yes';
								}else{
									last;
								}
							
							}
							
						}
					}
					my $gt="$ref".'/'."$alt";
					$id_rec_adjust=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,$recs[9],$strand,$adjust_or_not));
				}elsif($strand eq '+'){
					my $raw_alt_seq=substr($sequence,0,1)."$ins_alt_seq".substr($sequence,1,99);
					for(my $i=1;$i<100;$i++){
						my $new_alt_seq=substr($sequence,0,($i+1))."$ins_alt_seq".substr($sequence,($i+1),(100-$i-1));
						my $new_ref=substr($sequence,$i,1);
						my $new_alt="$new_ref"."$ins_alt_seq";
						if($new_alt_seq eq $raw_alt_seq){
							$start=$start+1;
							$ref=$new_ref;
							$alt=$new_alt;
							$adjust_or_not='yes';
						}else{
							#do not move to the left
							last;
						}
					}
					if($adjust_or_not ne 'yes'){#'chr7	55249002	.	C	CCAGCGTGGA' need to be 'chr7	55249011	.	A	ACAGCGTGGA'
						my $alt_seq_with_dup=substr($sequence,0,$alt_length)."$ins_alt_seq".substr($sequence,$alt_length,(100-$alt_length));
						my $ins_length=length($ins_alt_seq);
						my $ref_dup_seq=substr($sequence,1,$ins_length);
						if(($alt_seq_with_dup eq $raw_alt_seq)&&($ref_dup_seq eq $ins_alt_seq)){
							$start=$start+$ins_length;
							$ref=substr($ref_dup_seq,-1);
							$alt=$ref.$ins_alt_seq;
							for(my $i=0;$i<$ins_length;$i++){
								if(substr($ref_dup_seq,$i,1) eq substr($sequence,($alt_length+$i),1)){
									$ref=substr($ref_dup_seq,$i,1);
									$alt=$ref.substr($ref_dup_seq,($i+1),($ins_length-1-$i)).substr($ref_dup_seq,0,($i+1));
									$start=$start+1;
									$adjust_or_not='yes';
								}else{
									last;
								}
							
							}
						}
					}
					my $gt="$ref".'/'."$alt";
					$id_rec_adjust=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,$recs[9],$strand,$adjust_or_not));
				}else{
					$id_rec_adjust=join("\t",(@recs,'NA','no'));
				}
			}else{#complex indel
				$id_rec_adjust=join("\t",(@recs,'NA','no'));
			}
		}else{
			$id_rec_adjust=join("\t",(@recs,'NA','no'));
		}
		printf OUT "$id_rec_adjust\n";
	}
}
close IN;
close OUT;

sub sequence_from_bed{
	my ($func_chr,$func_start,$func_ref,$func_strand,$func_ins_or_del,$func_reference)=@_;
	#print "$func_ins_or_del\n";
	my($bed_start,$bed_end);
	if($func_strand eq '+'){
		$bed_start=$func_start;
		$bed_end=$func_start+199;
	}elsif($func_strand eq '-'){
		if($func_ins_or_del eq 'ins'){
			$bed_end=$func_start;
			$bed_start=$bed_end-200;
		}elsif($func_ins_or_del eq 'del'){
			$bed_end=$func_start+length($func_ref)-1;
			$bed_start=$bed_end-200;
		}
	}else{
		print "error: no strand\n";
	}
	my $seq=$db->seq($func_chr, $bed_start => $bed_end);
	#print "$func_chr|$bed_start|$bed_end|$seq\n";
	chomp ($seq);
	return($seq);
}

sub storeBedInfor{
	open(IN,"<$in_bed_file") or die "can not open $in_bed_file:$!";
	while(<IN>){
		chomp $_;
		my @recs=split('\t',$_);
		$recs[0]=~s/chr//;
		$bed_hash{$recs[0]}{$recs[3]}{$recs[1]}=$recs[2];
	}
	close IN;
}