#!/usr/bin/perl
#chrom start id ref alt strand sequence
my $in_file=$ARGV[0];
my $raw_occur_left_right=$ARGV[1];
my $out_file=$ARGV[2];

my(@recs,$start,$ref,$alt,$sequence,$ref_end,$alt_end,$ref_length,$alt_length);
my($ref_str,$right_pos,$right_str,$id,$strand);
my ($ref_end,$alt_end,$ref_ancor,$ref_end,$ref_ancor,$ref_ancor_pos,$ref_end_pos,$left_pos,$new_ref,$new_alt,$new_alt_end_base);
my $new_ref_end_base;

open(IN,"<$in_file") or die "can not open $in_file\n";
open(OUT,">$out_file") or die "can not open $out_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	@recs=split('\t',$_);
	if($_=~m/CHROM/g){
		#print "$_\n";
		$idx=$#recs-2;
		#print "idx:$idx\n";
		$id=join("\t",@recs[0..$idx]);
		printf OUT "$id\n";
		#print "$id\n";
	}else{
		$chr=$recs[0];
		$start=$recs[1];
		$ref=$recs[3];
		$alt=$recs[4];
		$strand=$recs[10];
		$sequence=$recs[11];
		$ref_end=substr($ref,-1);
		$alt_end=substr($alt,-1);
		$ref_length=length($ref);
		$alt_length=length($alt);
		if($raw_occur_left_right eq 'left'){
			if($strand eq '+'){
				if($ref_end eq $alt_end){
					#the variants occure at the most right positaion, do nothing for the variants located in the positive strand
					$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
					printf OUT "$id\n";
				}else{ #the variants occure at the most left positaion, need change to occure at the most right positaion.
					if(($ref_length > $alt_length)&&($alt_length==1)&&($sequence ne 'pass')){ #deletion
						for(my $i=1;$i<100;$i++){
							$ref_str=substr($ref,1,1);
							#print "ref:ref_str|i:$ref|$ref_str|$i\n";
							#print "1:$ref_str\n";
							$right_pos=$ref_length-1+$i;
							$right_str=substr($sequence,$right_pos,1);
							#print "ref|right|right_pos:$ref_str|$right_str|$right_pos\n";
							if($ref_str eq $right_str){#the new sequences are the same after move to the right
								$start=$start+1;
								$ref=substr($sequence,$i,$ref_length);
								$alt=substr($sequence,$i,1);
								#print "2:$ref\n";
							}else{
								# do not move to right 
								last;
							}
						}
						$gt="$ref".'/'."$alt";
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,@recs[9..$idx]));
						printf OUT "$id\n";
					}elsif(($ref_length < $alt_length)&&($ref_length==1)&&($sequence ne 'pass')){ #insertion
						$alt_str=substr($alt,1,($alt_length-1));
						for(my $i=1;$i<100;$i=$i+$alt_length-1){#if the insertion base/bases is/are the minimum repeat unit(refseq:GATATATC ref:G Alt:GAT)
							$right_pos=$i;
							$right_str=substr($sequence,$right_pos,($alt_length-1));
							#print "alt:alt_str|i:$alt|$alt_str|$i\n";
							$ins_ref=join("",($alt_str,$right_str));
							$ref_ins=join("",($right_str,$alt_str));
							#print "ins_ref|ref_ins:$ins_ref|$ref_ins\n";
							if($ins_ref eq $ref_ins){
								$start=$start+$alt_length-1;
								$ref=substr($sequence,$i,1);
								$alt=join("",($ref,$alt_str));
							}else{
								# do not move to right 
								last;
							}
						}
						$gt="$ref".'/'."$alt";
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,@recs[9..$idx]));
						printf OUT "$id\n";
					}else{
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
						printf OUT "$id\n";
					}
				}
			}else{
				$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
				printf OUT "$id\n";
				#do nothing because of standard format or do not know which strand located in
			}
		}elsif($raw_occur_left_right eq 'right'){
			if($strand eq "+"){
				#the variants located in negative strand occure at the most right positaion, do nothing for the variants located in the negative strand
				$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
				printf OUT "$id\n";
			}else{
				if($ref_end ne $alt_end){
					#the variants occure at the most left positaion, do nothing for the variants located in negative strand
					$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
					printf OUT "$id\n";
				}else{ #the variants occure at the most right positaion, need change to occure at the most left positaion.
					if(($ref_length > $alt_length)&&($alt_length==1)&&($sequence ne 'pass')){ #deletion
						for(my $i=1;$i<100;$i++){
							$left_pos=-$i-1;
							$left_base=substr($sequence,$left_pos,1);
							$new_ref=substr($sequence,$left_pos,$ref_length);
							$new_ref_end_base=substr($new_ref,-1,1);
							$new_alt=substr($new_ref,0,1);
							if($new_ref_end_base eq $new_alt){
								$start=$start-1;
								$ref=$new_ref;
								$alt=$new_alt;
							}else{
								$start=$start;
								$ref=$new_ref;
								$alt=$new_alt;
								# do not move to right 
								last;
							}
						}
						$gt="$ref".'/'."$alt";
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,@recs[9..$idx]));
						printf OUT "$id\n";
					}elsif(($ref_length < $alt_length)&&($ref_length==1)&&($sequence ne 'pass')){ #insertion
						$ins_base=substr($alt,1,($alt_length-1));
						for(my $i=$ref_length;$i<100;$i=$i+$alt_length-1){
							$left_pos=-$i-$ref_length;
							$new_alt=substr($sequence,$left_pos,$alt_length);
							$new_alt_end_base=substr($new_alt,-1,1);
							$new_ref=substr($new_alt,0,1);
							if($new_ref eq $new_alt_end_base){
								$start=$start-$alt_length+1;
								$ref=$new_ref;
								$alt=$new_alt;
							}else{
								$start=$start-$alt_length+1;
								$ref=$ref=$new_ref;
								$alt=$new_alt;
								# do not move to right 
								last;
							}
						}
						$gt="$ref".'/'."$alt";
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..7],$gt,@recs[9..$idx]));
						printf OUT "$id\n";
					}else{
						$id=join("\t",($chr,$start,$recs[2],$ref,$alt,@recs[5..$idx]));
						printf OUT "$id\n";
					}
				}
			}
		}
	}
}
close IN;
close OUT;
