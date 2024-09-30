#filter variants:
use strict;
use warnings;


my $in_hgmd_db=$ARGV[0];
my $in_clinvar_db=$ARGV[1];
my $infile=$ARGV[2];
my $outfile=$ARGV[3];
my $out_not_pass=$ARGV[4];

my (%hgmd_variant_category,%hgmd_variant_id,%clinvar_variant_category,%clinvar_variant_id, %clinvar_variant_category_info, %clinvar_latest_evidence)=();
open(IN,"<$in_hgmd_db") or die "cannot open $in_hgmd_db";
while(my $line=<IN>){
	chomp $line;
	if($line !~m/^#/){
		my @recs=split(/\t/,$line);
		if($recs[0]=~m/chr/){
			$recs[0]=~s/chr//;
		}
		my $id_key=join("\_",@recs[0,1,3,4]);
		my @infors=split(/\;/,$recs[7]);
		if($infors[0]=~m/CLASS\=/){
			my $variant_category=(split(/\=/,$infors[0]))[1];
			$hgmd_variant_category{$id_key}=$variant_category;
			$hgmd_variant_id{$id_key}=$recs[2];
		}
	}
}
close IN;

open(IN,"<$in_clinvar_db") or die "cannot open $in_clinvar_db";
while(my $line=<IN>){
	chomp $line;
	if($line !~m/^#/){
		my @recs=split(/\t/,$line);
		if($recs[0]=~m/chr/){
			$recs[0]=~s/chr//;
		}
		my $id_key=join("\_",@recs[0,1,3,4]);
		my @infors=split(/\;/,$recs[7]);
		for(my $i=0;$i<=$#infors;$i++){
			if($infors[$i]=~m/CLNSIG\=/){
				my $variant_category=(split(/\=/,$infors[$i]))[1];
				$clinvar_variant_category{$id_key}=$variant_category;
				$clinvar_variant_id{$id_key}=$recs[2];
			}
			if($infors[$i]=~m/CLNSIGCONF\=/){
				my $variant_category=(split(/\=/,$infors[$i]))[1];
				$clinvar_variant_category_info{$id_key}=$variant_category;
			}
			if($infors[$i]=~m/CLASS\_Latest\_Evidence\=/){
				my $variant_category=(split(/\=/,$infors[$i]))[1];
				$clinvar_latest_evidence{$id_key}=$variant_category;
			}
		}
	}
}
close IN;

open(IN,"<$infile") or die "cannot open $infile";
open(OUT,">$outfile") or die "cannot open $outfile";
open(OUT2,">$out_not_pass") or die "cannot open $out_not_pass";
while(my $line=<IN>){
	chomp $line;
	if($line=~m/^Tag/){
		print OUT "$line\tHGMD_ID\tHGMD_category\tClinvar_ID\tClinvar_category\tClinvar_category_info\tClinvar_latest_evidence\n";
		print OUT2 "$line\tHGMD_ID\tHGMD_category\tClinvar_ID\tClinvar_category\tClinvar_category_info\tClinvar_latest_evidence\n";
		
	}else{
		my @recs=split(/\t/,$line);
		my $id_key=join("\_",@recs[6..9]);
		my $variant_category_hgmd='NA';
		my $hgmd_id='NA';
		my $variant_category_clinvar='NA';
		my $clinvar_id='NA';
		my $variant_category_clinvar_info='NA';
		my $clinvar_latest_one_evidence='NA';
		
		my $impact_func='filter';
		my $pop_tag='pass';
		my $pathogenic_sift_polyphen='no';
		my $func_tag='filter';
		my $pathogenic_clinvar='no';
		my $pathogenic_hgmd='no';
		my $pathogenic_sift='no';
		my $pathogenic_polyphen_hdiv='no';
		my $pathogenic_polyphen_hvar='no';
		if($recs[15] !~m/LOW/){
			$impact_func='pass';
		}
		for(my $i=25;$i<=30;$i++){
			if(($recs[$i] ne 'NA') and ($recs[$i] ne '.')){
				if($recs[$i] >= 0.01){
					$pop_tag='filter';
				}
			}
		}
		
		if(($impact_func eq 'pass') and ($pop_tag eq 'pass')){
			if(exists $hgmd_variant_category{$id_key}){
				$variant_category_hgmd=$hgmd_variant_category{$id_key};
				$hgmd_id=$hgmd_variant_id{$id_key};
			}
			if(exists $clinvar_variant_category{$id_key}){
				$variant_category_clinvar=$clinvar_variant_category{$id_key};
				$clinvar_id=$clinvar_variant_id{$id_key};
				if(exists $clinvar_variant_category_info{$id_key}){
					$variant_category_clinvar_info=$clinvar_variant_category_info{$id_key};
				}
				if(exists $clinvar_latest_evidence{$id_key}){
					$clinvar_latest_one_evidence=$clinvar_latest_evidence{$id_key};
				}
			}
			if(($variant_category_clinvar=~m/pathogenic/i) and ($variant_category_clinvar!~m/conflicting/i) and ($variant_category_clinvar ne 'NA')){#clinvar pathogenic, or likely pathogenic
				$pathogenic_clinvar='yes';
			}
			if($variant_category_hgmd ne 'NA'){#HGMD database
				$pathogenic_hgmd='yes';
			}
			if(($recs[31] ne 'NA') and ($recs[31] ne '.')){
				if($recs[31]<0.05){
					#$pathogenic_sift_polyphen='yes';
					$pathogenic_sift='yes';
				}
			}else{
				$pathogenic_sift='NA';
			}
			if(($recs[33] ne 'NA') and ($recs[33] ne '.')){
				if(($recs[33] eq 'D') or ($recs[33] eq 'P')){#polyphen
					#$pathogenic_sift_polyphen='yes';
					$pathogenic_polyphen_hdiv='yes';
				}
			}else{
				$pathogenic_polyphen_hdiv='NA';
			}
			if(($recs[35] ne 'NA') and ($recs[35] ne '.')){
				if(($recs[35] eq 'D') or ($recs[35] eq 'P')){#polyphen
					#$pathogenic_sift_polyphen='yes';
					$pathogenic_polyphen_hvar='yes';
				}
			}else{
				$pathogenic_polyphen_hvar='NA';
			}
			if((($pathogenic_sift eq 'yes') or ($pathogenic_polyphen_hdiv eq 'yes') or ($pathogenic_polyphen_hvar eq 'yes')) or (($pathogenic_sift eq 'NA') and ($pathogenic_polyphen_hdiv eq 'NA') and ($pathogenic_polyphen_hvar eq 'NA'))){
				$pathogenic_sift_polyphen='yes'
			}
			
			### filter
			if(($variant_category_clinvar=~m/pathogenic/i) and ($variant_category_clinvar!~m/conflicting/i) and ($variant_category_clinvar ne 'NA')){#clinvar: pathogenic, likely pathogenic
				$func_tag='pass';
				if($clinvar_latest_one_evidence eq 'NA'){
					$clinvar_latest_one_evidence=$variant_category_clinvar;
				}
			}elsif(($variant_category_clinvar=~m/benign/i) and ($variant_category_clinvar!~m/pathogenic/i)){#benign or likely benign
				$func_tag='filter';
				if($clinvar_latest_one_evidence eq 'NA'){
					$clinvar_latest_one_evidence=$variant_category_clinvar;
				}
			}else{#cinvar: conflicting, NA, uncertain, not_provided,
				if(($variant_category_clinvar_info!~m/benign/i) and($variant_category_clinvar_info=~m/pathogenic/i)){#likely pathogenic|uncertain
					$func_tag='pass';
					if($clinvar_latest_one_evidence eq 'NA'){
						$clinvar_latest_one_evidence=$variant_category_clinvar_info;
					}
				}elsif(($variant_category_clinvar_info=~m/benign/i) and($variant_category_clinvar_info!~m/pathogenic/i)){#likely benign|uncertain
					$func_tag='filter';
					if($clinvar_latest_one_evidence eq 'NA'){
						$clinvar_latest_one_evidence=$variant_category_clinvar_info;
					}
				}elsif($clinvar_latest_one_evidence=~m/benign/i){#likely pathogenic(2021)|likely benign(2023)
					$func_tag='filter';
					if($clinvar_latest_one_evidence eq 'NA'){
						$clinvar_latest_one_evidence=$variant_category_clinvar_info;
					}
				}elsif($clinvar_latest_one_evidence=~m/pathogenic/i){#likely pathogenic(2023)|likely benign(2021)
					$func_tag='pass';
					if($clinvar_latest_one_evidence eq 'NA'){
						$clinvar_latest_one_evidence=$variant_category_clinvar_info;
					}
				}else{#clinvar: NA, uncertain, not_provided
					if($pathogenic_hgmd eq 'yes'){#HGMD: in HGMD database
						$func_tag='pass';
					}elsif($recs[$#recs]=~m/polymorphism/){#MutationTaster: polymorphism, polymorphism|disease causing, polymorphism|no prediction, et al.
						$func_tag='filter';
					}else{#mutationTaster: disease causing, Not Found, no prediction
						if($recs[$#recs] eq 'disease causing'){
							$func_tag='pass';
						}else{#mutationTaster: Not Found, no prediction
							#print "$pathogenic_sift_polyphen\n";
							if($pathogenic_sift_polyphen eq 'yes'){
								$func_tag='pass';
							}else{
								$func_tag='filter';
							}
						}
					}
					if($clinvar_latest_one_evidence eq 'NA'){
						$clinvar_latest_one_evidence=$variant_category_clinvar;
					}
				}
			}
			# if($id_key eq '20_5294679_A_G'){
				# print "1:yes\n";
			# }
			if($func_tag eq 'pass'){
				#print "1variant_category_clinvar_info:$variant_category_clinvar_info\n";
				# if($id_key eq '20_5294679_A_G'){
					# #print "2:yes\n";
					# #print OUT "$line\t$hgmd_id\t$variant_category_hgmd\t$clinvar_id\t$variant_category_clinvar\t$variant_category_clinvar_info\t$clinvar_latest_one_evidence\n";
				# }
				print OUT "$line\t$hgmd_id\t$variant_category_hgmd\t$clinvar_id\t$variant_category_clinvar\t$variant_category_clinvar_info\t$clinvar_latest_one_evidence\n";
			}else{
				#print "2variant_category_clinvar_info:$variant_category_clinvar_info\n";
				# if($id_key eq '20_5294679_A_G'){
					# print "3:yes\n";
				# }
				print OUT2 "$line\t$hgmd_id\t$variant_category_hgmd\t$clinvar_id\t$variant_category_clinvar\t$variant_category_clinvar_info\t$clinvar_latest_one_evidence\n";
			}
		}else{
			if(exists $hgmd_variant_category{$id_key}){
				$variant_category_hgmd=$hgmd_variant_category{$id_key};
				$hgmd_id=$hgmd_variant_id{$id_key};
			}
			if(exists $clinvar_variant_category{$id_key}){
				$variant_category_clinvar=$clinvar_variant_category{$id_key};
				$clinvar_id=$clinvar_variant_id{$id_key};
				if(exists $clinvar_latest_evidence{$id_key}){
					$clinvar_latest_one_evidence=$clinvar_latest_evidence{$id_key};
				}else{
					if(exists $clinvar_variant_category_info{$id_key}){
						$clinvar_latest_one_evidence=$clinvar_variant_category_info{$id_key};
					}
				}
			}
			
			print OUT2 "$line\t$hgmd_id\t$variant_category_hgmd\t$clinvar_id\t$variant_category_clinvar\t$variant_category_clinvar_info\t$clinvar_latest_one_evidence\n";
			# if($id_key eq '20_5294679_A_G'){
				# print "no\n";
			# }
			
		}
	}
}
close IN;
close OUT;
close OUT2;