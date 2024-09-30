#!/usr/bin/perl
use strict;
use warnings;

my $in_file=$ARGV[0];
my $tumor_sid=$ARGV[1];
my $normal_sid=$ARGV[2];
if(defined $tumor_sid){
	chomp $tumor_sid;
}
if(defined $normal_sid){
	chomp $normal_sid;
}
my $header_line_idx=0;

open(IN,"<$in_file") or die "can not open $in_file\n";
open(OUT_ADJUST,">temp_input.vcf") or die "can not open temp_input.vcf\n";
open(OUT_TUMOR,">tumor_for_snpEff.vcf") or die "can not open tumor_for_snpEff.vcf\n";
printf OUT_TUMOR "#CHROM	ID	Database_infor	REF	ALT	VarFreq_precent	Cov	Var_Cov	Genome_Type	Source\n";
if(defined $normal_sid){
	#print "yes\n";
	open(OUT_NORMAL,">normal_for_snpEff.vcf") or die "can not open normal_for_snpEff.vcf\n";
	open(OUT_COMPARE,">compare_for_snpEff.vcf") or die "can not open compare_for_snpEff.vcf\n";
	printf OUT_NORMAL "#CHROM	ID	Database_infor	REF	ALT	VarFreq_precent	Cov	Var_Cov	Genome_Type	Source\n";
	printf OUT_COMPARE "#CHROM	ID	Database_infor	REF	ALT	VarFreq_precent	Cov	Var_Cov	Genome_Type	Source\n";
}
#change the windows format to the linux format
while(<IN>){
	$_=~s/\r\n/\n/ if ($_=~/\r\n/);
	$_=~s/\r/\n/ if ($_=~/\r/);
	my $rec=$_;
	chomp $rec;
	$rec=~s/：/:/g;
	print OUT_ADJUST "$rec\n";
}
close IN;
close OUT_ADJUST;

#find the header line
my $header_index=0;

open(IN,"<temp_input.vcf") or die "can not open temp_input.vcf \n";
while(<IN>){
	#print "$_";
	chomp $_;
	$header_index++;
	if($_=~/#CHROM/){
		$header_line_idx=$header_index;
	}
}
close IN;
#print "header index:$header_line_idx\n";
#find the columns of the INFO FORMAT and tumor_sid, and justify is it a single sample or compare sample
my $numb=0;
my ($filter_idx,$tumor_idx,$normal_idx,$info_idx,$format_idx,$single_or_compare);
open(IN,"<temp_input.vcf") or die "can not open temp_input.vcf\n";
while(<IN>){
	$numb=$numb+1;
	if($numb==$header_line_idx){
		$_=~s/( )+//g;
		my @recs=split('\t',$_);
		if(scalar(@ARGV)==3){
			print "compare sample\n";
			if(($_=~/$tumor_sid/)&&($_=~/$normal_sid/)){
				for(my $i=0;$i<scalar(@recs);$i++){
					#print "$recs[$i]\n";
					if($recs[$i]=~/$tumor_sid/g){
						$tumor_idx=$i;
					}elsif($recs[$i] =~ /$normal_sid/){
						$normal_idx=$i;
					}elsif($recs[$i]=~/FILTER/){
						$filter_idx=$i;
					}
				}
			}else{
				die "ERROR: The name of the sample id is different with the iCMDB\n ";
			}
		}else{
			print "single sample\n";
			for(my $i=0;$i<scalar(@recs);$i++){
				#print "$recs[$i]\n";
				if($recs[$i]=~/FORMAT/){
					$tumor_idx=$i+1;
					#print "$tumor_idx\n";
				}elsif($recs[$i]=~/FILTER/){
					$filter_idx=$i;
				}
			}
		}
	}elsif($numb==$header_line_idx+1){
		$_=~s/( )+//g;
		my @recs=split('\t',$_);
		if(($_=~m/\=/)||($_=~m/\:/)||($_=~m/\：/)){
			for(my $i=0;$i<scalar(@recs);$i++){
				if($recs[$i]=~/\=/g){
					$info_idx=$i;
					#print "infor:$info_idx\n";
				}elsif(($recs[$i]=~/GT\:/)||($recs[$i]=~m/GT\：/)){
					$format_idx=$i;
					#print "format:$format_idx\n";
				}
			}
			if((defined $info_idx)&&(defined $format_idx)){
				if((scalar(@recs)-$format_idx-1)==1){
					$single_or_compare='single';
				}else{
					$single_or_compare='compare';
				}
			}elsif(defined $info_idx){
				if((scalar(@recs)-$info_idx-1)==1){
					$single_or_compare='single';
				}else{
					$single_or_compare='compare';
				}
			}elsif(defined $format_idx){
				if((scalar(@recs)-$format_idx-1)==1){
					$single_or_compare='single';
				}else{
					$single_or_compare='compare';
				}
			}
		}else{
			die "ERROR: can not find the column: INFO and FORMAT\n";
		}
	}

}
close IN;
#standard for GT DP AD and decompress the variant with more than 2 allele
$numb=0;
open(IN,"<temp_input.vcf") or die "can not open temp_input.vcf\n";
while(<IN>){
	chomp $_;
	next if (m/^\s*$/);
	$numb++;
	my ($filter,$dp_tumor,$ad_tumor);
	my ($gt_tumor,$varfreq_tumor,$ad_normal,$dp_normal,$gt_normal,$varfreq_normal,$id_tumor,$id_normal,$id_compare);
	my (@normal_formats,@alts,$alt,@ads_tumor,@ads_normal,@ads,);
	my ($varfreq,$dp,$ad,$gt,$dp_tag,$ad_tag);
	$_=~s/\r|\n//g;
	$_=~s/( )+//g;
	if($_!~m/^#/){
		my @recs=split('\t',$_);
		
		my @formats=split(':|：',$recs[$format_idx]);
		my @tumor_formats=split(':|：',$recs[$tumor_idx]);
		my (%tumor_format_rec,%tag_hash,%normal_format_rec)=();
		my @tags=split('\:',$recs[$format_idx]);
		for(my $i=0;$i<scalar(@tags);$i++){
			$tag_hash{$tags[$i]}=$tags[$i];
		}
		if((exists $tag_hash{'FDP'})&&((exists $tag_hash{'FAO'})||(exists $tag_hash{'FAD'}))){
			if(exists $tag_hash{'FAO'}){
				$ad_tag='FAO';
			}else{
				$ad_tag='FAD';
			}
			$dp_tag='FDP';
		}elsif((exists $tag_hash{'DP'})&&((exists $tag_hash{'AO'})||(exists $tag_hash{'AD'})||(exists $tag_hash{'CLCAD2'}))){
			if(exists $tag_hash{'AO'}){
				$ad_tag='AO';
			}elsif(exists $tag_hash{'AD'}){
				$ad_tag='AD';
			}else{
				$ad_tag='CLCAD2';
			}
			$dp_tag='DP';
		}else{
			$ad_tag='NA';
			$dp_tag='NA';
		}
		if(($ad_tag eq 'NA') &&($dp_tag eq 'NA')){
			if(exists $tag_hash{'AD'}){
				$ad_tag='AD';
			}
		}
		#print "$tumor_formats[1]\n";
		if((($filter_idx ne "")&&(($recs[$filter_idx] eq "PASS")||($recs[$filter_idx] eq '.')||($recs[$filter_idx] eq "")))||($filter_idx eq "")){
			$filter='yes';
		}else{
			$filter='no';
		}
		for(my $i=0;$i<scalar(@formats);$i++){#tumor format values
			$tumor_format_rec{$formats[$i]}=$tumor_formats[$i];
			# if(($recs[1] eq '15685950')){
				# print "$format_idx|$recs[1]|$formats[$i]|$tumor_formats[$i]\n";
			# }
		}
		
		#print"$dp_tag|$ad_tag\n";
		if($dp_tag eq 'NA'){
			$dp_tumor=0;
			if($ad_tag ne 'NA'){
				$ad_tumor=$tumor_format_rec{$ad_tag};
				#print"$ad_tumor\n";
				if($ad_tumor =~m/\,/){
					my @ad_tumors=split(',',$ad_tumor);
					for(my $i=0;$i<scalar(@ad_tumors);$i++){
						$dp_tumor=$dp_tumor+$ad_tumors[$i];
					}
				}
			}else{
				$ad_tumor='-1';
				$dp_tumor='-1';
			}
		}else{
			$dp_tumor=$tumor_format_rec{$dp_tag};
			$ad_tumor=$tumor_format_rec{$ad_tag};
		}
		#print "$dp_tag|$ad_tag|$dp_tumor|$ad_tumor\n";
		if($single_or_compare eq 'compare'){#normal format values
			@normal_formats=split(':|：',$recs[$normal_idx]);
			#print "$normal_idx\n";
			for(my $i=0;$i<scalar(@formats);$i++){
				$normal_format_rec{$formats[$i]}=$normal_formats[$i];
			}
			if($dp_tag eq 'NA'){
				$dp_normal=0;
				if($ad_tag ne 'NA'){
					$ad_normal=$normal_format_rec{$ad_tag};
					if($ad_normal =~m/\,/){
						my @ad_normals=split(',',$ad_normal);
						for(my $i=0;$i<scalar(@ad_normals);$i++){
							$dp_normal=$dp_normal+$ad_normals[$i];
						}
					}
				}else{
					$ad_normal='-1';
					$dp_normal='-1';
				}
				
			}else{
				$dp_normal=$normal_format_rec{$dp_tag};
				$ad_normal=$normal_format_rec{$ad_tag};
			}
		}
		if(($recs[4]=~m/\,/)||($recs[4]=~m/\//)){#more than 2 allele in one record
			my $ad_tumor_new=$ad_tumor;
			if($recs[4]=~m/\,/){
				@alts=split('\,',$recs[4]);
			}else{
				@alts=split('\/',$recs[4]);
			}
			if($ad_tumor ne '-1'){
				my @tumor_ads=split(',',$ad_tumor_new);
				if($dp_tumor==-1){
					$dp_tumor=0;
					if(((scalar(@tumor_ads))-(scalar(@alts)))==1){
						for(my $i=0;$i<scalar(@alts);$i++){
							$dp_tumor=$dp_tumor+$tumor_ads[$i];
						}
					}else{
						$dp_tumor=-1;
					}
				}
			}
			if($single_or_compare eq 'compare'){
				if($ad_normal!= -1){
					my @normal_ads=split(',',$ad_normal);
					if($dp_normal==-1){
						$dp_normal=0;
						if(((scalar(@normal_ads))-(scalar(@alts)))==1){
							for(my $i=0;$i<scalar(@alts);$i++){
								$dp_normal=$dp_normal+$normal_ads[$i];
							}
						}else{
							$dp_normal=-1;
						}
					}
				}
			}
			for(my $i=0;$i<scalar(@alts);$i++){
				$alt=$alts[$i];
				if($ad_tumor ne '-1'){
					@ads_tumor=split(',',$ad_tumor_new);
					if(scalar(@ads_tumor)==scalar(@alts)){
						$ad_tumor=$ads_tumor[$i];
					}elsif(((scalar(@ads_tumor))-(scalar(@alts)))==1){
						$ad_tumor=$ads_tumor[$i+1];
					}
				}
				if(($dp_tumor eq '-1')||($ad_tumor eq '-1')){
					$dp_tumor=-1;
					$ad_tumor=-1;
					$varfreq_tumor=-1;
				}else{
					if($dp_tumor==0){
						$varfreq_tumor=0;
					}else{
						$varfreq_tumor=sprintf("%.2f",($ad_tumor/$dp_tumor*100));
					}
				}
				if(($tumor_format_rec{'GT'} ne '0/0')&&($tumor_format_rec{'GT'} ne '0/1')&&($tumor_format_rec{'GT'} ne '1/1')){#it is './.' or '1/2' etc
					if($dp_tumor ne '-1'){
						if($varfreq_tumor >90){
							$gt_tumor='1/1';
						}elsif($varfreq_tumor <20){
							$gt_tumor='0/0';
						}else{
							$gt_tumor='0/1';
						}
					}else{
						$gt_tumor='NA';
					}
				}else{
					$gt_tumor=$tumor_format_rec{'GT'};
				}
				my $tumor_ref_length=length($recs[3]);
				my $tumor_alt_length=length($recs[4]);
				my $tumor_ref=$recs[3];
				my $tumor_alt=$alt;
				if($tumor_ref_length>50){
					$tumor_ref=$tumor_ref_length;
				}
				if($tumor_alt_length>50){
					$tumor_alt=$tumor_ref_length;
				}
				if($gt_tumor eq '0/0'){
					my $ref_str=substr($recs[3],0,1);
					$gt_tumor=join("\/",($ref_str,$ref_str));
				}else{
					$gt_tumor=~s/0/$tumor_ref/g;
					$gt_tumor=~s/1/$tumor_alt/g;
				}
				if($single_or_compare eq 'compare'){#there are two samples in the vcf file
					if($ad_normal ne '-1'){
						@ads_normal=split(',',$ad_normal);
						if(scalar(@ads_normal)==scalar(@alts)){
							$ad_normal=$ads_normal[$i];
						}elsif(((scalar(@ads_normal))-(scalar(@alts)))==1){
							$ad_normal=$ads_normal[$i+1];
						}
					}
					if(($dp_normal eq '-1')||($ad_normal eq '-1')){
						$dp_normal=-1;
						$ad_normal=-1;
						$varfreq_normal=-1;
					}else{
						if($dp_normal==0){
							$varfreq_normal=0;
						}else{
							$varfreq_normal=sprintf("%.2f",($ad_normal/$dp_normal*100));
						}
					}
					if(($normal_format_rec{GT} ne '0/0')&&($normal_format_rec{GT} ne '0/1')&&($normal_format_rec{GT} ne '1/1')){#it is './.' or '1/2' etc
						if($dp_normal ne '-1'){
							if($varfreq_normal >90){
								$gt_normal='1/1';
							}elsif($varfreq_normal <20){
								$gt_normal='0/0';
							}else{
								$gt_normal='0/1';
							}
						}else{
							$gt_normal='NA';
						}
					}else{
						$gt_normal=$normal_format_rec{GT};
					}
					$gt_normal=~s/0/$tumor_ref/g;
					$gt_normal=~s/1/$tumor_alt/g;
					$varfreq=join("_",($varfreq_tumor,$varfreq_normal));
					$dp=join("_",($dp_tumor,$dp_normal));
					$ad=join("_",($ad_tumor,$ad_normal));
					$gt=join("_",($gt_tumor,$gt_normal));
				}
				if($alt eq ""){
					$alt=$recs[3];
				}
				
				if((($recs[3] ne $alt)&&($ad_tumor!=0)&&($filter eq 'yes'))||(($recs[3] eq $alt)&&($filter eq 'yes'))){
				#if((($recs[3] ne $alt)&&($ad_tumor!=0))||(($recs[3] eq $alt))){
					$recs[0]=~s/chr//g;
					$id_tumor=join("\t",(@recs[0..3],$alt,$varfreq_tumor,$dp_tumor,$ad_tumor,$gt_tumor,'Unknown'));
					printf OUT_TUMOR "$id_tumor\n";
					if($single_or_compare eq 'compare'){
						$id_normal=join("\t",(@recs[0..3],$alt,$varfreq_normal,$dp_normal,$ad_normal,$gt_normal,'Unknown'));
						$id_compare=join("\t",(@recs[0..3],$alt,$varfreq,$dp,$ad,$gt,'Unknown'));
						if((($recs[3] ne $alt)&&($ad_normal!=0))||($recs[3] eq $alt)){
							printf OUT_NORMAL "$id_normal\n";
						}
						if($ad_tumor != 0){
							printf OUT_COMPARE "$id_compare\n";
						}
					}
				}
			}
		}else{ #only two allele in one record
			$alt=$recs[4];
			if($ad_tumor=~m/\,/){
				@ads_tumor=split(',',$ad_tumor);
				$ad_tumor=$ads_tumor[1];
				if($dp_tumor ne '-1'){
					$dp_tumor=$ads_tumor[0]+$ads_tumor[1];
				}
			}
			if(($dp_tumor eq '-1')||($ad_tumor eq '-1')){
				$dp_tumor=-1;
				$ad_tumor=-1;
				$varfreq_tumor=-1;
			}else{
				if($dp_tumor==0){
					$varfreq_tumor=0;
				}else{
					$varfreq_tumor=sprintf("%.2f",($ad_tumor/$dp_tumor*100));
				}
			}
			#print "$recs[1]|GT|$tumor_format_rec{'GT'}\n";
			if(($tumor_format_rec{'GT'} ne '0/0')&&($tumor_format_rec{'GT'} ne '0/1')&&($tumor_format_rec{'GT'} ne '1/1')){#it is './.' or '1/2' etc
				if($varfreq_tumor!=0){
                    $varfreq_tumor=sprintf("%.2f",($ad_tumor/$dp_tumor*100));
                }
				if($dp_tumor ne '-1'){
					if($varfreq_tumor >90){
						$gt_tumor='1/1';
					}elsif($varfreq_tumor <20){
						$gt_tumor='0/0';
					}else{
						$gt_tumor='0/1';
					}
				}else{
					$gt_tumor='NA';
				}
			}else{
				$gt_tumor=$tumor_format_rec{GT};
			}
			my $tumor_ref_length=length($recs[3]);
			my $tumor_alt_length=length($recs[4]);
			my $tumor_ref=$recs[3];
			my $tumor_alt=$alt;
			if($tumor_ref_length>50){
				$tumor_ref=$tumor_ref_length;
			}
			if($tumor_alt_length>50){
				$tumor_alt=$tumor_ref_length;
			}
			if($gt_tumor eq '0/0'){
				my $ref_str=substr($recs[3],0,1);
				$gt_tumor=join("\/",($ref_str,$ref_str));
			}else{
				$gt_tumor=~s/0/$tumor_ref/g;
				$gt_tumor=~s/1/$tumor_alt/g;
			}
			if($single_or_compare eq 'compare'){
				$alt=$recs[4];
				if($ad_normal=~/,/g){
					@ads_normal=split(',',$ad_normal);
					$ad_normal=$ads_normal[1];
				}
				if(($dp_normal eq '-1')||($ad_normal eq '-1')){
					$dp_normal=-1;
					$ad_normal=-1;
					$varfreq_normal=-1;
				}else{
					if($dp_normal==0){
						$varfreq_normal=0;
					}else{
						$varfreq_normal=sprintf("%.2f",($ad_normal/$dp_normal*100));
					}
				}
				if(($normal_format_rec{'GT'} ne '0/0')&&($normal_format_rec{'GT'} ne '0/1')&&($normal_format_rec{'GT'} ne '1/1')){#it is './.' or '1/2' etc
					#$varfreq_normal=$ad_normal/$dp_normal;
					if($dp_normal!= -1){
						if($varfreq_normal >90){
							$gt_normal='1/1';
						}elsif($varfreq_normal <20){
							$gt_normal='0/0';
						}else{
							$gt_normal='0/1';
						}
					}else{
						$gt_normal='NA';
					}
				}else{
					$gt_normal=$normal_format_rec{'GT'};
				}
				if($gt_normal eq '0/0'){
					my $ref_str=substr($recs[3],0,1);
					$gt_normal=join("\/",($ref_str,$ref_str));
				}else{
					$gt_normal=~s/0/$tumor_ref/g;
					$gt_normal=~s/1/$tumor_alt/g;
				}
				$varfreq=join("_",($varfreq_tumor,$varfreq_normal));
				$dp=join("_",($dp_tumor,$dp_normal));
				$ad=join("_",($ad_tumor,$ad_normal));
				$gt=join("_",($gt_tumor,$gt_normal));
			}
			# if($recs[1] eq "96741053"){
				# $alt=$recs[3];
				# print "alt:$alt\n$_\n";
			# }
			if($alt eq ""){
				$alt=$recs[3];
			}
			if((($recs[3] ne $alt)&&($ad_tumor!=0)&&($filter eq 'yes'))||(($recs[3] eq $alt)&&($filter eq 'yes'))){
			#if((($recs[3] ne $alt)&&($ad_tumor!=0))||(($recs[3] eq $alt))){
				$recs[0]=~s/chr//g;
				$id_tumor=join("\t",(@recs[0..3],$alt,$varfreq_tumor,$dp_tumor,$ad_tumor,$gt_tumor,'Unknown'));
				printf OUT_TUMOR "$id_tumor\n";
				if($single_or_compare eq 'compare'){
					$id_normal=join("\t",(@recs[0..3],$alt,$varfreq_normal,$dp_normal,$ad_normal,$gt_normal,'Unknown'));
					$id_compare=join("\t",(@recs[0..3],$alt,$varfreq,$dp,$ad,$gt,'Unknown'));
					if((($recs[3] ne $alt)&&($ad_normal!=0))||($recs[3] eq $alt)){
						printf OUT_NORMAL "$id_normal\n";
					}
					if($ad_tumor != 0){
						printf OUT_COMPARE "$id_compare\n";
					}
				}
			}
		}
	}
}
close IN;
close OUT_TUMOR;
if($single_or_compare eq 'compare'){
	close OUT_NORMAL;
	close OUT_COMPARE;
}