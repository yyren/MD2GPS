#!/usr/bin/perl
use strict;

my $in_normal=$ARGV[0];
my $in_tumor=$ARGV[1];
my $out_somatic=$ARGV[2];

my %normal_rec; 
open(IN,"<$in_normal") or die "can not open $in_normal:$!";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//;
	if($_ !~ m/#CHROM/){
		my @recs=split('\t',$_);
		my $id_rec=join("_",@recs[0,1,3,4]);
		$normal_rec{$id_rec}=$_;
	}
}
close IN;

my ($tumor_query,$normal_query,$normal_recs);
open(IN,"<$in_tumor") or die "can not open $in_tumor:$!";
open(OUT,">$out_somatic") or die "can not open $out_somatic:$!";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//;
	if($_=~m/#CHROM/){
		printf OUT "$_\n";
	}else{
		my @recs=split('\t',$_);
		my $id_rec=join("_",@recs[0,1,3,4]);
		if(exists $normal_rec{$id_rec}){
			my @normal_recs=split('\t',$normal_rec{$id_rec});
			my $normal_ad=$normal_recs[7];
			my $normal_gt=$normal_recs[8];
			if(($recs[7]>=10)&&($normal_ad>=10)&&($recs[8] eq $normal_gt)){# support varant reads>10(both case and control sample) and with the same genotype
				#it is germline
				#print "$_\n";
			}else{
				my $varfreq=join("_",($recs[5],$normal_recs[5]));
				my $dp=join("_",($recs[6],$normal_recs[6]));
				my $ad=join("_",($recs[7],$normal_recs[7]));
				my $gt=join("_",($recs[8],$normal_recs[8]));
				my $id_rec=join("\t",(@recs[0..4],$varfreq,$dp,$ad,$gt,$recs[9]));
				printf OUT "$id_rec\n";
				#print "ttt:$id_rec\n";
			}
		}else{
			my $varfreq=join("_",($recs[5],'0'));
			my $dp=join("_",($recs[6],'0'));
			my $ad=join("_",($recs[7],'0'));
			my $normal_ref=substr($recs[3],0,1);
			my $gt_normal=join("\/",($normal_ref,$normal_ref));
			my $gt=join("_",($recs[8],$gt_normal));
			my $id_rec=join("\t",(@recs[0..4],$varfreq,$dp,$ad,$gt,$recs[9]));
			
			printf OUT "$id_rec\n";
		}
	}
}
close IN;
close OUT;
