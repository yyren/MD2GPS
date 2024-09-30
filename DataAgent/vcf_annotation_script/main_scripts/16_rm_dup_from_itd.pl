#/usr/bin/perl
use strict;
my $in_merge_file=$ARGV[0];
my $in_itd_file=$ARGV[1];
my $out_merge_rm_dup=$ARGV[2];

my @headers;
my @recs;
my %rec_id=();
my $id;
my %itd_rec=();
open(IN, "<$in_itd_file") or die "can not open $in_itd_file\n";
open(OUT, ">$out_merge_rm_dup") or die "can not open $out_merge_rm_dup\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	
	if($_=~m/CHROM/g){
		@headers=split('\t',$_);
	}else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$rec_id{$headers[$i]}=$recs[$i];
		}
		$id=join("-",($rec_id{"CHROM"}, $rec_id{"ID"}, $rec_id{"REF"},$rec_id{"ALT"}));
		$itd_rec{$id}=$_;
	}
}
close IN;
open(IN, "<$in_merge_file") or die "can not open $in_merge_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	
	if($_=~m/CHROM/g){
		@headers=split('\t',$_);
		printf OUT "$_\n";
	}else{
		@recs=split('\t', $_);
		for(my $i=0;$i<scalar(@headers);$i++){
			$rec_id{$headers[$i]}=$recs[$i];
		}
		$id=join("-",($rec_id{"CHROM"}, $rec_id{"ID"}, $rec_id{"REF"},$rec_id{"ALT"}));
		if(! exists $itd_rec{$id}){
			printf OUT "$_\n";
		}
	}
}
close IN;
close OUT;