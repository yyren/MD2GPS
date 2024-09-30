use strict;
use warnings;

my $in_hg19=$ARGV[0];
my $in_hg38_given_region=$ARGV[1];
my $in_hg38=$ARGV[2];
my $out_file=$ARGV[3];

my %hg38_given_region_hg19=();
my %hg38_given_region_hg38=();
my %hg38_recs_hg19=();

open(IN,"<$in_hg38_given_region") or die "cannot open $in_hg38_given_region\n";
open(OUT,">$out_file") or die "cannot open $out_file\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/^#/){
        print OUT "$line\n";
    }else{
        my @recs=split(/\t/,$line);
		my $id_hg19=join("\_",@recs[0,2]);
		my $id_hg38=join("\_",@recs[0,1]);
		if($line!~m/WARNING_REF_DOES_NOT_MATCH_GENOME/){
			$hg38_given_region_hg19{$id_hg19}=1;
			$hg38_given_region_hg38{$id_hg38}=1;
			$recs[1]=$recs[2];
			$recs[2]='.';
			my $new_line=join("\t",@recs);
			print OUT "$new_line\n";
		}
    }
}
close IN;
close OUT;


open(IN,"<$in_hg38") or die "cannot open $in_hg38\n";
open(OUT,">>$out_file") or die "cannot open $out_file\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/^#/){
        #print OUT "$line\n";
    }else{
        my @recs=split(/\t/,$line);
        my $id_hg19=join("\_",@recs[0,2]);
		my $id_hg38=join("\_",@recs[0,1]);
		$recs[1]=$recs[2];
        $recs[2]='.';
		if((not exists $hg38_given_region_hg38{$id_hg38}) and ($line!~m/WARNING_REF_DOES_NOT_MATCH_GENOME/)){
			$hg38_recs_hg19{$id_hg19}=1;
			my $new_line=join("\t",@recs);
			print OUT "$new_line\n";
		}
    }
}
close IN;
close OUT;

open(IN,"<$in_hg19") or die "cannot open $in_hg19\n";
open(OUT,">>$out_file") or die "cannot open $out_file\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/^#/){
		#do nothing
	}else{
		my @recs=split(/\t/,$line);
		my $id=join("\_",@recs[0,1]);
		if((not exists $hg38_recs_hg19{$id}) and (not exists $hg38_given_region_hg19{$id})){
			print OUT "$line\n";
		}
	}
}
close IN;
close OUT;
