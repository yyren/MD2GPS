#!/usr/bin/perl
my $in_file=$ARGV[0];#bed file
my $in_file2=$ARGV[1];#snpeff input file
my $raw_occur_left_right=$ARGV[2];#Illumina or Ion_torrent
my $reference=$ARGV[3];#fasta file with path
my $out_file=$ARGV[4];#out put file with two additional column:strand and ref seq(200 bases)

my (@recs,$id,%rec_id,$id2,$ref_length,$alt_length,$start,$end,$bed_id,$seq);
open(IN,"<$in_file") or die "can not open $in_file\n";
open(OUT,">$out_file") or die "can not open $out_file\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	@recs=split('\t',$_);
	for(my $i=$recs[1];$i<=$recs[2];$i++){
		$id=join("_",($recs[0],$i));
		$rec_id{$id}=$recs[4];
		#print "$id\n";
	}
}
close IN;
open(IN,"<$in_file2") or die "can not open $in_file2\n";
while(<IN>){
	chomp $_;
	$_=~s/\r|\n//g;
	@recs=split('\t',$_);
	if($_=~m/#/g){
		printf OUT "$_\tstrand\tRefseq\n";
	}else{
		$id2=join("_",($recs[0],$recs[1]));
		$ref_length=length($recs[3]);
		$alt_length=length($recs[4]);
		if(exists $rec_id{$id2}){
			if((($raw_occur_left_right eq 'left')&&($rec_id{$id2} eq '-'))||(($raw_occur_left_right eq 'right')&&($rec_id{$id2} eq '+'))){
				printf OUT "$_\t$rec_id{$id2}\tpass\n";
			}elsif((($ref_length==1)&&($alt_length !=1))||(($ref_length!=1)&&($alt_length==1))){
				if(($raw_occur_left_right eq 'left')&&($rec_id{$id2} eq '+')){
					$start=$recs[1]-1;
					$end=$start+200;
					$bed_id=join("\t",($recs[0],$start,$end));
					`echo '$bed_id' >temp.bed`;
					`fastaFromBed -fi $reference -bed temp.bed -fo temp_seq.txt -tab`;
					$seq=`sed -n "1, 1p" temp_seq.txt |awk '{print \$2}'`;
					chomp($seq);
					printf OUT "$_\t$rec_id{$id2}\t$seq\n";
				}elsif(($raw_occur_left_right eq 'right')&&($rec_id{$id2} eq '-')){
					$end=$recs[1]+$ref_length-1;
					$start=$end-200;
					$bed_id=join("\t",($recs[0],$start,$end));
					`echo '$bed_id' >temp.bed`;
					`fastaFromBed -fi $reference -bed temp.bed -fo temp_seq.txt -tab`;
					$seq=`sed -n "1, 1p" temp_seq.txt |awk '{print \$2}'`;
					chomp($seq);
					printf OUT "$_\t$rec_id{$id2}\t$seq\n";
				}else{
					printf OUT "$_\t$rec_id{$id2}\tpass\n";
					#print "$rec_id{$id2}|$_\n";
				}
			}else{
				printf OUT "$_\t$rec_id{$id2}\tpass\n";
			}
		}else{
			printf OUT "$_\tunknown\tpass\n";
		}
	}
}
close IN;
close OUT;