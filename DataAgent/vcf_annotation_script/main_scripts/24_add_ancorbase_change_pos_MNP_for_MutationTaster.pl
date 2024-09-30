use strict;
use warnings;
use Bio::DB::Fasta;



my $in_reference=$ARGV[0];
my $in_file=$ARGV[1];
my $out_file=$ARGV[2];

my $db = Bio::DB::Fasta->new($in_reference);

open(IN,"<$in_file") or die "can not open $in_file\n";
open(OUT,">$out_file") or die "can not open $out_file\n";
my ($chr_idx, $pos_idx, $ref_idx, $alt_idx, $gt_idx);
while(my $line=<IN>){
    chomp $line;
    my @recs=split(/\t/,$line);
    if($line=~m/CHROM/){
        print OUT "$line\n";
        for(my $i=0;$i<=$#recs;$i++){
            if($recs[$i]=~m/CHROM/i){
                $chr_idx=$i;
            }elsif($recs[$i]=~m/^(ID|POS)$/i){
                $pos_idx=$i;
                #print "pos_idx: $pos_idx\n";
            }elsif($recs[$i]=~m/^REF$/i){
                $ref_idx=$i;
            }elsif($recs[$i]=~m/^ALT$/i){
                $alt_idx=$i;
            }elsif($recs[$i]=~m/^Genome\_Type$/i){
                $gt_idx=$i;
            }
        }
    }else{
        if((length($recs[$ref_idx]) == length($recs[$alt_idx])) and (length($recs[$ref_idx])>=2)){#mnp variant
            #print "1: $recs[$chr_idx]|$recs[$pos_idx]|$recs[$ref_idx]|$recs[$alt_idx]\n";
            
            $recs[$pos_idx]=$recs[$pos_idx]-1;
            my $bed_chr=$recs[$chr_idx];
            $bed_chr=~s/chr//;
            my $ancor_base=$db->seq($bed_chr, $recs[$pos_idx] => $recs[$pos_idx]);
            my $new_ref=$ancor_base.$recs[$ref_idx];
            my $new_alt=$ancor_base.$recs[$alt_idx];
            my @gts=split(/\//,$recs[$gt_idx]);
            my $gt=$new_ref.'/'.$new_alt;
            if($gts[0] eq $gts[1]){
                if($gts[0] eq $recs[$ref_idx]){
                    $gt=$new_ref.'/'.$new_ref;
                }else{
                    $gt=$new_alt.'/'.$new_alt;
                }
            }
            $recs[$ref_idx]=$new_ref;
            $recs[$alt_idx]=$new_alt;
            $recs[$gt_idx]=$gt;
            #print "2: $recs[$chr_idx]|$recs[$pos_idx]|$recs[$ref_idx]|$recs[$alt_idx]\n";
        }
        my $new_line=join("\t",@recs);
        print OUT "$new_line\n";
    }
}
close IN;
close OUT;
