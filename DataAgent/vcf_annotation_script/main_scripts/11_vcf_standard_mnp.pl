use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];


my (%pos, %ref, %alt, %gt, %remove_rec)=();
open(IN,"<$infile") or die "can not open $infile\n";
open(OUT,">$outfile") or die "can not open $outfile\n";
while(my $line=<IN>){
    chomp $line;
    if($line!~m/#/){
        my @recs=split(/\t/,$line);
        my $id=join("\_",@recs[0..1]);
        $pos{$id}=$recs[1];
        $ref{$id}=$recs[3];
        $alt{$id}=$recs[4];
        $gt{$id}=$recs[8];
    }
}
close IN;

open(IN,"<$infile") or die "can not open $infile\n";
while(my $line=<IN>){
    chomp $line;
    if($line=~m/#/){
        print OUT "$line\n";
    }else{
        my @recs=split(/\t/,$line);
        my $id=join("\_",@recs[0..1]);
        my $alt_pos=length($recs[4])+$recs[1]-1;
        my $query_start=$alt_pos+1;
        my $query_id=join("\_",($recs[0],$query_start));
        if((length($recs[3])==length($recs[4])) and (exists $pos{$query_id})){
            if((length($ref{$query_id})==length($alt{$query_id}))){#change 1 A T A/T; 2 C G C/G; -> 1 AC TG AC/TG, and combine into MNP only if the variants with the same genotype
                my @now_gt=split(/\//,$recs[8]);
                my @next_gt=split(/\//,$gt{$query_id});
                my $gt_check_tag='diff';
                if($now_gt[0] eq $now_gt[1]){
                    if($next_gt[0] eq $next_gt[1]){
                        $gt_check_tag='cosistent_homo';
                    }
                }else{
                    if($next_gt[0] ne $next_gt[1]){
                        $gt_check_tag='cosistent_het';
                    }
                }
                if($gt_check_tag ne 'diff'){# only combine the variants with the same genotype
                    my $new_ref=$recs[3].$ref{$query_id};
                    my $new_alt=$recs[4].$alt{$query_id};
                    if($gt_check_tag eq 'cosistent_homo'){
                        if($now_gt[0] eq $recs[3]){
                            $recs[8]=$new_ref.'/'.$new_ref;
                        }else{
                            $recs[8]=$new_alt.'/'.$new_alt;
                        }
                    }else{#cosistent_het
                        $recs[8]=$new_ref.'/'.$new_alt;
                    }
                    $recs[3]=$new_ref;
                    $recs[4]=$new_alt;
                    $remove_rec{$query_id}=$query_id;
                }
            }
        }
        my $new_line=join("\t",@recs);
        if(not exists($remove_rec{$id})){
            print OUT "$new_line\n";
        }
        
    }
}
close IN;
close OUT;