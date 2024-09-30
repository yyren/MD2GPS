#!/usr/bin/perl -w
use strict;
my $in_target=$ARGV[0];
my $in_rd=$ARGV[1];
my $out_rd=$ARGV[2];
my $chr;
my $loci;
my $id;
my %id_target=();
my $start;
my $end;

open(IN_target, "<", $in_target) or die "can not open < $in_target: $!";
open(IN, "<", $in_rd) or die "can not open < $in_rd: $!";
open(OUT, ">", $out_rd) or die "cannot open < $out_rd: $!";

while (<IN_target>) {
    chomp $_;
    $_=~s/\r|\n//g;
    @_=split('\t', $_);
    if ($_[0] eq "TargetA") {
        #code
    }else{
        $chr=$_[0];
        #$chr=~s/chr//g;
        $start=$_[1];
        $end=$_[2];
        my $i=0;
        for($i=$start; $i<$end; $i++){
            $loci=$i;
            $id=join("\_", ($chr, $loci));
            if (exists $id_target{$id}) {
                $id_target{$id}=$id_target{$id}+1;
            }
            else{
                $id_target{$id}=1;
            }
            
        }
    }    
}
close IN_target;

while (<IN>) {
    chomp $_;
    $_=~s/\r|\n//g;
    @_=split('\t', $_);
    $chr=$_[0];
    $loci=$_[1];
    $id=join("\_", ($chr, $loci));
    if (exists $id_target{$id}) {
        printf OUT "$chr\t$loci\t$_[3]\t$id_target{$id}\n";
    }
    else{
        printf OUT "$chr\t$loci\t$_[3]\t0\n";
    }
}
close IN_target;
close OUT;