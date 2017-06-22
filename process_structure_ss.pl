#!/bin/perl
use warnings;
use strict;
#This script takes a results_ss file from structure sitebysite and sums columns to make the maternal and paternal values. It then averages them.
#This requires using the sitebysite= 1, linkage=1, phased=0 and markovphased=0.

print "sample\tloci\tM1\tM2\tM3\tM4\tM5\tP1\tP2\tP3\tP4\tP5\tT1\tT2\tT3\tT4\tT5";
while(<STDIN>){
    chomp;
    s/^\s+//;  # remove leading whitespace
    s/\s+$//; # remove trailing whitespace
    next unless length; # next rec unless anything left
    my @a = split(/ /, $_);
    unless ($#a == 26){ next;}
    my $sample = $a[0];
    unless (($sample == 60) or ($sample == 12)){ #For filter out individuals that were reference and don't vary.
        next;
    }
    my $loci = $a[1];
    my $M1 = ($a[2] + $a[3] + $a[4] + $a[5] + $a[6]);
    my $M2 = ($a[7] + $a[8] + $a[9] + $a[10] + $a[11]);
    my $M3 = ($a[12] + $a[13] + $a[14] + $a[15] + $a[16]);
    my $M4 = ($a[17] + $a[18] + $a[19] + $a[20] + $a[21]);
    my $M5 = ($a[22] + $a[23] + $a[24] + $a[25] + $a[26]);
    my $P1 = ($a[2] + $a[7] + $a[12] + $a[17] + $a[22]);
    my $P2 = ($a[3] + $a[8] + $a[13] + $a[18] + $a[23]);
    my $P3 = ($a[4] + $a[9] + $a[14] + $a[19] + $a[24]);
    my $P4 = ($a[5] + $a[10] + $a[15] + $a[20] + $a[25]);
    my $P5 = ($a[6] + $a[11] + $a[16] + $a[21] + $a[26]);
    my $T1 = ($M1 + $P1) / 2;
    my $T2 = ($M2 + $P2) / 2;
    my $T3 = ($M3 + $P3) / 2;
    my $T4 = ($M4 + $P4) / 2;
    my $T5 = ($M5 + $P5) / 2;
    print "\n$sample\t$loci\t$M1\t$M2\t$M3\t$M4\t$M5\t$P1\t$P2\t$P3\t$P4\t$P5\t$T1\t$T2\t$T3\t$T4\t$T5";

}
