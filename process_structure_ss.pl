#!/bin/perl
use warnings;
use strict;
#This script takes a results_ss file from structure sitebysite and sums columns to make the maternal and paternal values. It then averages them.
#This requires using the sitebysite= 1, linkage=1, phased=0 and markovphased=0.

my $k;
while(<STDIN>){
    chomp;
    s/^\s+//;  # remove leading whitespace
    s/\s+$//; # remove trailing whitespace
    next unless length; # next rec unless anything left
    my @a = split(/ /, $_);
    if($k){
      if ((sqrt($#a-1)) != $k){
	print STDERR "Wrong number of columns!";
      }
    }else{      
      $k = (sqrt($#a-1));
      print  "sample\tloci";
      foreach my $i (1..$k){
        print "\tT$i";
      }
    }
    my $sample = $a[0];
    my $loci = $a[1];
    my %m;
    my %p;
    my %t;
    foreach my $i (0..($k-1)){
      foreach my $j (0..($k-1)){
        $m{$i+1} += $a[($i * $k) + ($j) + 2];
        $p{$i+1} += $a[($i) + ($j * $k) + 2];
      }
    }
    foreach my $i (1..$k){
      $t{$i} = ($m{$i} + $p{$i}) / 2;
    }
    print "\n$sample\t$loci";
    foreach my $i (1..$k){
      print "\t$t{$i}";
    }
}
