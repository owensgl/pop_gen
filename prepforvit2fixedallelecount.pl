#!/bin/perl
use strict;
use warnings;

my $cut_off = 0.98; #minimum difference in allele frequency between parents
my %name;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (5..$#a){
      $name{$i} = $a[$i];
    }
    print "name\tchr\tbp\tcm\tpet_alleles";
  }else{
    my $chr = $a[0];
    my $bp = $a[1];
    my $cm = $a[2];
    my $p1 = $a[3];
    my $p2 = $a[4];
    my $dif = abs($p1 - $p2);
    if ($dif < $cut_off){
      next;
    }
    my $state = "P2";
    if ($p1 > $p2){
      $state = "P1";
    }
    foreach my $i (5..$#a){
      if ($a[$i] eq "N"){
        next;
      }
      if (($a[$i] eq "0") and ($state eq "P1")){
        print "\n$name{$i}\t$chr\t$bp\t$cm\t2";
      }elsif (($a[$i] eq "0") and ($state eq "P2")){
        print "\n$name{$i}\t$chr\t$bp\t$cm\t0";
      }elsif (($a[$i] eq "2") and ($state eq "P1")){ 
         print "\n$name{$i}\t$chr\t$bp\t$cm\t0";
      }elsif (($a[$i] eq "2") and ($state eq "P2")){ 
        print "\n$name{$i}\t$chr\t$bp\t$cm\t2";
      }elsif ($a[$i] eq "1"){
        print "\n$name{$i}\t$chr\t$bp\t$cm\t1";
      }
    }
  }
}
