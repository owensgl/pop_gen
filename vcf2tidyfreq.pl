#!/bin/perl
use warnings;
use strict;
#Calculates allele frequency and outputs in a tidy format. Only works for biallelic snps.
my $popfile = $ARGV[0]; #Two columns, column 1 = name, column 2 = population. No header;

open POP, $popfile;

my %pop;
my %poplist;
while(<POP>){
  chomp;
  my @a = split(' ',$_);
  my $name = $a[0];
  my $pop = $a[1];
  $pop{$name} = $pop;
  $poplist{$pop}++;
}
close POP;

print "chr\tpos\tgroup\tn\talt_freq";
my %sample;
while(<STDIN>){
  my %altcount;
  my %totalcount;
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/##/){next;}
  if ($_ =~ m/#/){
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $pos = $a[1];
    foreach my $i (9..$#a){
      unless ($pop{$sample{$i}}){next;} #Skip non-population assigned samples
      my @info = split(/:/,$a[$i]);
      if (($info[0] eq '.') or ($info[0] eq './.')){next;}
      my @bases = split("\/",$info[0]);
      foreach my $j (0..1){
        if ($bases[$j] >= 2){
          die "Non-biallelic bases found in $chr-$pos\n";
        }
        $altcount{$pop{$sample{$i}}}+=$bases[$j];
        $totalcount{$pop{$sample{$i}}}++;
      }
    }
    foreach my $pop (sort keys %totalcount){
      unless($altcount{$pop}){
        $altcount{$pop} = 0;
      }
      my $freq = $altcount{$pop}/$totalcount{$pop};
      print "\n$chr\t$pos\t$pop\t$totalcount{$pop}\t$freq";
    }
  }
}
