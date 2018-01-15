#!/bin/perl
use strict;
use warnings;

my $popfile = $ARGV[0];

my %info;
open POP, $popfile;
my %groups;
while(<POP>){
  chomp;
  my @a = split(' ',$_);
  $info{$a[0]} = $a[1];
  $groups{$a[1]}++;
}
close POP;

my %gene;
my %totalcounts;
my %counts;
while(<STDIN>){
  chomp;
  my @a = split(' ',$_);
  if ($. == 1){
    foreach my $i (1..$#a){
      $gene{$i} = $a[$i];
    }
  }else{
    my $sample = $a[0];
    unless ($info{$sample}){
#print STDERR "Skipping $sample\n";
      next;
     }
    foreach my $i  (1..$#a){
      $counts{$info{$sample}}{$gene{$i}}+=$a[$i];
      $totalcounts{$info{$sample}}{$gene{$i}}++;
    }
  } 
}
print "gene";
foreach my $group (sort keys %groups){
  print "\t$group";
}
foreach my $gene (sort values %gene){
  print "\n$gene";
  foreach my $group (sort keys %groups){
    unless ($counts{$group}{$gene}){
      $counts{$group}{$gene} = 0;
    }
    my $percent_present = $counts{$group}{$gene}/$totalcounts{$group}{$gene};
    print "\t$percent_present";
  }
}
