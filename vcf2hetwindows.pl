#!/bin/perl 
use warnings;
use strict;
use POSIX;

#This script takes a vcf and calculates observed heterozygosity per window. 

my %sample;
my $window_size = 100000;
my %het;
my %site_count;
my $min_depth = 5;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){
    next; 
  }
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    next;
  }
  
  my $chr = $a[0];
  my $pos = $a[1];
  my $window = floor($pos/$window_size)*$window_size;
  foreach my $i (9..$#a){
    my @fields =split(/:/,$a[$i]);
    my $dp = $fields[2];
    my $genotype = $fields[0];
    if (($genotype eq './.') or ($genotype eq '.')){next;}
    if ($dp < $min_depth){next;}
    if ($genotype eq '0/1'){
      $het{$sample{$i}}{$chr}{$window}++;
      $site_count{$sample{$i}}{$chr}{$window}++;
    }else{
      $site_count{$sample{$i}}{$chr}{$window}++;
    }
  }
}
print "name\tchr\tstart\tend\thets\tcounts\tperc_het";

foreach my $sample (sort keys %site_count){
  foreach my $chr (sort keys %{$site_count{$sample}}){
    foreach my $window (sort {$a <=> $b} keys %{$site_count{$sample}{$chr}}){
      my $window_end = $window+ $window_size;
      unless($het{$sample}{$chr}{$window}){
        $het{$sample}{$chr}{$window} = 0;
      }
      my $perc_het = $het{$sample}{$chr}{$window}/$site_count{$sample}{$chr}{$window};
      print "\n$sample\t$chr\t$window\t$window_end\t$het{$sample}{$chr}{$window}\t";
      print "$site_count{$sample}{$chr}{$window}\t$perc_het";
    }
  }
}
