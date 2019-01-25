#!/bin/perl
use warnings;
use strict;

#This takes the vcf file, and group file and outputs the observed heterozygosity per site per group

my $popfile = $ARGV[0];
my $min_depth = 5;

open POP, $popfile;
my %pop;
while(<POP>){
  chomp;
  my @a = split(/\t/,$_);
  my $sample = $a[0];
  my $pop = $a[1];
  $pop{$sample} = $pop;
}

close POP;

my %sample;
my %het;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    next;
  }
  my $chr = $a[0];
  my $pos = $a[1];
  foreach my $i (9..$#a){
    unless($pop{$sample{$i}}){next;}
    my @fields =split(/:/,$a[$i]);
    my $dp = $fields[2];
    my $genotype = $fields[0];
    if (($genotype eq './.') or ($genotype eq '.')){next;}
    if ($dp < $min_depth){next;}
    if ($genotype eq '0/1'){
      $het{$pop{$sample{$i}}}{$chr}{$pos}++;
      $sample_count{$pop{$sample{$i}}}{$chr}{$pos}++;
    }else{
      $sample_count{$pop{$sample{$i}}}{$chr}{$pos}++;
    }
  }
}
print "pop\tchr\tpos\thet_percent\tn_samples";
foreach my $pop (sort keys %sample_count){
  foreach my $chr (sort keys %{$sample_count{$pop}}){
    foreach my $pos (sort {$a<=>$b} keys %{$sample_count{$pop}{$chr}}){
      unless($het{$pop}{$chr}{$pos}){
        $het{$pop}{$chr}{$pos} = 0;
      }
      unless($sample_count{$pop}{$chr}{$pos}){
        next;
      }
      my $het_percent = $het{$pop}{$chr}{$pos} / $sample_count{$pop}{$chr}{$pos};
      print "\n$pop\t$chr\t$pos\t$het_percent\$sample_count{$pop}{$chr}{$pos}";
    }
  }
}
