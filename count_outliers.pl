#!/bin/perl
use strict;
use warnings;

#This script takes a piped in fst file, it picks out the sites that are above an Fst threshold and records what alleles are for each group. It then uses a vcf file to record how many of each type of allele each sample has. 
#Assumes that the vcf is biallelic. 

my $fst_threshold = 0.95;

my %allele_id;
while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  if ($a[9] == "Inf"){next;}
  if ($a[9] > $fst_threshold){
    if ($a[13] > 0.5){ #The 0 allele is reference
      $allele_id{"$a[0]_$a[1]"}{0} = $a[2];
      $allele_id{"$a[0]_$a[1]"}{1} = $a[3];
    }else{
      $allele_id{"$a[0]_$a[1]"}{0} = $a[3];
      $allele_id{"$a[0]_$a[1]"}{1} = $a[2];
    } 
#    print STDERR "Outlier found at $a[0]_$a[1]\n"; 
  }
}

my $vcf = $ARGV[0];
open(IN, "gunzip -c $vcf |");

my %sample;
my %counts;
my $site_count = 0;
while(<IN>){
  chomp;
  if ($_ =~ m/^##/){
    next;
  }
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    print STDERR "Loading vcf file now...\n";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    next;
  }
  my $chr = $a[0];
  my $pos = $a[1];
  $site_count++;
  if ($site_count % 100000 == 0){print STDERR "$chr $pos processed...\n";}
  unless ($allele_id{"${chr}_$pos"}{1}){next;}
  my $ref = $a[3];
  my $alt = $a[4];
#  print "Found outlier at $chr $pos\n";
  if ($alt =~ m/,/){print "Multiallelic sites at $chr $pos.";exit; }
  foreach my $i (9..$#a){
    my @fields = split(/:/,$a[$i]);
    my $genotype = $fields[0];
    my @calls = split(/\//,$genotype);
    foreach my $j(0..1){
      $calls[$j] =~ s/0/$ref/g;
      $calls[$j] =~ s/1/$alt/g;
    }
    my $allele_0_count;
    my $allele_1_count;
    foreach my $j(0..1){
      if ($calls[$j] eq $allele_id{"${chr}_$pos"}{0}){
        $allele_0_count++;
      }elsif( $calls[$j] eq $allele_id{"${chr}_$pos"}{1}){
        $allele_1_count++;
      }

      if ($allele_0_count){
	if ($allele_0_count == 2){
          $counts{$sample{$i}}{'00'}++;
	  next;
        }
      }
      if ($allele_1_count){
        if ($allele_1_count == 2){
          $counts{$sample{$i}}{'11'}++;
          next;
        }
      }
      if (($allele_0_count) and ($allele_1_count)){
        if (($allele_0_count == 1) and ($allele_1_count == 1)){
          $counts{$sample{$i}}{'01'}++;
        }
      }
    }
  }  
}
print "name\t00\t01\t11";
foreach my $sample (sort keys %counts){
  my @genotypes = qw( 00 01 11 );
  foreach my $g (@genotypes){
    unless ($counts{$sample}{$g}){
      $counts{$sample}{$g} = 0;
    }
  }
  print "\n$sample\t$counts{$sample}{'00'}\t$counts{$sample}{'01'}\t$counts{$sample}{'11'}";
}
