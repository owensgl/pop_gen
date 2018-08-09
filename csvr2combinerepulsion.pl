#!/bin/perl
use warnings;
use strict;

#This script takes two csvr files for two haplotypes in repulsion. It combines the information for both haplotypes and keeps the phase of the first script. It does not impute anything.

my $csvr1 = $ARGV[0];
my $csvr2 = $ARGV[1];

my $output_csvr1 = $csvr1;
$output_csvr1 =~ s/.csvr/.combined.csvr/g;

my %reverse;
$reverse{"H"} = "A";
$reverse{"A"} = "H";
$reverse{"-"} = "-";

my %header;
my $headercounter = 0;
my %data;
my $nsamples;
my %loci_label;
#Load in data for the first haplotype, not reversed.
open CSVR1, $csvr1;
while(<CSVR1>){
  chomp;
  unless ($_ =~ m/^HanXRQ/){
    $header{$headercounter} = $_;
    $headercounter++;
  }else{
    my @a = split(/,/,$_);
    my $locus = $a[0];
    my $chr = $a[1];
    my $cm = $a[2];
    $nsamples = $#a;
    $loci_label{$chr}{$cm} = $locus;
    foreach my $i (3..$#a){
      if ($a[$i] eq "-"){next;}
      $data{$chr}{$cm}{$i} = $a[$i];
    }
  }
}
close CSVR1;

#Load in data for the second haplotype. This one reverse the genotypes
open CSVR2, $csvr2;
while(<CSVR2>){
  chomp;
  unless ($_ =~ m/^HanXRQ/){
  }else{
    my @a = split(/,/,$_);
    my $locus = $a[0];
    my $chr = $a[1];
    my $cm = $a[2];
    $loci_label{$chr}{$cm} = $locus;
    foreach my $i (3..$#a){
      if ($a[$i] eq "-"){next;}
      $data{$chr}{$cm}{$i} = $reverse{$a[$i]};
    }
  }
}
close CSVR2;
#Print out combined file 1.
open(my $out1, '>', $output_csvr1);
#print header
foreach my $n (0..($headercounter-2)){
  print $out1 "$header{$n}\n";
}
print $out1 "$header{($headercounter-1)}\n";

foreach my $chr (sort keys %data){
  my @loci = (sort keys %{$data{$chr}});
  foreach my $loci (@loci){
    print $out1 "\n$loci_label{$chr}{$loci},$chr,$loci";
    foreach my $sample (3..$nsamples){
      if ($data{$chr}{$loci}{$sample}){
        print $out1 ",$data{$chr}{$loci}{$sample}";
      }else{
        print $out1 ",-";
      }
    }
  }
}

