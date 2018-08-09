#!/bin/perl
use strict;
use warnings;
#This counts the percent of 1/(n-1) heterozygotes that may come from barcode switching at different depths. From 5 to 15 total depth at the site.
my %sample;
my $counter;
my %counts;
my %hets;
while(<STDIN>){
  my $line = "$_";
  chomp $line;
  my @fields = split /\t/,$line;
  if($line=~m/^##/){
   next;
  }
  if ($line =~m/^#CHROM/){
   foreach my $i (9..$#fields){
    $sample{$i} = $fields[$i];
   }
  }
  else{
   $counter++;
   if ($counter % 500000 == 0){print STDERR "Currently processed $counter loci\n";}
   foreach my $i (9..$#fields){
    if ($fields[$i] ne '.'){
     my @info = split(/:/,$fields[$i]);
     if ($info[0] eq './.'){next;}
     my @genotype = split(/\//,$info[0])
;
     if ($genotype[0] ne $genotype[1]){
       $hets{$sample{$i}}++;
     }
     $counts{$sample{$i}}++;
    }
   }
  }
}
print "name\ttotal_sites\thets\tpercent_het";
foreach my $sample (sort values %sample){
 my $percent_het = "NA";
 if ($counts{$sample} > 0){
   $percent_het = $hets{$sample}/$counts{$sample};
 }
 print "\n$sample\t$counts{$sample}\t$hets{$sample}\t$percent_het";
}
