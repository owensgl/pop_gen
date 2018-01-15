#!/bin/perl
use warnings;
use strict;

#This will summarize Fst and heterozygosity by genes.

my $geneloc = $ARGV[0];
my %genes_by_chr;
my %gene_start;
my %gene_end;
my %chrlist;
open GENE, $geneloc;
while(<GENE>){
  chomp;
  my @a = split(/\t/,$_);
  my $gene = $a[0];
  my $chr = $a[1];
  my $start = $a[2];
  my $stop = $a[3];
  push(@{$genes_by_chr{$chr}},$gene);
  $gene_start{$gene}=$start;
  $gene_end{$gene} = $stop;
  $chrlist{$chr}++;
}
#pipe in the output of SNPtable2fst.pl
my %fst_num;
my %fst_denom;
my %hexp1;
my %hexp2;
my %sitecount;
while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[1];
  my $pos= $a[2];
  my $num = $a[6];
  my $denom = $a[7];
  my $hexp1 = $a[9];
  my $hexp2 = $a[10];
  foreach my $gene (@{$genes_by_chr{$chr}}){
    if (($gene_start{$gene} <= $pos) and ($gene_end{$gene} >= $pos)){
      $fst_num{$gene} +=$num;
      $fst_denom{$gene} += $denom;
      $hexp1{$gene} += $hexp1;
      $hexp2{$gene} += $hexp2;
      $sitecount{$gene}++;
    }
  }
}
print "Chr\tStartPos\tEndPos\tGene\tTotalSites\tFst\tHexp1\tHexp2";
foreach my $chr (sort keys %chrlist){
  foreach my $gene (@{$genes_by_chr{$chr}}){
    if ($sitecount{$gene}){
      my $fst = $fst_num{$gene} / $fst_denom{$gene};
      my $hexp1 = $hexp1{$gene}/$sitecount{$gene};
      my $hexp2 = $hexp2{$gene}/$sitecount{$gene};
      print "\n$chr\t$gene_start{$gene}\t$gene_end{$gene}\t";
      print "$gene\t$sitecount{$gene}\t$fst\t$hexp1\t$hexp2";
    }
  }
}
