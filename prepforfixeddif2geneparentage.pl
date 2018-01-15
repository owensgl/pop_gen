#!/bin/perl
use warnings;
use strict;

#This will summarize genes to see if they have SNPs from both parents

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
#
my %sample;
my %data;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (4..$#a){
      $sample{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $pos = $a[1];
    foreach my $gene (@{$genes_by_chr{$chr}}){
      if (($gene_start{$gene} <= $pos) and ($gene_end{$gene} >= $pos)){
        foreach my $i (4..$#a){
          if ($a[$i] eq "N"){next;}
          if ($a[$i] eq "0"){
	    $data{$sample{$i}}{$gene}{2}+=2;
	  }elsif ($a[$i] eq "2"){
            $data{$sample{$i}}{$gene}{1}+=2;
          }elsif ($a[$i] eq "1"){
            $data{$sample{$i}}{$gene}{2}++;
            $data{$sample{$i}}{$gene}{1}++;
	  }
        }
      }
    }
  }
}

print "sample\tChr\tStartPos\tEndPos\tGene\tP1alleles\tP2alleles\ttotal_alleles\tmixed";
foreach my $chr (sort keys %chrlist){
  foreach my $gene (@{$genes_by_chr{$chr}}){
    foreach my $sample (sort keys %data){
      foreach my $n (1..2){
        unless($data{$sample}{$gene}{$n}){
          $data{$sample}{$gene}{$n} = 0;
        }
      }
      my $mixed = 0;
      if (( $data{$sample}{$gene}{1}) and ($data{$sample}{$gene}{2})){
        $mixed = 1;
      }
      my $total = $data{$sample}{$gene}{1} + $data{$sample}{$gene}{2};
      if ($total == 0){next;}
      
      print "\n$sample\t$chr\t$gene_start{$gene}\t$gene_end{$gene}\t";
      print "$gene\t$data{$sample}{$gene}{1}\t$data{$sample}{$gene}{2}\t$total\t$mixed";
    }
  }
}
