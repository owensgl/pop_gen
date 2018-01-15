#!/bin/perl
use warnings;
use strict;
use POSIX;
#This will summarize emmax pvalues by gene.

#my $geneloc = $ARGV[0];
my $geneloc = "/home/owens/ref/HanXRQr1.0-20151230-EGN-r1.1.genelocations.txt";
#my $geneloc = "/home/owens/ref/HanXRQr1.0-20151230_genes_ubc_v1.0.genes.bed";
my %genes_by_chr;
my %genes_by_block;
my %gene_start;
my %gene_end;
my %chrlist;
open GENE, $geneloc;
while(<GENE>){
  chomp;
  my @a = split(/\t/,$_);
  my $gene = $a[0];
  my $chr = $a[1];
  $chr =~ s/HanXRQChr//g;
  my $start = $a[2];
  my $stop = $a[3];
  my $mb_block = floor($start / 1000000);
  push(@{$genes_by_block{$chr}{$mb_block}},$gene);
  push(@{$genes_by_chr{$chr}},$gene);
  $gene_start{$gene}=$start;
  $gene_end{$gene} = $stop;
  $chrlist{$chr}++;
}
#pipe in the output of emmax
my %pvalues;
my %sitecount;
my $sig;
my $prog; #ANGSD or EMMAX
while(<STDIN>){
  chomp;
  if ($. == 1){
    my @a = split(/\t/,$_);
    if ($#a == 3){
      $sig = "pvalue";
      $prog = "emmax";
    }elsif ($#a == 7){
      $sig = "lrt";
      $prog = "angsd";
    }else{
      die "Unrecognized format based on column number\n";
    }
    next;
  }
  my @a = split(/\t/,$_);
  my $pos;
  my $chr;
  my $pvalue;
  if ($prog eq "emmax"){
    my $locus = $a[0];
    my @info = split('_',$locus);
    $chr = $info[0];
    $pos = $info[1];
    $pvalue = $a[3];
  }elsif ($prog eq "angsd"){
    $chr = $a[0];
    $pos = $a[1];
    $pvalue = $a[6];
    if ($pvalue eq "-999.000000"){
      next; #Skip snps without enough data;
    }
    $chr =~ s/HanXRQChr//g;
  }
  my $mb_block = floor($pos / 1000000);
  foreach my $gene (@{$genes_by_block{$chr}{$mb_block}}){
    if (($gene_start{$gene} <= $pos) and ($gene_end{$gene} >= $pos)){
      $pvalues{$gene}{$pvalue}++;
      $sitecount{$gene}++;
      goto NEXTLINE;
    }elsif ($gene_start{$gene} > $pos){
      goto TRYPREVIOUSBLOCK;
    }
  }
  TRYPREVIOUSBLOCK:  
  if ($mb_block > 0){
    my $previous_block = $mb_block - 1;
    foreach my $gene (@{$genes_by_block{$chr}{$previous_block}}){
      if (($gene_start{$gene} <= $pos) and ($gene_end{$gene} >= $pos)){
        $pvalues{$gene}{$pvalue}++;
        $sitecount{$gene}++;
        goto NEXTLINE;
      }elsif ($gene_start{$gene} > $pos){
        goto NEXTLINE;
      }
    }
  }
  NEXTLINE:
}
if ($prog eq "emmax"){
  print "Chr\tStartPos\tEndPos\tGene\tTotalSites\tlowestp\tsecondp";
}elsif ($prog eq "angsd"){
  print "Chr\tStartPos\tEndPos\tGene\tTotalSites\thighest_lrt\tsecond_lrt";
}
foreach my $chr (sort keys %chrlist){
  foreach my $gene (@{$genes_by_chr{$chr}}){
    if ($sitecount{$gene}){
      my @pvalues;
      if ($prog eq "emmax"){    
        @pvalues = sort {$a<=>$b} keys %{$pvalues{$gene}};
      }elsif ($prog eq "angsd"){
        @pvalues = sort {$b<=>$a} keys %{$pvalues{$gene}};
      }
      my $lowest = $pvalues[0];
      my $secondlowest = "NA";
      if ($pvalues[1]){
        $secondlowest = $pvalues[1];
      }
      print "\n$chr\t$gene_start{$gene}\t$gene_end{$gene}\t";
      print "$gene\t$sitecount{$gene}\t$lowest\t$secondlowest";
    }
  }
}
