#!/bin/perl
use strict;
use warnings;
#This script correspondence measures for each pair of SNPs within a distance.

my $min_dif = 0.5; #minimum allele frequency difference between the parents
my $max_dist = 100000; # Max distance between markers to measure correspondence
my %sample;
my %data;
my %cmhash;
my %chrhash;
my %bphash;
my $counter1 = 0;
my @samplelist;
#Load in all data;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (5..$#a){
      $sample{$i} = $a[$i];
      push(@samplelist,$a[$i]);
    }
    next;
  }
  my $chr = $a[0];
  my $bp = $a[1];
  my $cm = $a[2];
  my $p1 = $a[3];
  my $p2 = $a[4];
  my $dif = abs($p1 - $p2);
  if ($dif < $min_dif){next;}
  $chrhash{$counter1} = $chr;
  $bphash{$counter1} = $bp;
  $cmhash{$chr}{$bp} = $cm;
  $counter1++;
  #Pick major allele for each parent
  my $major;
  my $minor;
  if ($p1 > $p2){
    $major = "p1";
    $minor = "p2";
  }else{
    $major = "p2";
    $minor = "p1";
  }
  foreach my $i (5..$#a){
    if ($a[$i] ne "N"){
      if ($a[$i] eq "2"){
        $data{$sample{$i}}{$chr}{$bp} = $major;
      }elsif ($a[$i] eq "0"){
        $data{$sample{$i}}{$chr}{$bp} = $minor;
      }
    }
  }
}
print "sample\tchr\tbp\tdist\tmatch";
foreach my $sample (@samplelist){

  foreach my $counter_current (0..($counter1-1)){
    foreach my $counter_match (0..($counter1-1)){
      if ($counter_current eq $counter_match){next;}      
      if ($chrhash{$counter_current} eq $chrhash{$counter_match}){
        my $chr1 = $chrhash{$counter_current};
        my $chr2 = $chrhash{$counter_match};
        my $bp1 = $bphash{$counter_current};
        my $bp2 = $bphash{$counter_match};
        my $dist = abs($bp1 - $bp2);
        if ($dist < $max_dist){
          if ($data{$sample}{$chr1}{$bp1} and $data{$sample}{$chr2}{$bp2}){
            if ($data{$sample}{$chr1}{$bp1} eq $data{$sample}{$chr2}{$bp2}){
              print "\n$sample\t$chr1\t$bp1\t$dist\t1";
            }else{
              print "\n$sample\t$chr1\t$bp1\t$dist\t0";
            }
          }
        }
      }  
    }
  }
  
}
