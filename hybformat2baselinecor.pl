#!/bin/perl
use strict;
use warnings;
#This script creates baseline levels of correlation.

my $min_dif = 0.5; #minimum allele frequency difference between the parents
my $max_dist = 100000; # Max distance between markers to measure correspondence
my $n_sites = 1000; #Number of target sites to test
my $x_sites = 1000; #Number of sites to test correlation with for target site.
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
      push(@samplelist, $a[$i]);
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
print "sample\tpercent";
foreach my $sample (@samplelist){
  my $counter2 = 0;

  until ($counter2 == $n_sites){
    my $counter3 = 0;
    my $matches;
    until ($counter3 == $x_sites){
        
      my $n1 = int(rand($counter1));
      my $n2 = int(rand($counter1));
      if ($chrhash{$n1} ne $chrhash{$n2}){
        my $chr1 = $chrhash{$n1};
        my $chr2 = $chrhash{$n2};
        my $bp1 = $bphash{$n1};
        my $bp2 = $bphash{$n2};
	if ($data{$sample}{$chr1}{$bp1} and $data{$sample}{$chr2}{$bp2}){
          if ($data{$sample}{$chr1}{$bp1} eq $data{$sample}{$chr2}{$bp2}){
	    $matches++;
	  }
	}else{
	  next;
	}
      }else{
        next;
      }
      $counter3++;
    }
    my $percent = $matches/$x_sites;
    print "\n$sample\t$percent";
    $counter2++;
    
  }
}
