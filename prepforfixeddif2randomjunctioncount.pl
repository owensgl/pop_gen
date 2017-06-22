#!/bin/perl
use warnings;
use strict;
use List::Util qw(shuffle);
#This script takes the output from "/home/owens/bin/reformat/ctab2prepforparentblocks.pl" and calculates the number of junctions under random assortment of  
my $current_chr;
my %current_state;
my %name;
my %junctions;
my $total_columns;
my $counter;
my %data;
my %chrom_hash;
my %ancestry_mix;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (4..$#a){
      $name{$i} = $a[$i];
      $total_columns = $#a;
    }
  }else{
    my $chr = $a[0];
    if ($chr =~ m/Chr00/){next;}
    $counter++;
    $chrom_hash{$counter} = $chr;
    foreach my $i (4..$#a){
      $data{$i}{$counter} = $a[$i];
    }
  }
}
print "sample\trep\tjunctions";
foreach my $i (4..$total_columns){
  foreach my $rep (1..1000){
    my $junc_counter;
    my $current_state;
    my $current_chr;
    my $n;
    foreach my $key (shuffle( keys %{$data{$i}}  )){
      $n++;
      unless($current_chr){
	$current_chr = $chrom_hash{$n};
      }
      if ($current_chr ne $chrom_hash{$n}){
	undef($current_state);
	$current_chr = $chrom_hash{$n};
      }
      if ($data{$i}{$key} ne "N"){
	if(defined $current_state){ #If it has a previous state.
          my $state_dif = abs($current_state - $data{$i}{$key})/2;
          $junc_counter+=$state_dif;
          $current_state = $data{$i}{$key};
        }else{
          $current_state = $data{$i}{$key};
        }
      }
    }
    print "\n$name{$i}\t$rep\t$junc_counter";
  }
}
