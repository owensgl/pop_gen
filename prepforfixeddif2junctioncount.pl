#!/bin/perl
use warnings;
use strict;

#This script takes the output from "/home/owens/bin/reformat/ctab2prepforparentblocks.pl" and calculates the number of junctions. 
my $current_chr;
my %current_state;
my %name;
my %junctions;
my $total_columns;
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
    unless($current_chr){
      $current_chr = $chr;
    }
    if ($current_chr ne $chr){
      undef(%current_state);
      $current_chr = $chr;
    }
    foreach my $i (4..$#a){
      if ($a[$i] ne "N"){
        if(defined $current_state{$i}){ #If it has a previous state.
          my $state_dif = abs($current_state{$i} - $a[$i])/2;
          $junctions{$i}+=$state_dif;
          $current_state{$i} = $a[$i];
          
        }else{
          $current_state{$i} = $a[$i];
        }
      }
    }
  }
}
print "sample\tjunctions";
foreach my $i (4..$total_columns){
  print "\n$name{$i}\t$junctions{$i}";
}
