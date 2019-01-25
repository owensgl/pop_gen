#!/bin/perl
use warnings;
use strict;

while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $name = $a[0];
  my $count00 = $a[1];
  my $count01 = $a[2];
  my $count11 = $a[3];
  my $percent_1 = (($count00 * 2) + $count01) / (($count00 + $count01 + $count11)*2);
  my $call;
  if ($percent_1 <= 0.25){
    $call = "p0";
  }elsif ($percent_1 >= 0.75){
    $call = "p2";
  }else{
    $call = "p1";
  }
  print "$name\t$call\n";
}
