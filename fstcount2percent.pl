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
  print "$name\t$percent_1\n";
}
