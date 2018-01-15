#!/bin/perl
use strict;
use warnings;
my $max_dist = 1;
my %sample;
my $current_chr = 0;
my %previous;
my %previous_data;
print "sample\tdist\tjunction";
while(<STDIN>){
  chomp;
  my @a = split(/\t/);
  if ($. == 1){
    foreach my $i (4..$#a){
      $sample{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $pos = $a[1];
    my $cm = $a[2];
    if ($cm eq "NA"){next};
    if ($chr ne $current_chr){
      undef(%previous_data);
      print STDERR "Now processing $chr...\n";
      $current_chr = $chr;
    }
    foreach my $i (4..$#a){
      my $sample = $sample{$i};
      my $current_state = $a[$i];
      if ($current_state eq "N"){next;}
      if ($current_state eq "1"){next;}
      my $current_cm = $cm;
      if ($previous_data{$sample}){
	foreach my $cm (keys %{$previous_data{$sample}}){
          my $current_dist = $current_cm - $cm;
          if ($current_dist > $max_dist){
	    delete($previous_data{$sample}{$cm});
	    next;
          }
	  my $previous_state = $previous_data{$sample}{$cm};
          if ($previous_state eq $current_state){
            print "\n$sample\t$current_dist\t0";
          }else{
            print "\n$sample\t$current_dist\t1";
          }
          $previous_data{$sample}{$current_cm} = $current_state;
        }
      }else{
        $previous_data{$sample}{$current_cm} = $current_state;
      }
    }
  }
}
