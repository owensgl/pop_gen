#!/bin/perl
use warnings;
use strict;

my $current_chr;
my $window_size = $ARGV[0];
my $start = 0;
my $end = $start+$window_size;

my $sitecounter = 0;
my %h;
my %counts;
my $verbose = "TRUE"; #Change to TRUE to turn on printing all counts
while(<STDIN>){
	   chomp;
       my @a = split(/\t/,$_);
       if ($. == 1){
           print "chr\tstart\tend\tmid\tsites";
           foreach my $i (2..$#a){
               print "\t$a[$i]";
           }
           next;
	      }else{
              my $chr = $a[0];
              my $pos = $a[1];
              unless($current_chr){
                  $current_chr = $chr;
              }
              if ($. == 2){
                  until ($pos < $end){
                      $start += $window_size;
                      $end = $start+$window_size;
                  }
              }
              if (($current_chr ne $chr) or ($pos >= $end)){
                  my $current_start = $start;
                  my $current_end = $end -1;
                  my $current_midpoint = $end - ($window_size/2);

                  print "\n$current_chr\t$current_start\t$current_end\t$current_midpoint";
                  print "\t$sitecounter";
                  foreach my $i (2..$#a){
                      print "\t";
                      foreach my $n(0..1){
			 my $sample = "$i.$n";
                         my @parents_tmp = (sort { $counts{$sample}{$b} <=> $counts{$sample}{$a}} keys %{$counts{$sample}});
                         if ($verbose ne "TRUE"){
                             #Check to make sure that the best choice is better than the next best
                             if ($counts{$sample}{$parents_tmp[0]} >$counts{$sample}{$parents_tmp[1]}){
                                 print "$parents_tmp[0]";
                             }else{
                                 print "NA";
                             }
                             if ($n == 0){
                                 print "|";
                             }
                         }else{ #While in verbose mode
                             my $firstprint;
                             foreach my $parent (@parents_tmp){
                                 unless ($firstprint){
                                     print "$parent=$counts{$sample}{$parent}";
                                     $firstprint++;
                                 }else{
                                     print ",$parent=$counts{$sample}{$parent}";
                                 }
                             }if ($n == 0){
				print "|";
			    }
                         }
                     }
                  }

			#Reset variables for new window
            %counts = ();
			$sitecounter = 0;
			if ($current_chr ne $chr){
				$current_chr = $chr;
				$start = 0;
				$end =  $start+$window_size;
			}
			until ($pos <$end){
				$start += $window_size;
				$end = $start+$window_size;
			}

		}
        foreach my $i (2..$#a){
            if ($a[$i] eq "NA"){
                next;
            }
            my @strands = split(/\|/,$a[$i]);
            foreach my $n (0..1){
                my @parents = split(/,/,$strands[$n]);
                foreach my $parent (@parents){
                    $counts{"$i.$n"}{$parent}++;
                }
            }
        }
		$sitecounter++;
	}
}
