#!/bin/perl
use warnings;
use strict;

my $current_chr;
my $window_size = $ARGV[0];
my $speciesfile = $ARGV[1];
#Loads up names to species
my %specieshash;
open SPECIES, $speciesfile;
while(<SPECIES>){
	chomp;
	my @a = split(/\t/,$_);
	$specieshash{$a[0]} = $a[1];
}
close SPECIES;

my $start = 0;
my $end = $start+$window_size;


my %h;


#starting variables

while(<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	if ($. == 1){
		print "Sample\tSpecies\tChr\tStartPos\tEndPos\tMidPos\tTotalSites\tP1_sites\tP2_sites\tUniqueSites";
		next;
	}else{
		my $chr = $a[0];
		my $pos = $a[1];
		my $sample = $a[2];
		my $type = $a[4];
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
			my @samplelist = sort keys %h;
			foreach my $current_sample(@samplelist){
				my $P1count = 0;
				my $P2count = 0;
				my $uniquecount = 0;
				if ($h{$current_sample}{"0"}){
					$uniquecount = $h{$current_sample}{"0"};
				}
				if ($h{$current_sample}{"1"}){
					$P1count = $h{$current_sample}{"1"};
				}
				if ($h{$current_sample}{"2"}){
					$P2count = $h{$current_sample}{"2"};
				}
				my $total = $P1count + $P2count + $uniquecount;
				print "\n$current_sample\t$specieshash{$current_sample}\t$current_chr\t$current_start\t$current_end\t$current_midpoint\t$total\t$P1count\t$P2count\t$uniquecount";
			}

			#Reset variables for new window
			undef(%h);
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
		$h{$sample}{$type}++;
	}
}
