#!/bin/perl
use warnings;
use strict;

my $current_chr;
my $window_size = $ARGV[0];
my $start = 0;
my $end = $start+$window_size;


my %h;
my $sitecounter = 0;
my $variablecounter = 0;

#starting variables
$h{"fstnum"}=0;
$h{"fstdenom"}=0;
$h{"Hexp1"}=0;
$h{"Hexp2"}=0;
while(<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	if ($. == 1){
		print "Chr\tStartPos\tEndPos\tMidPos\tTotalSites\tVariableSites\tFst\tHexp1\tHexp2";
		next;
	}else{
		my $chr = $a[1];
		my $pos = $a[2];
		my $fstnum = $a[6];
		if ($fstnum eq "NA"){
			next;
		}
		my $fstdenom = $a[7];
		my $Hexp1 = $a[9];
		my $Hexp2 = $a[10];
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
			my $current_fst;
			if ($h{"fstdenom"} eq 0){
				$current_fst = "NA";
			}else{
				$current_fst = $h{"fstnum"}/$h{"fstdenom"};
			}
			my $current_Hexp1 = $h{"Hexp1"}/$sitecounter;
			my $current_Hexp2 = $h{"Hexp2"}/$sitecounter;
			print "\n$current_chr\t$current_start\t$current_end\t$current_midpoint";
			print "\t$sitecounter\t$variablecounter\t$current_fst\t$current_Hexp1\t$current_Hexp2";

			#Reset variables for new window
			$h{"fstnum"}=0;
			$h{"fstdenom"}=0;
			$h{"Hexp1"}=0;
			$h{"Hexp2"}=0;
			$sitecounter = 0;
			$variablecounter = 0;
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
		$h{"fstnum"}+=$fstnum;
		$h{"fstdenom"}+=$fstdenom;
		$h{"Hexp1"}+=$Hexp1;
		$h{"Hexp2"}+=$Hexp2;
		$sitecounter++;
		if ($fstdenom ne 0){
			$variablecounter++;
		}
	}
}
