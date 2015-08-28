#!/bin/perl
use warnings;
use strict;

my $window_size = 75000;
my $current_chrom;
my $end = $window_size;
my $sitecount;

my @hexp1;
my @hexp2;
while (<STDIN>){
	chomp;
	if ($. == 1){
		print "chr\tstart\tend\tsites\tHexp1\tHexp2\tHs";
		next;
	}
	my @a = split(/\t/,$_);
	next if $a[9] eq "Inf";

	my $chrom = $a[1];
	my $pos = $a[2];	
	my $hexp1 = $a[10];
	my $hexp2 = $a[11];
	unless ($current_chrom){
		$current_chrom = $chrom;
	}
	if ($. == 2){
		until($end > $pos){
			$end += $window_size;
		}
	}
	if (($current_chrom ne $chrom) or ($end < $pos)){
		my $start = $end - $window_size;
		my $tmp_end = $end -1;
		unless ($sitecount ){
			print "\n$current_chrom\t$start\t$tmp_end\t0\tNA\tNA\tNA";
			goto SKIP;
		}
		my $hexp1_mean = average(@hexp1);
		my $hexp2_mean = average(@hexp2);
		push @hexp1, @hexp2;
		my $hexp_mean = average(@hexp1);
		print "\n$current_chrom\t$start\t$tmp_end\t$sitecount\t$hexp1_mean\t$hexp2_mean\t$hexp_mean";
		SKIP:
		undef @hexp1;
		undef @hexp2;
		$sitecount = 0;
		if ($current_chrom ne $chrom){
			$current_chrom = $chrom;
			$end = $window_size;
		}
		until ($end > $pos){
			$end += $window_size;
		}
	}
	push (@hexp1, $hexp1);
	push (@hexp2, $hexp2);
	$sitecount++;
}




sub average {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array 
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}
	
