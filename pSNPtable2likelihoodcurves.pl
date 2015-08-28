#!/bin/perl
use warnings;
use strict;
use Math::CDF;
use Statistics::Basic qw(:all);
use lib '/home/owens/bin/pop_gen/'; #For GObox server



#INPUT
#ptab piped in stdin
my $popfile = $ARGV[0]; #A population file.

my %pop;
my %popList;
my %samplepop;
my %samplename;


open POP, $popfile;
while (<POP>){ #Load in population information linked with sample name
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$popList{$a[1]}++;
}
close POP;
#Window variables
my $window_size = 1000000; #Number of bases in each sliding window
my $current_chrom;
my $end = $window_size;
my $min_sites = 1; #minimum number of sites used in a window
#Looping stats
my %likelihoodcount;
my $sitecount;
my %likelihood;

my $badcolumns="2";

while (<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	
	if ($. == 1){ #Load in sample names associated with column numbers, as well as population
		foreach my $i($badcolumns..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				$samplename{$i} = $a[$i];
			}
		}
		print "sample\tchrom\tstart\tend\tsites\tpercentP2\tlikelihood";
	}else{
		next if /^\s*$/;
		my $loc = "$a[0]\t$a[1]";
		my $chrom = $a[0];
		my $pos = $a[1];
		unless ($current_chrom){
			$current_chrom = $chrom;
		}
		if ($. == 2){
			until($end > $pos){ #move the end point of the window until it is after the current marker.
            		$end += $window_size;
            		}
		}
		if (($current_chrom ne $chrom) or ($end < $pos)){ #If it's starting a new window, do all the calculations
			my $start = $end - $window_size;
			foreach my $i($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($likelihoodcount{$i}){
							if ($likelihoodcount{$i} > $min_sites){
								foreach my $n (1..2){
									foreach my $percentP2 (0..100){
										$percentP2 = $percentP2/100;
										print "\n$samplename{$i}.$n\t$current_chrom\t$start\t$end\t$likelihoodcount{$i}\t$percentP2\t$likelihood{$i.$n}{$percentP2}";
									}
								}
							}
						}
					}
				}
			}
			#Reset the variables
			undef %likelihood;
			undef %likelihoodcount;
			if ($current_chrom ne $chrom){ #If on the next chromosome, reset the end point for the window
				$current_chrom = $chrom;
				$end = $window_size;
			}
			until($end > $pos){ #move the end point of the window until it is after the current marker.
				$end += $window_size;
			}
			
		}
		my %BC;
		my %hybrids;
		my %total_alleles;
		foreach my $i($badcolumns..$#a){
			if ($samplepop{$i}){
				$BC{"total"}{"total"}++;
				unless ($a[$i] eq "N|N"){
					my @bases = split(/\|/, $a[$i]);
					$total_alleles{$bases[0]}++;
					$total_alleles{$bases[1]}++;
				
					$BC{"total"}{$bases[0]}++;
		        		$BC{"total"}{$bases[1]}++;
					$BC{$samplepop{$i}}{$bases[0]}++;
		 			$BC{$samplepop{$i}}{$bases[1]}++;

					$BC{$i}{$bases[0]}++;
					$BC{$i}{$bases[1]}++;
					$BC{$i}{"Calls"}++;

					$BC{"total"}{"Calls"}++;
					$BC{$samplepop{$i}}{"Calls"}++;
				}
			}
		}
		foreach my $i($badcolumns..$#a){
			if ($samplepop{$i} eq "H"){
				unless ($a[$i] eq "N|N"){
					my @bases = split(/\|/, $a[$i]);
					$hybrids{$i.1}{$bases[0]}++;
					$hybrids{$i.2}{$bases[1]}++;
					$hybrids{$i.1}{"Calls"}++;
					$hybrids{$i.2}{"Calls"}++;
				}
			}
		}
		if (keys %total_alleles == 2){
			my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
			#Major allele
			my $b1 = $bases[1];
			#Minor allele
			my $b2 = $bases[0];
			unless (($BC{"P1"}{"Calls"}) and ($BC{"P2"}{"Calls"})){
				goto SKIP;
			}unless (($BC{"P1"}{"Calls"} >= 5) and ($BC{"P2"}{"Calls"} >= 5)){
				goto SKIP;
			}
			my $p1;
			my $p2;
			my $q1;
			my $q2;
			#Allele frequency of each allele in each population
			if ($BC{"P1"}{$b1}){
				$p1 = $BC{"P1"}{$b1}/($BC{"P1"}{"Calls"}*2);
			}else{
				$p1 = 0.01;
			}
			if ($BC{"P2"}{$b1}){
				$p2 = $BC{"P2"}{$b1}/($BC{"P2"}{"Calls"}*2);
			}else{
				$p2 = 0.01;
			}
			if ($BC{"P1"}{$b2}){
				$q1 = $BC{"P1"}{$b2}/($BC{"P1"}{"Calls"}*2);
			}else{
				$q1 = 0.01;
			}
			if ($BC{"P2"}{$b2}){
				$q2 = $BC{"P2"}{$b2}/($BC{"P2"}{"Calls"}*2);
			}else{
				$q2 = 0.01;
			}
			my %likelihoodP1;
			my %likelihoodP2;
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($hybrids{$i.1}{"Calls"}){
							$likelihoodcount{$i}++;
							foreach my $n (1..2){
								foreach my $percentP2(0..100){
									$percentP2 = $percentP2/100;
									my $percentP1 = 1 - $percentP2;
									if ($hybrids{$i.$n}{$b1}){
										$likelihood{$i.$n}{$percentP2} += log(($p1 * $percentP1) + ($p2 * $percentP2));
									}elsif($hybrids{$i.$n}{$b2}){
										$likelihood{$i.$n}{$percentP2} += log(($q1 * $percentP1) + ($q2 * $percentP2));
									}
	
								}
							}
						
						}
					}
				}
			}
		}
	}
	SKIP:
}
			
			
