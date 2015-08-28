#!/bin/perl
use warnings;
use strict;
use Math::CDF;
use Statistics::Basic qw(:all);
use lib '/home/owens/bin/pop_gen/'; #For GObox server


my %t; #Convert from IUPAC to normal
$t{"N"} = "NN";
$t{"A"} = "AA";
$t{"T"} = "TT";
$t{"G"} = "GG";
$t{"C"} = "CC";
$t{"W"} = "AT";
$t{"R"} = "AG";
$t{"M"} = "AC";
$t{"S"} = "CG";
$t{"K"} = "GT";
$t{"Y"} = "CT";

#INPUT
#Hapmap piped in stdin
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
my $window_size = 0.5; #in cM
my $current_chrom;
my $end = $window_size;
my $min_sites = 1; #minimum number of sites used in a window
#Looping stats
my %likelihoodcount;
my $sitecount;
my %likelihood;

#Variables guessed from file, set for hapmap without iupac
my $iupac_coding= "False";
my $badcolumns="12";

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
		print  "sample\tchrom\tstart\tend\tpercentP2\tlikelihood";
	}else{
		next if /^\s*$/;
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
		my $cM = $a[4];
		if ($cM eq "NA"){
			next;
		}
		unless ($current_chrom){
			$current_chrom = $chrom;
		}
		if ($. == 2){
			until($end > $cM){ #move the end point of the window until it is after the current marker.
            		$end += $window_size;
            		}
		}
		if (($current_chrom ne $chrom) or ($end < $cM)){ #If it's starting a new window, do all the calculations
			my $start = $end - $window_size;
			#print FINALOUT "\n$current_chrom\t$start\t$end\t$sitecount";
			foreach my $i($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($likelihoodcount{$i}){
							if ($likelihoodcount{$i} > $min_sites){
								foreach my $percentP2 (0..100){
									$percentP2 = $percentP2/100;
									print "\n$samplename{$i}\t$current_chrom\t$start\t$end\t$percentP2\t$likelihood{$i}{$percentP2}";
											
								}

							}else{	
								#print FINALOUT "\tNA:NA:NA:NA";
							}
						}else{
							# print FINALOUT "\tNA:NA:NA:NA";
						}
					}
				}
			}
			#Reset the variables
			undef %likelihood;
			undef %likelihoodcount;
			$sitecount = 0;
			if ($current_chrom ne $chrom){ #If on the next chromosome, reset the end point for the window
				$current_chrom = $chrom;
				$end = $window_size;
			}
			until($end > $cM){ #move the end point of the window until it is after the current marker.
				$end += $window_size;
			}
			
		}
		my %BC;
		my %total_alleles;
		foreach my $i($badcolumns..$#a){
			if ($samplepop{$i}){
				$BC{"total"}{"total"}++;
				if ($iupac_coding eq "TRUE"){
					$a[$i] = $t{$a[$i]};
				}
				unless (($a[$i] eq "NN")or($a[$i] eq "XX")){
					my @bases = split(//, $a[$i]);
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
		if (keys %total_alleles == 2){
			$sitecount++;
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
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($BC{$i}{"Calls"}){
							$likelihoodcount{$i}++;
							foreach my $percentP2 (0..100){
                                                        	$percentP2 = $percentP2/100;
								my $percentP1 = 1 - $percentP2;
								if ($BC{$i}{$b1}){
									if ($BC{$i}{$b1} == 2){
										$likelihood{$i}{$percentP2} += log((($p1*$percentP1)+($p2*$percentP2)) * (($p1*$percentP1)+($p2*$percentP2)));
									}elsif ($BC{$i}{$b1} == 1){
										$likelihood{$i}{$percentP2} += log(2 * (($p1*$percentP1)+($p2*$percentP2)) * (($q1*$percentP1)+($q2*$percentP2)));
									}
								}else{
									$likelihood{$i}{$percentP2} += log((($q1*$percentP1)+($q2*$percentP2)) * (($q1*$percentP1)+($q2*$percentP2)));
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
			
			
