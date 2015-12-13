#!/bin/perl
use warnings;
use strict;
use Math::CDF;
use Statistics::Basic qw(:all);
use lib '/home/owens/bin/pop_gen/'; #For GObox server
#This version uses a beta distribution to estimate allele frequency, then integrates that probability to estimate the likelihood for each parent

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
my $window_size = 1000000; #Number of bases in each sliding window
my $current_chrom;
my $end = $window_size;
my $min_sites = 10; #minimum number of sites used in a window
#Looping stats
my %likelihoodcount;
my $sitecount;
my %likelihood;

#Variables guessed from file, set for hapmap without iupac
my $iupac_coding= "False";
my $badcolumns="11";
my %good_number_hash;
my $counter;
while (<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	
	if ($. == 1){ #Load in sample names associated with column numbers, as well as population
		foreach my $i($badcolumns..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				$samplename{$i} = $a[$i];
				$good_number_hash{$i}++;
			}
		}
		print  "sample\tchrom\tstart\tend\tpercentP2\tlikelihood";
	}else{
		next if /^\s*$/;
		$counter++;
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
		if (($counter % 100000)== 0){
        	        print STDERR "Hyblik Processing $chrom $pos...\n";
        	}
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
			#print FINALOUT "\n$current_chrom\t$start\t$end\t$sitecount";
			foreach my $i($badcolumns..$#a){
				if ($samplepop{$i}){
#					if ($samplepop{$i} eq "H"){
					if (($samplepop{$i} eq "H") or ($samplepop{$i} eq "P1") or ($samplepop{$i} eq "P2")){
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
			until($end > $pos){ #move the end point of the window until it is after the current marker.
				$end += $window_size;
			}
			
		}
		my %BC;
		my %total_alleles;
		my $P1count = 0;
		my $P2count = 0;
		foreach my $i (keys %good_number_hash){ #Load up parental alleles
                	if ($a[$i] ne "NN"){
                        	if ($samplepop{$i} eq "P1"){
                        	        $P1count++;
                        	}elsif($samplepop{$i} eq "P2"){
                                	$P2count++;
                       		}
                	}
        	}
		unless(($P1count >=5) and ($P2count >= 5)){
                	next;
        	}
		my $min_count;
		if ($P1count < $P2count){
                	$min_count = $P1count;
        	}else{
                	$min_count = $P2count;
        	}
		my %P1alleles;
		my %P2alleles;
		my $P1count2 = 0;
		my $P2count2 = 0;
		foreach my $i(keys %good_number_hash){
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
					if (($samplepop{$i} eq "P1") and ($P1count2 < $min_count)){
						$P1alleles{$bases[0]}++;
						$P1alleles{$bases[1]}++;
						$P1alleles{"Calls"}++;
						$P1count2++;
					}elsif(($samplepop{$i} eq "P2") and ($P2count2 < $min_count)){
						$P2alleles{$bases[0]}++;
                                                $P2alleles{$bases[1]}++;
						$P2alleles{"Calls"}++;
                                                $P2count2++;
					}
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
			my $p1;
			my $p2;
			my $q1;
			my $q2;
			#Allele frequency of each allele in each population
			if ($P1alleles{$b1}){
				$p1 = $P1alleles{$b1}/($P1alleles{"Calls"}*2);
			}else{
				$p1 = 0.01;
			}
			if ($P2alleles{$b1}){
				$p2 = $P2alleles{$b1}/($P2alleles{"Calls"}*2);
			}else{
				$p2 = 0.01;
			}
			if ($P1alleles{$b2}){
				$q1 = $P1alleles{$b2}/($P1alleles{"Calls"}*2);
			}else{
				$q1 = 0.01;
			}
			if ($P2alleles{$b2}){
				$q2 = $P2alleles{$b2}/($P2alleles{"Calls"}*2);
			}else{
				$q2 = 0.01;
			}
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
#					if ($samplepop{$i} eq "H"){
					if (($samplepop{$i} eq "H") or ($samplepop{$i} eq "P1") or ($samplepop{$i} eq "P2")){
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
			
sub integrate {
	
