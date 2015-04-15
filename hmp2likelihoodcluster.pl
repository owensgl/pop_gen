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
my $out = $ARGV[1]; #outfile prefix

open (FINALOUT, "> $out.final.txt") or die "Could not open a file\n"; #open outfile
open (LLOUT, "> $out.likelihoods.txt") or die "Could not open a file\n"; #open population file
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
my %likelihoodP1;
my %likelihoodP2;
my %likelihoodDIF;
my %totalP1;
my %totalP2;
my %likelihoodcount;
my $sitecount;
my %totaLP1;
my %totaLP2;

#Variables guessed from file, set for hapmap without iupac
my $iupac_coding= "False";
my $badcolumns="11";

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
		print LLOUT "chom\tsite";
		print FINALOUT "chrom\tstart\tend\tn_sites";
		foreach my $i($badcolumns..$#a){
			if ($samplepop{$i}){
				if ($pop{$a[$i]} eq "H"){
					print LLOUT"\t$samplename{$i}";
					print FINALOUT "\t$samplename{$i}";
				}
			}
		}
	}else{
		next if /^\s*$/;
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
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
			print FINALOUT "\n$chrom\t$start\t$end\t$sitecount";
			foreach my $i($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($likelihoodcount{$i}){
							if ($likelihoodcount{$i} > $min_sites){
								my $LR = $totaLP1{$i} - $totaLP2{$i};
								my @likelihoodarray = @{$likelihoodDIF{$i}}; 
								my $std = stddev(@likelihoodarray);
								my $Nsquare = sqrt($likelihoodcount{$i});
								my $Zscore = $LR / ($Nsquare * $std);
								my $pvalue = Math::CDF::pnorm(abs($Zscore));
								my $outcome;
								if (($Zscore > 0) and ($pvalue > 0.95)){
									$outcome = "P1";
								}elsif (($Zscore < 0) and ($pvalue > 0.95)){
									$outcome = "P2";
								}else{
									$outcome = "NA";
								}
								print FINALOUT "\t$outcome:$Zscore:$pvalue";
							}else{	
								print FINALOUT "\tNA:NA:NA";
							}
						}else{
							 print FINALOUT "\tNA:NA:NA";
						}
					}
				}
			}
			#Reset the variables
			undef %totalP1;
			undef %totalP2;
			undef %likelihoodP1;
			undef %likelihoodP2;
			undef %likelihoodDIF;
			undef %likelihoodcount;
			undef $sitecount;
			undef %totaLP1;
			undef %totaLP2;
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
			print LLOUT "\n$chrom\t$pos";
			my %likelihoodP1;
			my %likelihoodP2;
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($BC{$i}{"Calls"}){
							$likelihoodcount{$i}++;
							if ($BC{$i}{$b1}){
								if ($BC{$i}{$b1} == 2){
									$likelihoodP1{$i} = log($p1 * $p1);
									$likelihoodP2{$i} = log($p2 * $p2);
									push @{ $likelihoodDIF{$i} }, (log($p1 * $p1) - log($p2 * $p2));
									$totaLP1{$i} += log($p1 * $p1);
									$totaLP2{$i} += log($p2 * $p2);
								}elsif ($BC{$i}{$b1} == 1){
									$likelihoodP1{$i} = log(2 * ($p1 * $q1));
									$likelihoodP2{$i} = log(2 * ($p2 * $q2));
									push @{ $likelihoodDIF{$i} }, (log(2 * ($p1 * $q1)) - log(2 * ($p2 * $q2)));
									$totaLP1{$i} += log(2 * ($p1 * $q1));
									$totaLP2{$i} += log(2 * ($p2 * $q2));
								}
							}else{
								$likelihoodP1{$i} =  log($q1 * $q1);
								$likelihoodP2{$i} = log($q2 * $q2);
								push @{ $likelihoodDIF{$i} }, (log($q1 * $q1) - log($q2 * $q2));
								$totaLP1{$i} += log($q1 * $q1);
								$totaLP2{$i} += log($q2 * $q2);		
							}
						}
						else{
							$likelihoodP1{$i} = "NA";
							$likelihoodP2{$i} = "NA";
						}
					}
				}
			}
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						print LLOUT "\t$likelihoodP1{$i}:$likelihoodP2{$i}";
					}
				}
			}
		}
	}
	SKIP:
}
			
			
