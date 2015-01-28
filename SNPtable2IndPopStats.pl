#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin/pop_gen/'; #For GObox server
my %t;
$t{"N"} = "NN";
$t{"A"} = "AA";
$t{"T"} = "TT";
$t{"G"} = "GG";
$t{"C"} = "CC";
$t{"W"} = "TA";
$t{"R"} = "AG";
$t{"M"} = "AC";
$t{"S"} = "CG";
$t{"K"} = "TG";
$t{"Y"} = "CT";

my %samples;
my @Good_samples;
my %Anc;
my %AncCount;
my %TotalSites;
my %pop;
my %Genotype;

my %samplepop;

my %poplist;
my %Pi;
my %Sites;
my %Polymorphic;
my %HetObs;
my %Invariant;


unless (@ARGV == 2) {die;}
my $in = $ARGV[0];
my $pop = $ARGV[1];

require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;

open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
}
close POP;


open IN, $in;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i ($badcolumns..$#a){ #Get sample names for each column
        		if ($pop{$a[$i]}){
        			$samplepop{$i} = $pop{$a[$i]};
        		}
        	}
	}
	else{
		next if /^\s*$/;
		my %BC;
		my %Counts;
		my %Hets;
		my %total_alleles;
		foreach my $i ($badcolumns..$#a){
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

					$BC{"total"}{"Calls"}++;
					$Counts{$samplepop{$i}}++;
					
					if($bases[0] ne $bases[1]){
						$BC{"total"}{"Het"}++;
						$Hets{$samplepop{$i}}++;
					}
				}
			} 
		}

		

		
		if (keys %total_alleles <= 2){
			foreach my $pop (sort(keys %poplist)){
				if ($Counts{$pop}){			
					$Sites{$pop}++;
				}
				#my $tmp = scalar keys %{$BC{$pop}};
				#my $tmp2 = scalar keys %total_alleles;
				if (scalar keys %{$BC{$pop}} > 1){
					$Polymorphic{$pop}++;
				}elsif (scalar keys %{$BC{$pop}} eq 1){
					$Invariant{$pop}++;
				}
				if ($Hets{$pop}){
					$HetObs{$pop} += ($Hets{$pop} / $Counts{$pop});
				}
				if ($Counts{$pop}){
					my @bases = sort { ${$BC{$pop}}{$a} <=> ${$BC{$pop}}{$b} } keys %{$BC{$pop}} ;
					if ($#bases == 1){
						my $b1 = $bases[1];
						my $b2 = $bases[0];
						
						my $p = $BC{$pop}{$b1}/($Counts{$pop}*2);
						my $q = $BC{$pop}{$b2}/($Counts{$pop}*2);
						$Pi{$pop} += ($p * $q * 2);
					}
				}
			
			}
		}
	}
}

print "Population\tSites\tPolymorphic\tInvariant\tHetObs\tPi\tFis";
foreach my $pop (sort(keys %poplist)){
	print "\n$pop\t$Sites{$pop}\t$Polymorphic{$pop}\t$Invariant{$pop}\t";
#	print "$HetObs{$pop}\t";
#	print "$Pi{$pop}\t";
	my $Hetobserved = ($HetObs{$pop} / $Sites{$pop});
	my $Piobserved = ($Pi{$pop} / $Sites{$pop});
	my $Fisobserved = (($Piobserved - $Hetobserved)/ $Piobserved);
	print "$Hetobserved\t$Piobserved\t$Fisobserved";
}
	
close IN;

