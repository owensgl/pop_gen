#!/usr/bin/perl

use warnings;
use strict;

my $POP = $ARGV[0];
my $HMP = $ARGV[1];




my %samples;
my @Good_samples;
my %Anc;
my %AncCount;
my %TotalSites;
my %pop;
my %Genotype;

my %samplepop;

my %poplist;
my $Npops = 2;
my $BiCount = 0;
my $TriCount= 0;
my $QuadCount = 0;
my $SingleTri =0;



my %pop2cols;
my %PPN2POP;

open IN, $POP;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
	$PPN2POP{$a[0]} = $a[1];
}
close IN;

open IN, $HMP;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i (2..$#a){ #Get sample names for each column
			my $pop = $PPN2POP{$a[$i]};
			push(@{$pop2cols{$pop}},$i);

        	}
        
		print "Chr\tpos\tpop1\tpop2\tContrast\tN1\tN2\tNTotal\tDxy\tFstNum\tFstDenom\tFst\tHexp1\tHexp2\tFreqDif\n";
	}
	else{
		my %done;
		foreach my $pop1 (sort keys %pop2cols){
			foreach my $pop2 (sort keys %pop2cols){
				if($pop1 ne $pop2){
					unless ($done{"$pop1-$pop2"}){
						$done{"$pop1-$pop2"}++;
						$done{"$pop2-$pop1"}++;
						print "$a[0]\t$a[1]\t$pop1\t$pop2\t$pop2-$pop1";

						my %BC;
						my %BS;
						my %total_alleles;
		
						#pop1
						foreach my $i (@{$pop2cols{$pop1}}){
							$BC{"total"}{"total"}++;
							my $tmp = $a[$i];
							unless (($tmp eq "NN")or($tmp eq "XX")){
								my @bases = split(//, $tmp);
								$total_alleles{$bases[0]}++;
								$total_alleles{$bases[1]}++;
								$BC{"total"}{$bases[0]}++;
							       	$BC{"total"}{$bases[1]}++;
								$BC{"1"}{$bases[0]}++;
							 	$BC{"1"}{$bases[1]}++;
								$BC{"total"}{"Calls"}++;
								$BC{"1"}{"Calls"}++;
								if($bases[0] ne $bases[1]){
									$BC{"total"}{"Het"}++;
									$BC{"1"}{"Het"}++;
								
								}
							}
						}
						#pop2
						foreach my $i (@{$pop2cols{$pop2}}){
							$BC{"total"}{"total"}++;
							my $tmp = $a[$i];
							unless (($tmp eq "NN")or($tmp eq "XX")){
								my @bases = split(//, $tmp);
								$total_alleles{$bases[0]}++;
								$total_alleles{$bases[1]}++;
								$BC{"total"}{$bases[0]}++;
							       	$BC{"total"}{$bases[1]}++;
								$BC{"2"}{$bases[0]}++;
							 	$BC{"2"}{$bases[1]}++;
								$BC{"total"}{"Calls"}++;
								$BC{"2"}{"Calls"}++;
								if($bases[0] ne $bases[1]){
									$BC{"total"}{"Het"}++;
									$BC{"2"}{"Het"}++;
				
								}
							}
						}
						my $pAll;
						my $qAll;
						my $HeAll;
						my $HoAll;
						my $CallRate;
						my $p1;
						my $q1;
						my $p2;
						my $q2;
						my $Ho1;
						my $Ho2;
						my $He1 ;
						my $He2;
						my $dxy;
						my $HsBar;
						my $H_bar;
						my $n_bar;
						my $n_1;
						my $n_2;
						my $n_total;
						my $sigma_squared;
						my $n_c;
						my $WC_a;
						my $WC_b;
						my $WC_c;
						my $WC_denom;
						my $WC_fst;
						my $pi;
						my $freq_dif;
		

						unless ($BC{"total"}{"Calls"}){
							$BC{"total"}{"Calls"} = 0;
						}
		
						$CallRate = $BC{"total"}{"Calls"}/ $BC{"total"}{"total"};

						#print "\t".keys %total_alleles;
						if (keys %total_alleles == 2){
		
							#Sort bases so p is the major allele and q is the minor allele
							my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
							#Major allele
							my $b1 = $bases[1];
							#Minor allele
							my $b2 = $bases[0];
			
							#Total number of samples
							$n_total = $BC{"total"}{"Calls"};
							#Major allele frequency in all samples
							$pAll = $BC{"total"}{$b1}/($BC{"total"}{"Calls"}*2);
							#Minor allele frequency in all samples
							$qAll = $BC{"total"}{$b2}/($BC{"total"}{"Calls"}*2);

							#Heterozygosity expected in all samples
							$HeAll = 2*($pAll * $qAll);

							#Heterozygosity observed in all samples
							if ($BC{"total"}{"Het"}){
								$HoAll = $BC{"total"}{"Het"}/($BC{"total"}{"Calls"}*2);
							}else{
								$HoAll = 0;
							}
							#Allele frequency of each allele in each population
							if ($BC{"1"}{$b1}){
								$p1 = $BC{"1"}{$b1}/($BC{"1"}{"Calls"}*2);
							}else{
								$p1 = 0;
							}
							if ($BC{"2"}{$b1}){
								$p2 = $BC{"2"}{$b1}/($BC{"2"}{"Calls"}*2);
							}else{
								$p2 = 0;
							}
							if ($BC{"1"}{$b2}){
								$q1 = $BC{"1"}{$b2}/($BC{"1"}{"Calls"}*2);
							}else{
								$q1 = 0;
							}
							if ($BC{"2"}{$b2}){
								$q2 = $BC{"2"}{$b2}/($BC{"2"}{"Calls"}*2);
							}else{
								$q2 = 0;
							}

							#Heterozygosity observed in each population
							if ($BC{"1"}{"Het"}){
								$Ho1 = $BC{"1"}{"Het"}/$BC{"1"}{"Calls"}
							}else{
								$Ho1 = 0;
							}
							if ($BC{"2"}{"Het"}){
								$Ho2 = $BC{"2"}{"Het"}/$BC{"2"}{"Calls"}
							}else{
								$Ho2 = 0;
							}

							#Amount of pairwise difference between population
							$dxy = (($p1 * $q2) + ($p2 * $q1));

							#Heterozygosity expected
							$He1 = 2*($p1 * $q1);
							$He2 = 2*($p2 * $q2);


							#Average expected heterozygosity in each population	
							$HsBar = (($He1+$He2)/2);
			
							#The difference in alleles frequency
							$freq_dif = abs($p1 - $p2);


							#Average sample size for populations
							$n_bar = ($BC{"total"}{"Calls"} / 2);
							#Sample size for population 1
							if ($BC{"1"}{"Calls"}){
								$n_1 = $BC{"1"}{"Calls"};
							}else{
								$n_1 = 0;
							}
							#Sample size for population 2
							if ($BC{"2"}{"Calls"}){
								$n_2 = $BC{"2"}{"Calls"};
							}else{
								$n_2 = 0;
							}
							#Average observed heterozygosity weighted by population (NEed to scale for sample size)
							$H_bar = ((($Ho1 * $n_1) + ($Ho2 * $n_2)) / $n_total);

							#Sigma squared. The sample variance of allele p frequencies over populations
							$sigma_squared = ((($n_1 * (($p1 - $pAll) ** 2)) / (($Npops - 1) * $n_bar)) + (($n_2 * (($p2 - $pAll) ** 2) / (($Npops - 1) * $n_bar))));
							#The squared coefficient of variation of sample sizes
							$n_c = ((($Npops * $n_bar) - ((($n_1 ** 2) / ($Npops * $n_bar)) + (($n_2 ** 2) / ($Npops * $n_bar)))) / ($Npops - 1));
							#Weir and Cockerham, the observed component of variance for between populations
							unless (($n_c eq 0) or ($n_bar <= 1)){
								$WC_a  = (($n_bar / $n_c) * ($sigma_squared - ((1 / ($n_bar - 1)) * (($pAll * $qAll) - ((($Npops - 1) / $Npops) * $sigma_squared) - (0.25 * $H_bar)))));
							}else{
								$WC_a = "NA";
							}
							#Weir and Cockerham, the observed component of variance for between individuals within a population
							unless ($n_bar <= 1){
								$WC_b = (($n_bar / ($n_bar - 1)) * (($pAll * $qAll) - ((($Npops - 1) / $Npops) * $sigma_squared) - (((2 * $n_bar) - 1) / (4 * $n_bar) * $H_bar)));
							}else{
								$WC_b = "NA";
							}
							#Weir and Cockerham, the observed component of variance for between gametes within individuals
							$WC_c = (0.5 * $H_bar);
							#Weir and Cockerham denominator in Fst calculation
							unless (($WC_a eq "NA") or ($WC_b eq "NA")){			
								$WC_denom = ($WC_a + $WC_b + $WC_c);
							}else{
								$WC_denom = "NA";
							}
							#Weir and Cockerham, Theta (Fst)
							unless ($WC_denom eq "NA"){
								$WC_fst = ($WC_a / $WC_denom);
							}else {
								$WC_fst = "NA";
							}
							#Diversity (pi) for single site.
							$pi = 2*($pAll * $qAll);
							$BiCount++;
		
						}elsif (keys %total_alleles eq 1){
							$pAll = 1;
							$qAll = 0;
							$HeAll = 0;
							$HoAll = 0;
							$p1 = 1;
							$q1 = 0;
							$p2 = 1;
							$q2 = 0;
							$Ho1 = 0;
							$Ho2 = 0;
							$He1 = 0;
							$He2 = 0;
							$dxy = 0;
							$HsBar = "NA";
							$sigma_squared = "NA";
							$H_bar = "0";
							#Average sample size for populations
							$n_bar = ($BC{"total"}{"Calls"} / 2);
							#Sample size for population 1
							if ($BC{"1"}{"Calls"}){
								$n_1 = $BC{"1"}{"Calls"};
							}else{
								$n_1 = 0;
							}
							if ($BC{"2"}{"Calls"}){
								$n_2 = $BC{"2"}{"Calls"};
							}else{
								$n_2 = 0;
							}
							#Total number of samples
							$n_total = $BC{"total"}{"Calls"};
							$WC_a = "0";
							$WC_b = "0";
							$WC_c = "0";
							$WC_denom = "0";
							$WC_fst = "0";
							$freq_dif = "0";
							$pi = "0";
						}
						elsif (keys %total_alleles eq 3){ #Need to account for three alleles in tri-allelic sites.
							$pAll = "NA";
							$qAll = "NA";
							$HeAll = "NA";
							$HoAll = "NA";
							$p1 = "NA";
							$q1 = "NA";
							$p2 = "NA";
							$q2 = "NA";
							$Ho1 = "NA";
							$Ho2 = "NA";
							$He1 = "NA";
							$He2 = "NA";
							$dxy = "NA"; ##Need to fix
							$HsBar = "NA";
							$n_total = $BC{"total"}{"Calls"}; 
							$sigma_squared = "NA"; 
							$H_bar = "NA";
							$n_1 = "NA";
							$n_2 = "NA";
							$n_bar = "NA";
							$WC_a = "NA"; ##
							$WC_b = "NA"; ##
							$WC_c = "NA"; ##
							$WC_denom = "NA"; ##
							$WC_fst = "NA"; ##
							$pi = "NA"; ###Need to fix
							$freq_dif = "NA";
						}elsif ((keys %total_alleles eq 4) or (keys %total_alleles eq 0)){ #If there are four alleles
							$pAll = "NA";
							$qAll = "NA";
							$HeAll = "NA";
							$HoAll = "NA";
							$p1 = "NA";
							$q1 = "NA";
							$p2 = "NA";
							$q2 = "NA";
							$Ho1 = "NA";
							$Ho2 = "NA";
							$He1 = "NA";
							$He2 = "NA";
							$dxy = "NA"; ##Need to fix
							$HsBar = "NA";
							$n_total = $BC{"total"}{"Calls"}; 
							$sigma_squared = "NA"; 
							$H_bar = "NA";
							$n_bar = "NA";
							$n_1 = "NA";
							$n_2 = "NA";
							$WC_a = "NA"; ##
							$WC_b = "NA"; ##
							$WC_c = "NA"; ##
							$WC_denom = "NA"; ##
							$WC_fst = "NA"; ##
							$pi = "NA"; ###Need to fix
							$freq_dif = "NA";
						}
						print "\t$n_1\t$n_2\t$n_total\t$dxy\t$WC_a\t$WC_denom\t$WC_fst\t$He1\t$He2\t$freq_dif\n";
					}
				}
			}
		}
	}	
}
