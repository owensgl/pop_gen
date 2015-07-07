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

my $linecount;
my $min_MAF = 0.05; #minimum total minor allele frequency
my $min_n = 1; #minimum number of samples called per population
my $max_Hobs = 0.6; #Maximum observed heterozygosity
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

unless (@ARGV == 4) {die;}
my $in = $ARGV[0]; #SNP table
my $pop = $ARGV[1]; #List of samples linked to population 
my $groups = $ARGV[2]; #List of all populations with populations selected (1 and 2)
my $out = $ARGV[3];
require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;
my $datestring = localtime();
open README, "> $out.outflank.readme.txt";
print README "This script was run on $datestring. It loaded the SNP table $in.\n";
print README "The population file was named $pop, and contained:\n";
open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
	print README "$_\n";
}
close POP;

print README "The group file was named $groups, and contained:\n";
my %group;
open GROUP, $groups;
while (<GROUP>){
        chomp;
        my @a = split (/\t/,$_);
        $group{$a[0]} = $a[1];
	print README "$_\n";
}
close GROUP;

open OUT, "> $out.outflankin.txt";

print README "The minimum minor allele frequency was >= $min_MAF.\n";
print README "The minimum number of samples per population was $min_n.\n";
print README "The maximum observed heterozygosity was $max_Hobs.\n";
my %samplegroup;
open IN, $in;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i ($badcolumns..$#a){ #Get sample names for each column
        		if ($pop{$a[$i]}){
        			$samplepop{$i} = $pop{$a[$i]};
				if ($group{$pop{$a[$i]}}){
					$samplegroup{$i} = $group{$pop{$a[$i]}};
				}
        		}
        	}
        	print OUT "LocusName\t$a[0]\t$a[1]";
		print OUT "\tN1\tN2\tNTotal\tT1\tT2\tFST\tT1NoCorr\tT2NoCorr\tFSTNoCorr\tHe";
	}
	else{
		next if /^\s*$/;
		my %BC;
		my %BS;
		my %total_alleles;
		foreach my $i ($badcolumns..$#a){
			if ($samplegroup{$i}){
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
					$BC{$samplegroup{$i}}{$bases[0]}++;
		 			$BC{$samplegroup{$i}}{$bases[1]}++;

					$BC{"total"}{"Calls"}++;
					$BC{$samplegroup{$i}}{"Calls"}++;
					
					if($bases[0] ne $bases[1]){
						$BC{"total"}{"Het"}++;
						$BC{$samplegroup{$i}}{"Het"}++;
					}
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
		my $WCun_a;
		my $WCun_b;
		my $WCun_c;
		my $WCun_denom;
		my $WCun_fst;

		unless ($BC{"total"}{"Calls"}){
			$BC{"total"}{"Calls"} = 0;
		}
		
		$CallRate = $BC{"total"}{"Calls"}/ $BC{"total"}{"total"};

		#print "\t".keys %total_alleles;
		unless (($BC{"1"}{"Calls"}) and ($BC{"2"}{"Calls"})){
			goto SKIP;
		}elsif (keys %total_alleles == 2){
		
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
			if ($qAll < $min_MAF){
				goto SKIP; #Skip line if minor allele freq is less than cut off
			}
			if ($BC{"2"}{"Calls"} < $min_n){ #Skip line if less than minimum number of samples sequenced
				goto SKIP;
			}elsif ($BC{"1"}{"Calls"} < $min_n){
				goto SKIP;
			}
			#Heterozygosity expected in all samples
			$HeAll = 2*($pAll * $qAll);

			#Heterozygosity observed in all samples
			if ($BC{"total"}{"Het"}){
				$HoAll = $BC{"total"}{"Het"}/($BC{"total"}{"Calls"}*2);
			}else{
				$HoAll = 0;
			}
			if ($HoAll > $max_Hobs){
				goto SKIP;
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
			#Uncorrect WC Fst for outflank
			unless ($n_c eq 0){
				$WCun_a = (($n_bar / $n_c) * ($sigma_squared));
			}else{
				$WCun_a = "NA";
			}
			$WCun_b = (($pAll * $qAll) - ((($Npops - 1) / $Npops) * $sigma_squared) - (((2 * $n_bar)) / (4 * $n_bar) * $H_bar));
			$WCun_c = (0.5 * $H_bar);
			
			unless ($WCun_a eq "NA"){
				$WCun_denom = ($WCun_a + $WCun_b + $WCun_c);
			}else{
				$WCun_denom = "NA";
			}
			
			unless ($WCun_denom eq "NA"){
				$WCun_fst = ($WCun_a / $WCun_denom);
			}else{
				$WCun_fst = "NA";
			}
		
		}elsif (keys %total_alleles eq 1){
			goto SKIP;
		}
		elsif (keys %total_alleles eq 3){ #Need to account for three alleles in tri-allelic sites.
			goto SKIP;
		}elsif ((keys %total_alleles eq 4) or (keys %total_alleles eq 0)){ #If there are four alleles
			goto SKIP;
		}
		$linecount++;
		print OUT "\n";
                print OUT $a[0]."-"."$a[1]\t$a[0]\t$a[1]";
		print OUT "\t$n_1\t$n_2\t$n_total\t$WC_a\t$WC_denom\t$WC_fst\t$WCun_a\t$WCun_denom\t$WCun_fst\t$HeAll";
	}
	SKIP:
}
close IN;
print README "There were $linecount SNPs printed";
close README;
