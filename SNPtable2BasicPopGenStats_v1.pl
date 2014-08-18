#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin'; #For GObox server
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
my $Npops = 2;
my $BiCount = 0;
my $TriCount= 0;
my $QuadCount = 0;
my $SingleTri =0;

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

my @tmp = keys %poplist;

if ($#tmp eq 2){
	print STDERR "This is really only designed for 2 pops\tyou have $#tmp\n";
}

my $pop1 = $tmp[0];
my $pop2 = $tmp[1];

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
        	print $a[0]."-"."$a[1]\t$a[0]";
		foreach my $i (1..($badcolumns-1)){
			print "\t$a[$i]";
		}
		print "\tpAll\tqAll\tHeAll\tHoAll\tCallRate\tp1\tq1\tp2\tq2\tHo1\tHo2\tHe1\tHe2\tFst\tlnRH\tprivate1\tprivate2\tprivatepresent1\tprivatepresent2\tFis1\tFis2\tshared\tdxy\tHsBar\tfixedDif\tn_1\tn_2\tH_bar\tsigma_squared\tWC_a\tWC_b\tWC_c\tWC_fst\tpi\n";
	}
	else{
		next if /^\s*$/;
		print $a[0]."-"."$a[1]\t$a[0]";
		foreach my $i (1..($badcolumns-1)){
			print "\t$a[$i]";
		}
		my %BC;
		my %BS;
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
					$BC{$samplepop{$i}}{"Calls"}++;
					
					if($bases[0] ne $bases[1]){
						$BC{"total"}{"Het"}++;
						$BC{$samplepop{$i}}{"Het"}++;
					}
				}
			} 
		}

		my $pAll;
		my $qAll;
		my $rAll;
		my $HeAll;
		my $HoAll;
		my $CallRate;
		my $p1;
		my $q1;
		my $r1;
		my $p2;
		my $q2;
		my $r2;
		my $Ho1;
		my $Ho2;
		my $He1 ;
		my $He2;
		my $Fst;
		my $lnRH;
		my $private1;
		my $private2;
		my $privatepresent1;
		my $privatepresent2;
		my $Fis1;
		my $Fis2;
		my $shared;
		my $dxy;
		my $HsBar;
		my $fixedDif;
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
			if ($BC{$pop1}{$b1}){
				$p1 = $BC{$pop1}{$b1}/($BC{$pop1}{"Calls"}*2);
			}else{
				$p1 = 0;
			}
			if ($BC{$pop2}{$b1}){
				$p2 = $BC{$pop2}{$b1}/($BC{$pop2}{"Calls"}*2);
			}else{
				$p2 = 0;
			}
			if ($BC{$pop1}{$b2}){
				$q1 = $BC{$pop1}{$b2}/($BC{$pop1}{"Calls"}*2);
			}else{
				$q1 = 0;
			}
			if ($BC{$pop2}{$b2}){
				$q2 = $BC{$pop2}{$b2}/($BC{$pop2}{"Calls"}*2);
			}else{
				$q2 = 0;
			}

			#Heterozygosity observed in each population
			if ($BC{$pop1}{"Het"}){
				$Ho1 = $BC{$pop1}{"Het"}/$BC{$pop1}{"Calls"}
			}else{
				$Ho1 = 0;
			}
			if ($BC{$pop2}{"Het"}){
				$Ho2 = $BC{$pop2}{"Het"}/$BC{$pop2}{"Calls"}
			}else{
				$Ho2 = 0;
			}

			#Amount of pairwise difference between population
			$dxy = (($p1 * $q2) + ($p2 * $q1));

			#Heterozygosity expected
			$He1 = 2*($p1 * $q1);
			$He2 = 2*($p2 * $q2);

			my $tmp = 1 / (1 - $He1);
			$tmp = $tmp * $tmp;
			my $tmpTOP = $tmp -1;
			
			$tmp = 1 / (1 - $He2);
			$tmp = $tmp * $tmp;
			my $tmpBOT = $tmp -1;
			if (($tmpTOP) and ($tmpBOT)){ 
				$lnRH = log( $tmpTOP / $tmpBOT );
			}else{
				$lnRH = "NA";
			}
			#print "$tmpTOP\t$tmpBOT\n\n";
			#$lnRH = "X";
			#Private alleles for each pop
			if (($p1 > 0) and ($p2 eq 0)){
				$private1 = $b1;
			}elsif(($q1 > 0) and ($q2 eq 0)){
				$private1 = $b2;
			}else{
				$private1 = "-";
			}
			if (($p2 > 0) and ($p1 eq 0)){
				$private2 = $b1;
			}elsif(($q2 > 0) and ($q1 eq 0)){
				$private2 = $b2;
			}else{
				$private2 = "-";
			}

			#Presence of private allele?
			if ($private1 eq "-"){
				$privatepresent1 = 0;
			}else{
				$privatepresent1 = 1;
			}
			if ($private2 eq "-"){
				$privatepresent2 = 0;
			}else{
				$privatepresent2 = 1;
			}

			#Presence of shared polymorphism?
			if (($p1 > 0) and ($p2 > 0) and ($q1 > 0) and ($q2 > 0)){
				$shared = 1;
			}else{
				$shared = 0;
			}	

			#Average expected heterozygosity in each population	
			$HsBar = (($He1+$He2)/2);

			#Fst 
			$Fst = ($HeAll - $HsBar)/ $HeAll;

			#Fis for each population 
			if ($He1 ne 0) {
				$Fis1 = ($He1 - $Ho1)/$He1 ;
			}else{
				$Fis1 = 0;
			}
			if ($He2 ne 0) {			
				$Fis2 = ($He2 - $Ho2)/$He2;
			}else{
				$Fis2 = 0;
			}
			#Proportion of third alleles in each population and total
			$rAll = 0;
			$r1 = 0;
			$r2 = 0;

			#Fixed differences?
			if (($p1 == 1) and ($q2 == 1)){
				$fixedDif = 1;
			}elsif (($p2 == 1) and ($q1 == 1)){
				$fixedDif = 1;
			}else{
				$fixedDif = 0;
			}

			#Average sample size for populations
			$n_bar = ($BC{"total"}{"Calls"} / 2);
			#Sample size for population 1
			if ($BC{$pop1}{"Calls"}){
				$n_1 = $BC{$pop1}{"Calls"};
			}else{
				$n_1 = 0;
			}
			#Sample size for population 2
			if ($BC{$pop2}{"Calls"}){
				$n_2 = $BC{$pop2}{"Calls"};
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
			$rAll = 0;
			$HeAll = 0;
			$HoAll = 0;
			$p1 = 1;
			$q1 = 0;
			$r1 = 0;
			$p2 = 1;
			$q2 = 0;
			$r2 = 0;
			$Ho1 = 0;
			$Ho2 = 0;
			$He1 = 0;
			$He2 = 0;
			$Fst = "NA";
			$lnRH = "NA";
			$private1 = "-";
			$private2 = "-";
			$privatepresent1 = 0;
			$privatepresent2 = 0;
			$Fis1 = "NA";
			$Fis2 = "NA";
			$shared = "NA";
			$dxy = "NA";
			$HsBar = "NA";
			$fixedDif = "NA";
			$sigma_squared = "NA";
			$H_bar = "0";
			#Average sample size for populations
			$n_bar = ($BC{"total"}{"Calls"} / 2);
			#Sample size for population 1
			if ($BC{$pop1}{"Calls"}){
				$n_1 = $BC{$pop1}{"Calls"};
			}else{
				$n_1 = 0;
			}
			if ($BC{$pop2}{"Calls"}){
				$n_2 = $BC{$pop2}{"Calls"};
			}else{
				$n_2 = 0;
			}
			#Total number of samples
			$n_total = $BC{"total"}{"Calls"};
			$WC_a = "0";
			$WC_b = "0";
			$WC_c = "0";
			$WC_denom = "0";
			$WC_fst = "NA";

			$pi = "0";
		}
		elsif (keys %total_alleles eq 3){ #Need to account for three alleles in tri-allelic sites.
                        $pAll = "NA";
                        $qAll = "NA";
                        $rAll = "NA";
                        $HeAll = "NA";
                        $HoAll = "NA";
                        $p1 = "NA";
                        $q1 = "NA";
                        $r1 = "NA";
                        $p2 = "NA";
                        $q2 = "NA";
                        $r2 = "NA";
                        $Ho1 = "NA";
                        $Ho2 = "NA";
                        $He1 = "NA";
                        $He2 = "NA";
                        $Fst = "NA";
                        $lnRH = "NA";
                        $private1 = "NA";
                        $private2 = "NA";
                        $privatepresent1 = "NA";
                        $privatepresent2 = "NA";
			$Fis1 = "NA";
			$Fis2 = "NA";
			$shared = "NA";
			$dxy = "NA"; ##Need to fix
			$HsBar = "NA";
			$fixedDif = "NA";
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
		}elsif ((keys %total_alleles eq 4) or (keys %total_alleles eq 0)){ #If there are four alleles
			$pAll = "NA";
			$qAll = "NA";
			$rAll = "NA";
			$HeAll = "NA";
			$HoAll = "NA";
			$p1 = "NA";
			$q1 = "NA";
			$r1 = "NA";
			$p2 = "NA";
			$q2 = "NA";
			$r2 = "NA";
			$Ho1 = "NA";
			$Ho2 = "NA";
			$He1 = "NA";
			$He2 = "NA";
			$Fst = "NA";
			$lnRH = "NA";
			$private1 = "NA";
			$private2 = "NA";
			$privatepresent1 = "NA";
			$privatepresent2 = "NA";
			$Fis1 = "NA";
			$Fis2 = "NA";
			$shared = "NA";
			$dxy = "NA"; ##Need to fix
			$HsBar = "NA";
			$fixedDif = "NA";
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
		}
		print "\t$pAll\t$qAll\t$HeAll\t$HoAll\t$CallRate\t$p1\t$q1\t$p2\t$q2\t$Ho1\t$Ho2\t$He1\t$He2\t$Fst\t$lnRH\t$private1\t$private2\t$privatepresent1\t$privatepresent2\t$Fis1\t$Fis2\t$shared\t$dxy\t$HsBar\t$fixedDif\t$n_1\t$n_2\t$H_bar\t$sigma_squared\t$WC_a\t$WC_b\t$WC_c\t$WC_fst\t$pi\n";
	}
}
close IN;

