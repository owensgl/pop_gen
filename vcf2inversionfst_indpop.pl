#!/usr/bin/perl
#This calculates weir and cockerhams Fst between homozygous inversion genotypes. Both across all samples and within populations. 
#Usage: cat snptable.tab | perl SNPtable2fst.pl sampleinfo.txt populations.txt > snptable.fst.txt
use warnings;
use strict;

my $min_MAF = 0.05; #minimum total minor allele frequency
my $min_n = 1; #minimum number of samples called per population
my %pop;
my %poplist;
my %sample;
my %inv_genotype;
my %pop_geno;
my @ind_pops;
my $Npops = 2; #Don't change unless there are more than two groups being compared in fst.

unless (@ARGV == 2) {die;}
my $pop_file = $ARGV[0]; #List of samples linked to population 
my $inversion_file = $ARGV[1]; #List of all samples with genotypes for a single "inversion"

open POP, $pop_file;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
}
close POP;


open INV, $inversion_file;
while (<INV>){
	chomp;
	my @a = split (/\t/,$_);
	$inv_genotype{$a[0]} = $a[1];
	$pop_geno{$pop{$a[0]}}{$a[1]}++;
}
close INV;



#Determine how many populations will be used. 
foreach my $pop (sort keys %poplist){
	#Check to make sure it has more than 0 of each 0 and 2 genotypes;
	if(($pop_geno{$pop}{0}) and ($pop_geno{$pop}{2})){
		if (($pop_geno{$pop}{0} >= 1) and ($pop_geno{$pop}{2} >= 1)){
			push(@ind_pops,$pop);
		}
	}
}

while (<STDIN>){
	chomp;
	my @a = split (/\t/,$_);
	if ($_ =~ m/^##/){next;}
	if ($_ =~ m/^#/){
		foreach my $i (9..$#a){
			$sample{$i} = $a[$i]
		}
		print "chr\tpos\tpop\tN1\tN2\tFstNum\tFstDenom\tFst\tINV_0_allele\tINV_2_allele";
		next;
	}
	my $ref = $a[3];
	my $alt = $a[4];	
	my %genotypes;
	my %BC;
	my %BS;
	my %total_alleles;
	foreach my $i (9..$#a){
		if ($pop{$sample{$i}}){
			my @info = split(/:/,$a[$i]);
			unless (($info[0] eq './.') or ($info[0] eq '.')){
				$genotypes{$sample{$i}} = $info[0];
			}
		}
	}
	
	#Calculate Fst for all populations together.
	my %counts;
	my %Ho;
	my %total_counts;
	foreach my $sample (sort keys %pop){
		if (exists($inv_genotype{$sample}) and exists($genotypes{$sample})){
			if ($inv_genotype{$sample} != 1){
				$total_counts{$inv_genotype{$sample}}+=2;
				if ($genotypes{$sample} eq '0/0'){
					$counts{$inv_genotype{$sample}}{1}+=2;
				}elsif ($genotypes{$sample} eq '0/1'){
					$counts{$inv_genotype{$sample}}{1}+=1;
					$counts{$inv_genotype{$sample}}{2}+=1;
					$Ho{$inv_genotype{$sample}}+=1;
				}elsif ($genotypes{$sample} eq '1/1'){
					$counts{$inv_genotype{$sample}}{2}+=2;
				}
			}
		}
	}
	#Put in blanks.
	foreach my $i (1..2){
		foreach my $j (0,2){
			unless($counts{$j}{$i}){
				$counts{$j}{$i} = 0;
			}
			unless($Ho{$j}){
				$Ho{$j} = 0;
			}
		}
	}
	
	
	my $b1;
	my $b2;
	my $pAll;
	my $qAll;
	my $HeAll;
	my $HoAll;
	my $p1;
	my $q1;
	my $p2;
	my $q2;
	my $Ho1;
	my $Ho2;
	my $He1 ;
	my $He2;
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
	
	unless (($total_counts{0}) and ($total_counts{2})) {
		next;
	}
	
	unless(($total_counts{0} >= $min_n) and ($total_counts{2} >= $min_n)){
		next;
	}
	if (($counts{0}{1} == 0) and ($counts{2}{1} == 0)){next;}
	if (($counts{0}{2} == 0) and ($counts{2}{2} == 0)){next;}
	
	if (($counts{0}{1} + $counts{2}{1}) > ($counts{0}{2} + $counts{2}{2})){
		$b1 = 1;
		$b2 = 2;
	}else{
		$b1 = 2;
		$b2 = 1;
	}
	
	#Total number of samples
	$n_total = ($counts{0}{1} + $counts{2}{1} + $counts{0}{2} + $counts{2}{2})/2;
	#Major allele frequency in all samples
	$pAll = ($counts{0}{$b1} + $counts{2}{$b1}) / ($n_total * 2);
	#Minor allele frequency in all samples
	$qAll = 1 -$pAll;
	if ($qAll < $min_MAF){
		goto SKIP; #Skip line if minor allele freq is less than cut off
	}
	#Heterozygosity expected in all samples
	$HeAll = 2*($pAll * $qAll);
	
	#Sample size for populations
	$n_1 = ($counts{0}{$b1} + $counts{0}{$b2})/2;
	$n_2 = ($counts{2}{$b1} + $counts{2}{$b2})/2;
	
	#Allele frequency of each allele in each population
	$p1 = $counts{0}{$b1}/($n_1 * 2);
	$p2 = $counts{2}{$b1}/($n_2 * 2);
	$q1 = 1 - $p1;
	$q2 = 1 - $p2;
	#Heterozygosity observed in each population
	$Ho1 = ($Ho{0}*2) / ($counts{0}{$b1} + $counts{0}{$b2});
	$Ho2 = ($Ho{2}*2) / ($counts{2}{$b1} + $counts{2}{$b2});
	
	#Heterozygosity expected
	$He1 = 2*($p1 * $q1);
	$He2 = 2*($p2 * $q2);
	
	#Average expected heterozygosity in each population	
	$HsBar = (($He1+$He2)/2);
	
	#Average sample size for populations
	$n_bar = ($n_total / 2);
	
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
	
	#Find which allele is more common in each version of the inversion;
	my $inv_1_allele = $ref;
	if ($counts{0}{1} < $counts{0}{2}){
		$inv_1_allele = $alt;
	}
	my $inv_2_allele = $ref;
	if ($counts{2}{1} < $counts{2}{2}){
		$inv_2_allele = $alt;
	}
	print "\n";
	print "$a[0]\t$a[1]\tall_pops";
	print "\t$n_1\t$n_2\t$WC_a\t$WC_denom\t$WC_fst";
	print "\t$inv_1_allele\t$inv_2_allele";
	
	SKIP:

	#Now try it for all the individual populations;

	foreach my $pop (@ind_pops){
		#Undefine a bunch of variables so that they can be reused. 
		undef(%counts);
		undef(%Ho);
		undef(%total_counts);
		undef(%counts);
		undef(%Ho);
		undef(%total_counts);
		undef($b1);
		undef($b2);
		undef($pAll);
		undef($qAll);
		undef($HeAll);
		undef($HoAll);
		undef($p1);
		undef($q1);
		undef($p2);
		undef($q2);
		undef($Ho1);
		undef($Ho2);
		undef($He1 );
		undef($He2);
		undef($HsBar);
		undef($H_bar);
		undef($n_bar);
		undef($n_1);
		undef($n_2);
		undef($n_total);
		undef($sigma_squared);
		undef($n_c);
		undef($WC_a);
		undef($WC_b);
		undef($WC_c);
		undef($WC_denom);
		undef($WC_fst);
		
		foreach my $sample (sort keys %pop){
			if (exists($inv_genotype{$sample}) and exists($genotypes{$sample})){
				if ($pop{$sample} ne $pop){next;}
				if ($inv_genotype{$sample} != 1){
					$total_counts{$inv_genotype{$sample}}+=2;
					if ($genotypes{$sample} eq '0/0'){
						$counts{$inv_genotype{$sample}}{1}+=2;
					}elsif ($genotypes{$sample} eq '0/1'){
						$counts{$inv_genotype{$sample}}{1}+=1;
						$counts{$inv_genotype{$sample}}{2}+=1;
						$Ho{$inv_genotype{$sample}}+=1;
					}elsif ($genotypes{$sample} eq '1/1'){
						$counts{$inv_genotype{$sample}}{2}+=2;
					}
				}
			}
		}
		#Put in blanks.
		foreach my $i (1..2){
			foreach my $j (0,2){
				unless($counts{$j}{$i}){
					$counts{$j}{$i} = 0;
				}
				unless($Ho{$j}){
					$Ho{$j} = 0;
				}
			}
		}
		unless (($total_counts{0}) and ($total_counts{2})) {
			next;
		}
		
		unless(($total_counts{0} >= $min_n) and ($total_counts{2} >= $min_n)){
			next;
		}
		if (($counts{0}{1} == 0) and ($counts{2}{1} == 0)){next;}
		if (($counts{0}{2} == 0) and ($counts{2}{2} == 0)){next;}
		
		if (($counts{0}{1} + $counts{2}{1}) > ($counts{0}{2} + $counts{2}{2})){
			$b1 = 1;
			$b2 = 2;
		}else{
			$b1 = 2;
			$b2 = 1;
		}
		
		#Total number of samples
		$n_total = ($counts{0}{1} + $counts{2}{1} + $counts{0}{2} + $counts{2}{2})/2;
		#Major allele frequency in all samples
		$pAll = ($counts{0}{$b1} + $counts{2}{$b1}) / ($n_total * 2);
		#Minor allele frequency in all samples
		$qAll = 1 -$pAll;
		if ($qAll < $min_MAF){
			next; #Skip line if minor allele freq is less than cut off
		}
		#Heterozygosity expected in all samples
		$HeAll = 2*($pAll * $qAll);
		
		#Sample size for populations
		$n_1 = ($counts{0}{$b1} + $counts{0}{$b2})/2;
		$n_2 = ($counts{2}{$b1} + $counts{2}{$b2})/2;
		
		#Allele frequency of each allele in each population
		$p1 = $counts{0}{$b1}/($n_1 * 2);
		$p2 = $counts{2}{$b1}/($n_2 * 2);
		$q1 = 1 - $p1;
		$q2 = 1 - $p2;
		#Heterozygosity observed in each population
		$Ho1 = ($Ho{0}*2) / ($counts{0}{$b1} + $counts{0}{$b2});
		$Ho2 = ($Ho{2}*2) / ($counts{2}{$b1} + $counts{2}{$b2});
		
		#Heterozygosity expected
		$He1 = 2*($p1 * $q1);
		$He2 = 2*($p2 * $q2);
		
		#Average expected heterozygosity in each population	
		$HsBar = (($He1+$He2)/2);
		
		#Average sample size for populations
		$n_bar = ($n_total / 2);
		
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
		
		#Find which allele is more common in each version of the inversion;
		my $inv_1_allele = $ref;
		if ($counts{0}{1} < $counts{0}{2}){
			$inv_1_allele = $alt;
		}
		my $inv_2_allele = $ref;
		if ($counts{2}{1} < $counts{2}{2}){
			$inv_2_allele = $alt;
		}
		print "\n";
		print "$a[0]\t$a[1]\t$pop";
		print "\t$n_1\t$n_2\t$WC_a\t$WC_denom\t$WC_fst";
		print "\t$inv_1_allele\t$inv_2_allele";

	}	


}



