#!/usr/bin/perl
#This calculates weir and cockerhams Fst between homozygous inversion genotypes. It also calculates the observed heterozygosity in 0/1 inversion samples; 
use warnings;
use strict;

#my $min_MAF = 0.005; #minimum total minor allele frequency
my $min_n = 5; #minimum number of samples called per population
my $min_MAF = 0.05;
my $min_het_depth = 5; #Min reads for observing heterozygosity in 0/1 heterozygotes
my %pop;
my %poplist;
my %sample;
my %inv_genotype;
my %pop_geno;
my @ind_pops;
my $Npops = 2; #Don't change unless there are more than two groups being compared in fst.
my $print_XRQ_location = "TRUE"; #Decides whether to include columns for the XRQ location of each marker, stored in the info column (For when we used BWA to remap SNPs).
my $inversion_file = $ARGV[0]; #List of all samples with genotypes for a single "inversion"

open INV, $inversion_file;
while (<INV>){
	chomp;
	my @a = split (/\t/,$_);
	$inv_genotype{$a[0]} = $a[1];
}
close INV;

my @samples;


while (<STDIN>){
	chomp;
	my @a = split (/\t/,$_);
	if ($_ =~ m/^##/){next;}
	if ($_ =~ m/^#/){
		foreach my $i (9..$#a){
			$sample{$i} = $a[$i];
			push (@samples, $a[$i]);
		}
		if ($print_XRQ_location eq "TRUE"){
			print "XRQchr\tXRQpos\t";
		}
		print "chr\tpos\tN1\tN2\tFstNum\tFstDenom\tFst\tfreq_dif\tINV_0_allele\tINV_2_allele\thet_het_observed\thets_sampled";
		next;
	}
	my $xrq_chr;
	my $xrq_pos;
	if ($print_XRQ_location eq "TRUE"){
		my @infos = split(/\;/,$a[7]);
		my @xrq_tag = split(/\./,$infos[$#infos]);
		$xrq_chr = $xrq_tag[0];
		$xrq_chr =~ s/XRQ=//g;
		$xrq_pos = $xrq_tag[1];
	}
	my $ref = $a[3];
	my $alt = $a[4];	
	my %genotypes;
	my %BC;
	my %BS;
	my %total_alleles;
	my %depth;
	foreach my $i (9..$#a){
		my @info = split(/:/,$a[$i]);
		unless (($info[0] eq './.') or ($info[0] eq '.')){
			$genotypes{$sample{$i}} = $info[0];
			$depth{$sample{$i}} = $info[2];
		}
	}
	
	#Calculate Fst for all populations together.
	my %counts;
	my %Ho;
	my %total_counts;
	my $het_samples = 0;
	my $het_het_counts = 0;
	foreach my $sample (@samples){
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
			}else{
				#Counts heterozygosity in inversion heterozygotes
				if ($depth{$sample} < $min_het_depth){next;}
				$het_samples++;
				if ($genotypes{$sample} eq '0/1'){
					$het_het_counts++;
				}
			}
		}
	}
	my $het_het_percent = "NA";
	if ($het_samples > 0){
		$het_het_percent = $het_het_counts / $het_samples;
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
	my $freq_dif = abs($p1 - $p2);
	print "\n";
	if ($print_XRQ_location eq "TRUE"){
		print "$xrq_chr\t$xrq_pos\t";
	}	
	print "$a[0]\t$a[1]";
	print "\t$n_1\t$n_2\t$WC_a\t$WC_denom\t$WC_fst";
	print "\t$freq_dif\t$inv_1_allele\t$inv_2_allele";
	print "\t$het_het_percent\t$het_samples";
	
	SKIP:

	#Now try it for all the individual populations;


}



