#!/usr/bin/perl
#This calculates weir and cockerhams Fst between two populations
#Usage: cat snptable.vcf | perl SNPtable2fst.pl sampleinfo.txt populations.txt > snptable.fst.txt
use warnings;
use strict;

my $min_MAF = 0.05; #minimum total minor allele frequency
my $min_n = 1; #minimum number of samples called per population, requires more than min_n
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

unless (@ARGV == 2) {die;}
my $pop = $ARGV[0]; #List of samples linked to population 
my $groups = $ARGV[1]; #List of all populations with populations selected (1 and 2)

open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
}
close POP;

my %group;
open GROUP, $groups;
while (<GROUP>){
	chomp;
	my @a = split (/\t/,$_);
	$group{$a[0]} = $a[1];
}
close GROUP;

my $lastsample;
my %samplegroup;
while (<STDIN>){
	chomp;
	my @a = split (/\t/,$_);
	if ($_ =~ m/^##/){next;}
	if ($_ =~ m/^#/){
		foreach my $i (9..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				if ($group{$pop{$a[$i]}}){
					$samplegroup{$i} = $group{$pop{$a[$i]}};
				}
			}
		}
		print "chr\tpos\tref\talt";
		print  "\tN1\tN2\tNTotal\tFstNum\tFstDenom\tFst\tHexp1\tHexp2\tf1\tf2";
		next;
	}
	my $ref = $a[3];
	my $alt = $a[4];	
	my %BC;
	my %BS;
	my %total_alleles;
	foreach my $i (9..$#a){
		if ($samplegroup{$i}){
			my @info = split(/:/,$a[$i]);
			unless ($info[0] eq './.'){
				my @bases = split(/\//, $info[0]);
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
	my $allele_1;
	my $allele_2;
	
	unless ($BC{"total"}{"Calls"}){
		goto SKIP;
	}
	
	#print "\t".keys %total_alleles;
	unless (($BC{"1"}{"Calls"}) and ($BC{"2"}{"Calls"})){
		goto SKIP;
	}unless(($BC{"1"}{"Calls"} >= $min_n) and ($BC{"2"}{"Calls"} >= $min_n)){
		goto SKIP;
		
	}elsif (keys %total_alleles == 2){
		
		#Sort bases so p is the major allele and q is the minor allele
		my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
		#Major allele
		my $b1 = $bases[1];
		#Minor allele
		my $b2 = $bases[0];
		if ($b1 == 0){
			$allele_1 = $ref;
			$allele_2 = $alt;
		}elsif ($b1 == 1){
			$allele_1 = $alt;
			$allele_2 = $ref;
		}
		#Total number of samples
		$n_total = $BC{"total"}{"Calls"};
		#Major allele frequency in all samples
		$pAll = $BC{"total"}{$b1}/($BC{"total"}{"Calls"}*2);
		#Minor allele frequency in all samples
		$qAll = $BC{"total"}{$b2}/($BC{"total"}{"Calls"}*2);
		if ($qAll < $min_MAF){
			goto SKIP; #Skip line if minor allele freq is less than cut off
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
		
	}elsif (keys %total_alleles eq 1){
		goto SKIP;
	}
	elsif (keys %total_alleles eq 3){ 
		goto SKIP;
	}elsif ((keys %total_alleles eq 4) or (keys %total_alleles eq 0)){ #If there are four alleles
	goto SKIP;
	}
	print "\n";
	print "$a[0]\t$a[1]\t$ref\t$alt";
	print "\t$n_1\t$n_2\t$n_total\t$WC_a\t$WC_denom\t$WC_fst\t$He1\t$He2";
	if ($allele_1 = $ref){
		print "\t$p1\t$p2";
	}else{
		print "\t$q1\t$q2";
	}
	SKIP:
}


