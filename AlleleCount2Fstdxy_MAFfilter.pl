#!/bin/perl

use warnings;
use strict;

my $in = $ARGV[0];
my $minMAF = 0.05;
open IN, $in;

while (<IN>){
	chomp;
	if ($. == 1){
		print "CHROM-POS\tCHROM\tPOS\tN1\tN2\tNTotal\tDxy\tFstNum\tFstDenom\tFst\tHexp1\tHexp2\tFreqDif";
	}else{
		my @a = split(/\t/,$_);
		my $chrom = $a[0];
		my $n1;
		my $n2;
		my %pop1;
		my %pop2;
		my %alleles;
		my $pos = $a[1];
		my @nuc_cols = (4, 5, 10, 12);
		foreach my $i (3,5){
			unless ($a[$i] eq "."){
				$pop1{$a[$i]} = $a[$i+1];
				$n1 += $a[$i+1];
				$alleles{$a[$i]} += $a[$i+1];
			}
		}foreach my $i (10,12){
			unless ($a[$i] eq "."){
				$pop2{$a[$i]} = $a[$i+1];
				$n2 += $a[$i+1];
				$alleles{$a[$i]} += $a[$i+1];
			}
		}
		unless (($n1) and ($n2)){
			print "Don't have data for both pops\n";
			next;
		}
		my $count1 = ($n1/2);
		my $count2 = ($n2/2);
		my $countAll = $count1 + $count2;
		my $dxy;
		my $N;
		my $D;
		my $F;
		my $Hexp1;
		my $Hexp2;
		my $FreqDif;
		my $p1;
		my $p2;
		my $q1;
		my $q2;	
		
		if (keys %alleles == 1){
			$dxy = "0";
			$N = "0";
			$D = "0";
			$F = "Inf";
			$Hexp1 = "0";
			$Hexp2 = "0";
			$FreqDif = "0";
		}elsif (keys %alleles == 2){
			my @bases = sort { $alleles{$a} <=> $alleles{$b} } keys %alleles ;
			#Major allele
			my $b1 = $bases[1];
			#Minor allele
			my $b2 = $bases[0];
			if ($pop1{$b1}){
					$p1 = $pop1{$b1} / $n1; 
			}else{
					$p1 = 0;
			}
			$q1 = (1-$p1);
			if ($pop2{$b1}){
					$p2 = $pop2{$b1} / $n2; 
			}else{
					$p2 = 0;
			}
			$q2 = (1-$p2);
			#Check MAF
			my $maf = ($q1 + $q2)/2;
			if ($maf < $minMAF){
				next;
			}
			$dxy = (($p1 *$q2) + ($p2 * $q1));
			$N = ($p1 * ($q2 - $q1)) + ($p2 * ($q1 - $q2));
			$D = ($p1 * $q2) + ($q1 * $p2);
			if ($D > 0){
				$F = $N/$D;
			}else{
				$F = "Inf";
			}
			$Hexp1 = ($p1 * $q1);
			$Hexp2 = ($p2 * $q2);
			$FreqDif = abs($p1 - $p2);
		}else{
		#	print "Too many alleles\t";
			next;
		}
		print "\n";
		print "$chrom-$pos\t$chrom\t$pos\t$count1\t$count2\t$countAll\t$dxy\t$N\t$D\t$F\t$Hexp1\t$Hexp2\t$FreqDif";
	}
}

close IN;
