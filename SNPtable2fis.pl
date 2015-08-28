#!/usr/bin/perl

use warnings;
use strict;
use List::MoreUtils qw(uniq);
#Usage: cat SNPTABLE.tab | SNPtable2fis.pl POPFILE.txt > outfis.txt
my $POP = $ARGV[0]; #population file

my %pophash;
my %samplepop;
open IN, $POP;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
	$pophash{$a[0]} = $a[1];
}
close IN;
my %fis;

my @tmppoplist =  sort values %pophash;
my @poplist = uniq @tmppoplist;
while (<STDIN>){
	chomp;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i (2..$#a){ #Get sample names for each column
			if ($pophash{$a[$i]}){
				$samplepop{$i} = $pophash{$a[$i]};
			}

        	}
		print "chr\tloc";
		foreach my $pop (@poplist){
			print "\t${pop}_Hobs\t${pop}_Hexp";
		}
	}else{
		my %BC;
		my %total_alleles;
		foreach my $i (2..$#a){
			$BC{"total"}{"total"}++;
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
		if (keys %total_alleles == 2){
			print "\n$a[0]\t$a[1]";
			#Sort bases so p is the major allele and q is the minor allele
			my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
			#Major allele
			my $b1 = $bases[1];
			#Minor allele
			my $b2 = $bases[0];
			
			foreach my $pop (@poplist){
				my $Hobs;
				my $Hexp;
				if($BC{$pop}{"Calls"}){
					my $p;
					if ($BC{$pop}{$b1}){
						$p = $BC{$pop}{$b1}/($BC{$pop}{"Calls"}*2);
					}else{
						$p = 0;
					}	
					my $q = 1 -$p;
					$Hexp = 2*$p*$q;
					if ($BC{$pop}{"Het"}){
						$Hobs = $BC{$pop}{"Het"}/$BC{$pop}{"Calls"};
					}else{
						$Hobs = 0;
					}
					$fis{$pop}{"Hobs"} += $Hobs;
					$fis{$pop}{"Hexp"} += $Hexp;
				}else{
					$Hobs = "NA";
					$Hexp = "NA";
				}
				print "\t$Hobs\t$Hexp";
			}
		}
	}
}

foreach my $pop (@poplist){
	my $fis = 1- ($fis{$pop}{"Hobs"}/$fis{$pop}{"Hexp"});
	print STDERR "$pop\t$fis\n";
}
