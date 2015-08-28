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
my $fasta = $ARGV[1]; #Fasta file of each gene

my %genehash;
my %genestart;
my %geneend;
my %genechr;
my @genelist;
open FASTA, $fasta;
while (<FASTA>){
	chomp;
	if ($_ =~/^\>/g){
		my @a = split(/\ /, $_);
		my $begin = $a[2];
		$begin =~ s/begin=//g;
		my $end = $a[3];
		$end =~ s/end=//g;
		my $name = $a[0];
		$name =~ s/>//g;
		my $chr = $a[5];
		$chr =~ s/chr=//g;
		$genestart{$name} = ($begin - 200);
		$geneend{$name} = ($end + 200);
		$genechr{$name} = $chr;
		push (@genelist, $name);
	}
}

close FASTA;		 
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
my $current_gene;
my $min_sites = 1; #minimum number of sites used in a window
#Looping stats
my %likelihoodcount;
my $sitecount;
my %likelihood;

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
		print  "sample\tchrom\tstart\tend\tgene\tsites\tpercentP2\tlikelihood";
	}else{
		next if /^\s*$/;
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
		my $gene;
		foreach my $g (@genelist){
			if ($genechr{$g} eq $chrom){
				if ($genestart{$g} <= $pos){
					if ($geneend{$g} >= $pos){
						$gene = $g;
						goto FOUNDGENE; #If you found the gene, go on, or skip this line.
					}
				}
			}
		}
		next;
		FOUNDGENE:
		unless ($current_gene){
			$current_gene = $gene;
		}
		if ($current_gene ne $gene){
			foreach my $i($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($likelihoodcount{$i}){
							if ($likelihoodcount{$i} > $min_sites){
								foreach my $percentP2 (0..100){
									$percentP2 = $percentP2/100;
									print "\n$samplename{$i}\t$genechr{$current_gene}\t$genestart{$current_gene}\t$geneend{$current_gene}\t$current_gene\t$likelihoodcount{$i}\t$percentP2\t$likelihood{$i}{$percentP2}";
											
								}

							}else{	
								#print FINALOUT "\tNA:NA:NA:NA";
							}
						}else{
							# print FINALOUT "\tNA:NA:NA:NA";
						}
					}
				}
			}
			#Reset the variables
			undef %likelihood;
			undef %likelihoodcount;
			$sitecount = 0;
			$current_gene = $gene;
			
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
			foreach my $i ($badcolumns..$#a){
				if ($samplepop{$i}){
					if ($samplepop{$i} eq "H"){
						if ($BC{$i}{"Calls"}){
							$likelihoodcount{$i}++;
							foreach my $percentP2 (0..100){
                                                        	$percentP2 = $percentP2/100;
								my $percentP1 = 1 - $percentP2;
								if ($BC{$i}{$b1}){
									if ($BC{$i}{$b1} == 2){
										$likelihood{$i}{$percentP2} += log((($p1*$percentP1)+($p2*$percentP2)) * (($p1*$percentP1)+($p2*$percentP2)));
									}elsif ($BC{$i}{$b1} == 1){
										$likelihood{$i}{$percentP2} += log(2 * (($p1*$percentP1)+($p2*$percentP2)) * (($q1*$percentP1)+($q2*$percentP2)));
									}
								}else{
									$likelihood{$i}{$percentP2} += log((($q1*$percentP1)+($q2*$percentP2)) * (($q1*$percentP1)+($q2*$percentP2)));
								}
							}
						}
					}
				}
			}
		}
	}
	SKIP:
}
			
			
