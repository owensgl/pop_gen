#!/bin/perl
use warnings;
use strict;
use Math::CDF;
use Statistics::Basic qw(:all);
use lib '/home/owens/bin/pop_gen/'; #For GObox server
#This script filters sites based on parental allele frequency and only prints out sites that have two alleles and an allele frequency difference. It also only prints out the parental allele frequencies.
my $popfile = $ARGV[0]; #A population file.


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

my $current_chrom;
my $badcolumns = 12;
my @parent_number_list;
my @hyb_number_list;
my @good_number_list;
my %previous_likelihoodP1;
my %previous_likelihoodP2;
my $counter = 0;
my %previous_cm;
my $use_cm_transition = "TRUE";
my $cm_factor = 10;
my %previous_pos;
#print "#The transition change is $transition_change\n";
#print "#Use cm biased transitions = $use_cm_transition with factor = $cm_factor\n";
my $count_per_chrom;
while (<STDIN>){
	chomp;
	my @a = split(/\t/,$_);

	if ($. == 1){ #Load in sample names associated with column numbers, as well as population
		foreach my $i($badcolumns..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				$samplename{$i} = $a[$i];
				if (($samplepop{$i} eq "P1") or ($samplepop{$i} eq "P2")){
					push(@parent_number_list,$i);
				}else{
					push(@hyb_number_list, $i);
				}
				push(@good_number_list, $i);

			}
			if ($a[$i] =~ m/^sample_/){ #if its a replicate
				$samplename{$i} = $a[$i];
				$samplepop{$i} = "H";
				push (@good_number_list, $i);
				push (@hyb_number_list, $i);
			}
		}
		print "chrom\tbp\tcm\tp\tq\tp1\tq1\tp2\tq2";
	}else{
		next if /^\s*$/;
		$counter++;
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
		my $cm = $a[4];
		if ($cm eq "NA"){ #Skip cp, mt, rDNA
			next;
		}
		my %current_observationP1;
		my %current_observationP2;
		my %BC;
		my %total_alleles;
		if (($counter % 100000)== 0){
         	       print STDERR "Processing $chrom $pos...\n";
	        }
		
		foreach my $i(@good_number_list){ #Get parental allele frequencies
			if ($a[$i] ne "NN"){
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
		unless (($BC{"P1"}{"Calls"}) and ($BC{"P2"}{"Calls"})){
			next;
		}
		unless (($BC{"P1"}{"Calls"} >= 5) and ($BC{"P2"}{"Calls"} >= 5)){
			next;
		}
		if (keys %total_alleles == 2){ #There must be two alleles
			$count_per_chrom++;
			my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
			#Major allele
			my $b1 = $bases[1];
			#Minor allele
			my $b2 = $bases[0];
			unless (($BC{"P1"}{"Calls"}) and ($BC{"P2"}{"Calls"})){
				next;
			}unless (($BC{"P1"}{"Calls"} >= 5) and ($BC{"P2"}{"Calls"} >= 5)){
				next;
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
			if ($p1 eq $p2){
				next; #Skip sites with equal frequency from parents.
			}
			print "\n$chrom\t$pos\t$cm\t$b1\t$b2\t$p1\t$q1\t$p2\t$q2";
		}
	}
}
