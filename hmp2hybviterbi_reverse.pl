#!/bin/perl
use warnings;
use strict;
use Math::CDF;
use Statistics::Basic qw(:all);
use lib '/home/owens/bin/pop_gen/'; #For GObox server
#This script uses the viterbi algorithm to pick states for a hybrid based on parental allele frequencies
#This site calculates the probability at each site, based on it's observations and the previous site, it then chooses the state of the immediately prior site.
#It uses a transition probability that is multiplied by the distance in cm (with a multiplication factor). So the further things are the easier it is to transition between states.  
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
				$samplepop{$i} = "H";
                                $samplename{$i} = $a[$i];
                                push (@good_number_list, $i);
                                push (@hyb_number_list, $i);
                        }
		}
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
		unless($current_chrom){
			$current_chrom = $chrom;
			$count_per_chrom = 0;
			foreach my $i (@good_number_list){
				$previous_likelihoodP1{$i} = 1; #Start off with equal probabilities
				$previous_likelihoodP2{$i} = 1;
				$previous_cm{$i} = $cm + 0.00000001;
			}
		}
		if ($current_chrom ne $chrom){ #Starting a new chromosome
			#Reset the likelihoods.
			foreach my $i (@good_number_list){
				
				$previous_likelihoodP1{$i} = 1; #Start off with equal probabilities
                                $previous_likelihoodP2{$i} = 1;
                                $previous_cm{$i} = $cm + 0.00000001;
				$previous_pos{$i}= 0;
			}
			$current_chrom = $chrom;
			$count_per_chrom = 0;
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
			foreach my $i(reverse(@good_number_list)){
				if ($BC{$i}{"Calls"}){
					if ($BC{$i}{$b1}){
						if ($BC{$i}{$b1} == 2){
							$current_observationP1{$i} = $p1 * $p1;
							$current_observationP2{$i} = $p2 * $p2;

						}elsif ($BC{$i}{$b1} == 1){
							$current_observationP1{$i} = $p1 * $q1 * 2;
							$current_observationP2{$i} = $p2 * $q2 * 2;
						}
					}else{
						$current_observationP1{$i} = $q1 * $q1;
						$current_observationP2{$i} = $q2 * $q2;
					}
					#Calculate transition probability based on cm distance;
					my $current_transition_same;
					my $current_transition_change;
					if ($use_cm_transition eq "TRUE"){
						my $cm_distance = $previous_cm{$i} - $cm;
						$current_transition_change =  ($cm_distance * $cm_factor)/100;
						if ($current_transition_change > 0.5){
							$current_transition_change = 0.5;
						}
						$current_transition_same = 1 - $current_transition_change;
					}else{
						$current_transition_change = 0.5;
						$current_transition_same = 0.5;
					}
					#Calculate the probability of the previous state, and the four possible ways it can go
					my $p1_p1 = $previous_likelihoodP1{$i} * $current_transition_same * $current_observationP1{$i};
					my $p1_p2 = $previous_likelihoodP1{$i} * $current_transition_change * $current_observationP2{$i};
					my $p2_p2 = $previous_likelihoodP2{$i} * $current_transition_same * $current_observationP2{$i};
					my $p2_p1 = $previous_likelihoodP2{$i} * $current_transition_change * $current_observationP1{$i};
					my %state_hash;
					$state_hash{"p1_p1"} = $p1_p1;
					$state_hash{"p1_p2"} = $p1_p2;
					$state_hash{"p2_p2"} = $p2_p2;
					$state_hash{"p2_p1"} = $p2_p1;
					my @state = (sort { $state_hash{$b} <=> $state_hash{$a} } keys %state_hash); #Sort the states to get the most probable one
					if (($state[0] eq "p1_p1") or ($state[0] eq "p1_p2")){
						#The last state is P1
						#Normalize likelihood so it sums to 1.
						my $norm_factor = 1 / ($p1_p1 + $p1_p2);
						$previous_likelihoodP1{$i} = $p1_p1 * $norm_factor;
						$previous_likelihoodP2{$i} = $p1_p2 * $norm_factor;
						print "\n$samplename{$i}\t$samplepop{$i}\t$chrom\t$pos\t$cm\t$previous_likelihoodP1{$i}\t$previous_likelihoodP2{$i}";
					}else{
						#The last state is P2
                                                #Normalize likelihood so it sums to 1.
                                                my $norm_factor = 1 / ($p2_p1 + $p2_p2);
						$previous_likelihoodP1{$i} = $p2_p1 * $norm_factor;
						$previous_likelihoodP2{$i} = $p2_p2 * $norm_factor;
#The last state is P2
                                                print "\n$samplename{$i}\t$samplepop{$i}\t$chrom\t$pos\t$cm\t$previous_likelihoodP1{$i}\t$previous_likelihoodP2{$i}";
					}
					$previous_cm{$i} = $cm;
					$previous_pos{$i} = $pos;
				}
			}
		}
	}
}
print "\nsample\tspecies\tchrom\tbp\tcm\tstate\tlikelihoodP1\tlikelihoodP2\n";
