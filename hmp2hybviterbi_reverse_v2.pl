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
my $parentfile = $ARGV[1]; #A list of parental allele frequencies

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
my %parent_hash;
open PARENTS, $parentfile;
while (<PARENTS>){
	chomp;
	my @a = split(/\t/,$_);
	if ($. == 1){next;}
	my $chrom = $a[0];
	my $pos = $a[1];
	my $cm = $a[2];
	my $b1 = $a[3];
	my $b2 = $a[4];
	my $p1 = $a[5];
	my $q1 = $a[6];
	my $p2 = $a[7];
	my $q2 = $a[8];
	$parent_hash{$chrom.$pos}{"b1"} = $b1;
	$parent_hash{$chrom.$pos}{"b2"} = $b2;
	$parent_hash{$chrom.$pos}{"p1"} = $p1;
	$parent_hash{$chrom.$pos}{"q1"} = $q1;
	$parent_hash{$chrom.$pos}{"p2"} = $p2;
	$parent_hash{$chrom.$pos}{"q2"} = $q2;
}
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
my $cm_factor = 1000;
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
		print "sample\tspecies\tchrom\tbp\tcm\tstate\tlikelihoodP1\tlikelihoodP2";
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
				$previous_cm{$i} = 0;
			}
		}
		if ($current_chrom ne $chrom){ #Starting a new chromosome
			#Reset the likelihoods.
			foreach my $i (@good_number_list){
				
				$previous_likelihoodP1{$i} = 1; #Start off with equal probabilities
                                $previous_likelihoodP2{$i} = 1;
                                $previous_cm{$i} = 0;
				$previous_pos{$i}= 0;
			}
			$current_chrom = $chrom;
			$count_per_chrom = 0;
		}
		unless ($parent_hash{$chrom.$pos}{"p1"}){
			next;
		}
		foreach my $i(@good_number_list){ #Get parental allele frequencies
			if ($a[$i] ne "NN"){
				my @bases = split(//, $a[$i]);

				$BC{$i}{$bases[0]}++;
                                $BC{$i}{$bases[1]}++;
				$BC{$i}{"Calls"}++;

				$BC{"total"}{"Calls"}++;
				$BC{$samplepop{$i}}{"Calls"}++;
			}
		}
		if ($parent_hash{$chrom.$pos}{"p1"}){ #There must be two alleles
			$count_per_chrom++;
			#Major allele
			my $b1 = $parent_hash{$chrom.$pos}{"b1"};
			#Minor allele
			my $b2 = $parent_hash{$chrom.$pos}{"b2"};
			my $p1 = $parent_hash{$chrom.$pos}{"p1"};
			my $p2 = $parent_hash{$chrom.$pos}{"p2"};
			my $q1 = $parent_hash{$chrom.$pos}{"q1"};
			my $q2 = $parent_hash{$chrom.$pos}{"q2"};
			#Allele frequency of each allele in each population
			foreach my $i(@good_number_list){
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
						my $cm_distance = $cm - $previous_cm{$i};
						$current_transition_change = ($cm_distance * $cm_factor) / 100;
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
