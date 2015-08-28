#!/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin/pop_gen/'; #For GObox server

unless (@ARGV == 3) {die;}
my $in = $ARGV[0]; #SNP table
my $pop = $ARGV[1]; #List of samples linked to population
my $groups = $ARGV[2]; #List of all populations with populations selected (1-2 = young, 3-4 old)
my @comparisons = ('1','3');
require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;

my %samplepop;
my %samplegroup;
my %pop;
my %poplist;
my %group;

my $perm_n = 1000;
open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
}
close POP;

open GROUP, $groups;
while (<GROUP>){
        chomp;
        my @a = split (/\t/,$_);
        $group{$a[0]} = $a[1];
}
close GROUP;

open IN, $in;
while (<IN>){
	chomp;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i ($badcolumns..$#a){ #Get sample names for each column
        		if ($pop{$a[$i]}){
        			$samplepop{$i} = $pop{$a[$i]};
				if ($group{$pop{$a[$i]}}){
					$samplegroup{$i} = $group{$pop{$a[$i]}};
				}
			}
        	}
		print "chr\tpos\tt1_p1\tt1_p2\tt2_p1\tt2_p2\tt1_fst\tt2_fst\tt1_freqdif\tt2_freqdif\tperm_mean_fst\tperm_mean_freqdif\t";
		print "fst_gt_pvalue\tfst_lt_pvalue\tfreqdif_gt_pvalue\tfreqdif_lt_pvalue";
	}else{
		next if /^\s*$/;
		my $chr = $a[0];
		my $pos = $a[1];
		my %BC;
		my %total_alleles;
		foreach my $i ($badcolumns..$#a){
			if ($samplegroup{$i}){
				unless (($a[$i] eq "NN")or($a[$i] eq "XX")){
					my @bases = split(//, $a[$i]);
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
		if (keys %total_alleles != 2){
			next;
		}else{
			my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
			#Major allele
			my $b1 = $bases[1];
			#Minor allele
			my $b2 = $bases[0];
			my %hash;
			my %sizehash;
			my %p_hash;
			foreach my $num (@comparisons){ #Calculate basic stats and then pass that to fst/freq diff calculator
				my $pop1 = $num;
				my $pop2 = $num+1;
				unless(($BC{$pop1}{"Calls"}) and ($BC{$pop2}{"Calls"})){
					goto SKIP;
				}
				unless(($BC{$pop1}{"Calls"}>=10) and ($BC{$pop2}{"Calls"} >=10)){
					goto SKIP;
				}
				my $n_1 = $BC{$pop1}{"Calls"};
				my $n_2 = $BC{$pop2}{"Calls"};
				$sizehash{$pop1} = $n_1;
				$sizehash{$pop2} = $n_2;
				my $Ho1;
				my $Ho2;
				if ($BC{$pop1}{"Het"}){
					$Ho1 = $BC{$pop1}{"Het"} / $n_1;
				}else{
					$Ho1 = 0;
				}
				if ($BC{$pop2}{"Het"}){
					$Ho2 = $BC{$pop2}{"Het"} / $n_2;
				}else{
					$Ho2 = 0;
				}

				my $p1;
				my $p2;
				if ($BC{$pop1}{$b1}){
					$p1 = $BC{$pop1}{$b1}/ ($n_1*2);
				}else{
					$p1 = 0;
				}
				if ($BC{$pop2}{$b1}){
					$p2 = $BC{$pop2}{$b1}/($n_2*2);
				}else{
					$p2 = 0;
				}
				$p_hash{$pop1} = $p1;
				$p_hash{$pop2} = $p2;
				my @results = &calc_stats($Ho1,$Ho2,$n_1,$n_2,$p1,$p2);
				$hash{$num} = \@results; #This array contains Fst [0], and freq_dif [1]
			}
			my @perm_fst;
			my @perm_freqdif;
			foreach my $perm_number (1..$perm_n){
				my $count1 = 0;
				my @pop1_data;
				until ($count1 eq $sizehash{3}){
					my $rand = int(rand($#a-($badcolumns-1))) + $badcolumns; #generate a random value corresponding to a genotype
					if ($samplegroup{$rand}){
						if ($samplegroup{$rand} == 1){
							if ($a[$rand] ne "NN"){
								$count1++;
								push(@pop1_data, $a[$rand]);
							}
						}
					}
				}
				my @pop1_summed = &sum_alleles(\@pop1_data, $b1);
				my $count2 = 0;
				my @pop2_data;
				until ($count2 eq $sizehash{4}){
					my $rand = int(rand($#a-($badcolumns-1))) + $badcolumns; #generate a random value corresponding to a genotype
					if ($samplegroup{$rand}){
						if ($samplegroup{$rand} == 2){
							if ($a[$rand] ne "NN"){
								$count2++;
								push(@pop2_data, $a[$rand]);
							}
						}
					}
				}
				my @pop2_summed = &sum_alleles(\@pop2_data, $b1);
				my @stats = &calc_stats($pop1_summed[0],$pop2_summed[0],$pop1_summed[1],$pop2_summed[1],$pop1_summed[2],$pop2_summed[2]);
				push(@perm_fst, $stats[0]);
				push(@perm_freqdif, $stats[1]);
#				print "$stats[1]\n";
			}
			my $average_fst = &average(@perm_fst);
			my @sorted_perm_fst = sort { $a <=> $b } @perm_fst;
			my $fst_counter = 0;
			foreach my $value (@sorted_perm_fst){
				if ($hash{3}[0] > $value){
					$fst_counter++;
				}else{
					goto NEXT1;
				}
			}
			NEXT1:
			my $largerPvaluefst = (($perm_n - $fst_counter)+1) / ($perm_n+1);
			my @reverse_sorted_perm_fst = sort { $b <=> $a } @perm_fst;                  
			my $reverse_fst_counter = 0;
                        foreach my $value (@reverse_sorted_perm_fst){
                                if ($hash{3}[0] < $value){
                                        $reverse_fst_counter++;
                                }else{
                                        goto NEXT2;
                                }
                        }
                        NEXT2:
			my $smallerPvaluefst = (($perm_n - $reverse_fst_counter)+1) / ($perm_n+1);
			my $average_freqdif = &average(@perm_freqdif);
			my @sorted_perm_freqdif = sort { $a <=> $b } @perm_freqdif;
			my $freqdif_counter = 0;
			foreach my $value (@sorted_perm_freqdif){
				if ($hash{3}[1] > $value){
					$freqdif_counter++;
				}else{
					goto NEXT3;
				}
			}
			NEXT3:
			my $largerPvaluefreqdif = (($perm_n - $freqdif_counter)+1) / ($perm_n+1);
                        my @reverse_sorted_perm_freqdif = sort { $b <=> $a } @perm_freqdif;
                        my $reverse_freqdif_counter = 0;
                        foreach my $value (@reverse_sorted_perm_freqdif){
                                if ($hash{3}[1] < $value){
                                        $reverse_freqdif_counter++;
                                }else{
                                        goto NEXT4;
                                }
                        }
                        NEXT4:			
			my $smallerPvaluefreqdif = (($perm_n - $reverse_freqdif_counter)+1) / ($perm_n+1);

			print "\n$chr\t$pos\t$p_hash{1}\t$p_hash{2}\t$p_hash{3}\t$p_hash{4}\t";
			print "$hash{1}[0]\t$hash{3}[0]\t$hash{1}[1]\t$hash{3}[1]\t";
			print "$average_fst\t$average_freqdif\t$largerPvaluefst\t$smallerPvaluefst\t";
			print "$largerPvaluefreqdif\t$smallerPvaluefreqdif";
#			exit;

		}
	}
	SKIP:
}


sub average {
	my @array = @_; # save the array passed to this function
	my $sum; # create a variable to hold the sum of the array's values
	foreach (@array) { $sum += $_; } # add each element of the array 
	# to the sum
	return $sum/@array; # divide sum by the number of elements in the
	# array to find the mean
}




sub sum_alleles{
	my ($data, $b1) = @_;
	my @data = @{$data};
	my $Ho = 0;
	my $p1 = 0;
	my $n_1 = 0;
	foreach my $pair (@data){
		$n_1++;
		my @tmp = split(//,$pair);
		if ($tmp[0] ne $tmp[1]){
			$Ho++;
		}
		foreach my $n (0..1){
			if ($tmp[$n] eq $b1){
				$p1++;
			}
		}
	}
	$Ho = $Ho / $n_1;
	$p1 = $p1 / ($n_1 *2);
	my @results = ($Ho, $n_1, $p1);
	return (@results);
}



sub calc_stats{
	my $Ho1 = shift;
	my $Ho2 = shift;
	my $n_1 = shift;
	my $n_2 = shift;
	my $p1 = shift;
	my $p2 = shift;
	my $q1 = 1 - $p1;
	my $q2 = 1 - $q1;
	my $n_total = $n_1 + $n_2;
	my $Npops = 2;
	my $n_bar = $n_total / $Npops;
	my $pAll = (($p1 * $n_1) + ($p2 * $n_2)) / $n_total;
	my $qAll = 1- $pAll;
	my $freq_dif;
	my $WC_fst;
	if (($pAll == 1) or ($qAll == 1)){
		$WC_fst = 0;
		$freq_dif =0;
		goto SKIPLINE;
	}
		
				#Average observed heterozygosity weighted by population (NEed to scale for sample size)
	my $H_bar = ((($Ho1 * $n_1) + ($Ho2 * $n_2)) / $n_total);

			#Sigma squared. The sample variance of allele p frequencies over populations
	my $sigma_squared = ((($n_1 * (($p1 - $pAll) ** 2)) / (($Npops - 1) * $n_bar)) + (($n_2 * (($p2 - $pAll) ** 2) / (($Npops - 1) * $n_bar))));
			#The squared coefficient of variation of sample sizes
	my $n_c = ((($Npops * $n_bar) - ((($n_1 ** 2) / ($Npops * $n_bar)) + (($n_2 ** 2) / ($Npops * $n_bar)))) / ($Npops - 1));
			#Weir and Cockerham, the observed component of variance for between populations
	my $WC_a;
	unless (($n_c eq 0) or ($n_bar <= 1)){
		$WC_a  = (($n_bar / $n_c) * ($sigma_squared - ((1 / ($n_bar - 1)) * (($pAll * $qAll) - ((($Npops - 1) / $Npops) * $sigma_squared) - (0.25 * $H_bar)))));
	}else{
		$WC_a = "NA";
	}
	#Weir and Cockerham, the observed component of variance for between individuals within a population
	my $WC_b;
	unless ($n_bar <= 1){
		$WC_b = (($n_bar / ($n_bar - 1)) * (($pAll * $qAll) - ((($Npops - 1) / $Npops) * $sigma_squared) - (((2 * $n_bar) - 1) / (4 * $n_bar) * $H_bar)));
	}else{
		$WC_b = "NA";
	}
			#Weir and Cockerham, the observed component of variance for between gametes within individuals
	my $WC_c = (0.5 * $H_bar);
			#Weir and Cockerham denominator in Fst calculation
	my $WC_denom;
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
	$freq_dif = abs($p1 - $p2);
	SKIPLINE:
	my @results;
	push (@results, $WC_fst);
	push (@results, $freq_dif);
	return (@results);
}
