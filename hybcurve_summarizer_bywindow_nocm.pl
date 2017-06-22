#!/bin/perl
use warnings;
use strict;

#This script takes in the hybcurve likelihood outputs and summarizes each window for plotting.
#my $in = $ARGV[0];
my $current_chrom;
my $current_sample;
my $current_start;
my $current_end;
my $current_gene;
my $current_sites;
my $Lcutoff = 1.92;
my %hash;
#open IN, $in;
while (<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	if ($. == 1){
		print "sample\tchrom\tstart\tend\tmaxP\tlowP\thighP\twidth\tstatus";
	}else{
		my $sample = $a[0];
		my $chrom = $a[1];
		my $start = $a[2];
		my $end = $a[3];
		my $percentP2 = $a[4];
		my $L = $a[5];
		
		unless ($current_sample){
			$current_sample = $sample;
			$current_chrom = $chrom;
			$current_start = $start;
			$current_end = $end;
		}
		if (($sample eq $current_sample) and ($start eq $current_start)){
			$hash{$percentP2} = $L;
		}
		else{
			my $maxL;
			my $maxP;
			my $highP;
			my $lowP;
			foreach my $key (sort {$a<=>$b} keys %hash){
				unless ($maxL){
					$maxL = $hash{$key};
					$maxP = $key;
				}
				if ($maxL < $hash{$key}){
					$maxL = $hash{$key};
					$maxP = $key;
				}
			}
			my $lowPcounter;
			foreach my $key (sort {$a<=>$b} keys %hash){
				no warnings 'uninitialized';
				if ($hash{$key} >= ($maxL - $Lcutoff)){
					unless ($lowP){
						unless ($lowP eq "0"){
							$lowP = $key;
							$lowPcounter++;
						}
					}
				}else{
					if ($lowPcounter){
						unless ($highP){
							$highP = ($key - 0.01);
						}
					}
				}
			}
			unless ($highP){
				$highP = 1;
			}
			my $status;
			my $width = $highP - $lowP;
			if ($width > 0.6){
				$status = "NA"
			}
			elsif ($highP < 0.5){
				$status = "P1";
			}elsif ($lowP > 0.5){
				$status = "P2";
			}else{
				$status = "admixed";
			}
			print "\n$current_sample\t$current_chrom\t$current_start\t$current_end\t$maxP\t$lowP\t$highP\t$width\t$status";
			
			#Refresh variables
			$current_sample = $sample;
			$current_chrom = $chrom;
			$current_start = $start;
			$current_end = $end;
			undef(%hash);
			$hash{$percentP2} = $L;
		}
	}
}
		
