#!/bin/perl
use warnings;
use strict;
use List::Util 'shuffle';
#This script prints out each site where a snp is unique to one parent, or where it is unique to the hybrid.
#Requires parents to have 5 individuals sampled.
my $samplefile = $ARGV[0];
my %printlist;
my %specieslist;
open SAMPLEFILE, $samplefile;
while(<SAMPLEFILE>){
	chomp;
  	my @a = split(/\t/,$_);
	$printlist{$a[0]}++;
  	$specieslist{$a[0]} = $a[1]; #should have P1, P2 and any other name for hybrids
}
close SAMPLEFILE;

my $counter;
my %samplelist;
my @good_number_list;
my %good_number_hash;
my %species;
my $max_count = 20;
my @hybrid_species= qw(Ano Des Par);
while(<STDIN>){
  	$counter++;
	chomp;
	my $line = "$_";
	my @a = split(/\t/,$line);
 	if ($. == 1){
    		foreach my $i (2..$#a){
      			$samplelist{$i} = $a[$i];
      			if ($specieslist{$a[$i]}){
        			$species{$i} = $specieslist{$a[$i]};
        			push(@good_number_list, $i);
				$good_number_hash{$i}++;
      			}
    		}
		print "chrom\tpos\tvalue";
    		next;
  	}
	unless ($a[3]){
		next;
	}
  	my $pos = $a[1];
	my $chrom = $a[0];
 	my %P1alleles;
  	my %P2alleles;
  	my $P1count = 0;
  	my $P2count = 0;
	my %species_count;
  	if (($counter % 100000)== 0){
		print STDERR "Processing $chrom $pos...\n";
	}
  	foreach my $i (keys %good_number_hash){ #Load up parental alleles
    		if ($a[$i] ne "NN"){
      			if ($species{$i} eq "P1"){
			        $P1count++;
		      	}elsif($species{$i} eq "P2"){
			        $P2count++;
		      	}else{
				$species_count{$species{$i}}++;
			}
    		}
  	}
  	unless(($P1count >=5) and ($P2count >= 5)){
    		next;
  	}
	foreach my $key (@hybrid_species){ #Make sure each hybrid species has atleast two genotypes
		unless($species_count{$key}){
			goto SKIP;
		}elsif($species_count{$key} < 2){
			goto SKIP;
		}
	}
	my $min_count; #Samples are subset down to the lowest sample size of either parent
	if ($P1count < $P2count){
		$min_count = $P1count;
	}else{
		$min_count = $P2count;
	}
	if ($min_count > $max_count){
		$min_count = $max_count;
	}
        my $P1count2 = 0;
        my $P2count2 = 0;
	my %parental_alleles;
	
        foreach my $j ( shuffle keys %good_number_hash){ #Load up parental alleles
                if ($a[$j] ne "NN"){
                	if ((($species{$j} eq "P1") and ($P1count2 < $min_count)) or (($species{$j} eq "P2") and ($P2count2 < $min_count))){
                        	if ($species{$j} eq "P1"){
                                	$P1count2++;
                                }else{
                                	$P2count2++;
                                }
                               	my @parentalbases = split(//,$a[$j]);
				$parental_alleles{$parentalbases[0]}++;
				$parental_alleles{$parentalbases[1]}++;
                        }
		}
	}
	#
	my %hybrid_alleles;
	my %hybrid_count;
	foreach my $i ( shuffle keys %good_number_hash){
		if ($a[$i] ne "NN"){
			if (($species{$i} ne "P1") and ($species{$i} ne "P2")){
				if ($hybrid_count{$species{$i}}){ #This only selects two per hybrid species.
					if($hybrid_count{$species{$i}} == 2){ next;}
				}
				$hybrid_count{$species{$i}}++;
				my @bases = split(//,$a[$i]);
				foreach my $n(0..1){
					unless ($parental_alleles{$bases[$n]}){
						$hybrid_alleles{$bases[$n]}{$species{$i}}++;
					}
				}
			}
		}
	}
	unless(%hybrid_alleles){ #Go to next line if there are no new alleles in any of the hybrids
		next;
	}
	my $printvalue;
	foreach my $base (sort keys %hybrid_alleles){
		if ($hybrid_alleles{$base}{"Ano"}){
			$printvalue .= 1;
		}
		if ($hybrid_alleles{$base}{"Des"}){	
			$printvalue .= 2;
		}
		if ($hybrid_alleles{$base}{"Par"}){
			$printvalue .=3;
		}
		print "\n$chrom\t$pos\t$printvalue";
		undef($printvalue);
	}
	SKIP:
}
