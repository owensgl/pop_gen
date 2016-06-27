#!/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);
#Calculates genetic distance between all individuals (includes species ID), using genpofad distance. Allows no missing data. 
my $samplefile = $ARGV[0];

my %printlist;
my %specieslist;
my @species_array;
open SAMPLEFILE, $samplefile;
while(<SAMPLEFILE>){
	chomp;
	my @a = split(/\t/,$_);
	$printlist{$a[0]}++;
	$specieslist{$a[0]} = $a[1]; #should have real species names, don't include species you don't want compared.
	push(@species_array,$a[1]);
}
close SAMPLEFILE;

my @tmp = uniq @species_array;
@species_array = @tmp;
my %species;
my %samplelist;
my $counter;
my @good_number_list;
my %data;
my $sitecounter;
my %totalsites;
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
			}
		}
		next;
	}
	unless ($a[2]){
		next;
	}
  	my $pos = $a[1];
	my $chrom = $a[0];
	if (($counter % 1000)== 0){
		print STDERR "Processing $chrom $pos...\n";
	}
	my %species_counter;
	my %geno_hash;
	foreach my $i (@good_number_list){
		if ($a[$i] ne "NN"){
			$species_counter{$species{$i}}++;
			my @bases = split(//,$a[$i]);
			$geno_hash{$i}{0} = $bases[0];
			$geno_hash{$i}{1} = $bases[1];
		}
	}
	#This part is currently useless
#	foreach my $species (@species_array){ #Make sure there are atleast two samples genotyped per species
#		unless ($species_counter{$species}){
			#print "$species has no data\n";
#			goto NEXT_LINE;
#		}
#		unless($species_counter{$species} >=2){
			#print "$species has less than 2 samples\n";
#			goto NEXT_LINE;
#		}
#	}
	$sitecounter++;			
  	#Now to calculate the genetic distance
	my %already_done;
	foreach my $i (@good_number_list){
		foreach my $j (@good_number_list){
			unless(($geno_hash{$i}) and ($geno_hash{$j})){
				next;
			}
			if (($i == $j) or ($already_done{$j}{$i})){
				next;
			}
			$already_done{$i}{$j}++;
			#compare alleles
			my $distance = 0;
			foreach my $n (0..1){
				foreach my $m (0..1){
					if ($geno_hash{$i}{$n} ne $geno_hash{$j}{$m}){
						$distance += 0.25;
					}
				}
			}
			$data{$i}{$j} += $distance;
			$totalsites{$i}{$j}++;
		}
	}
}
print "sample1\tsample2\tspecies1\tspecies2\tspecies_comp\tdistance\tmarkers";
my %already_done;
foreach my $i (@good_number_list){
	foreach my $j (@good_number_list){
        	if (($i == $j) or ($already_done{$j}{$i})){
                	next;
                }
                $already_done{$i}{$j}++;
		print "\n$samplelist{$i}\t$samplelist{$j}\t$species{$i}\t$species{$j}";
		my $print1;
		my $species_comp;
		foreach my $species (@species_array){
			if ($species{$i} eq $species){
				if ($print1){
					$species_comp .= "-$species";
				}else{
					$species_comp = $species;
					$print1++;
				}
			}
			if ($species{$j} eq $species){
                                if ($print1){
                                        $species_comp .= "-$species";
                                }else{
                                        $species_comp = $species;
                                        $print1++;
                                }
                        }
		}
		print "\t$species_comp\t$data{$i}{$j}\t$totalsites{$i}{$j}";
	}
}
