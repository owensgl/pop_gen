#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin/pop_gen/'; #For GObox server
my %t;
$t{"N"} = "NN";
$t{"A"} = "AA";
$t{"T"} = "TT";
$t{"G"} = "GG";
$t{"C"} = "CC";
$t{"W"} = "TA";
$t{"R"} = "AG";
$t{"M"} = "AC";
$t{"S"} = "CG";
$t{"K"} = "TG";
$t{"Y"} = "CT";

my $in = $ARGV[0]; #SNP table
my $pop = $ARGV[1]; #Pop file with P1, P2, and H
require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;

my %pop;
my %pophash;
my %samplepop;
my %unique_alleles;

open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);
	$pop{$a[0]}=$a[1];
	$pophash{$a[1]}++;
}
close POP;
my @poplist = keys %pophash;

open IN, $in;
while (<IN>){
	chomp;
	my @a = split (/\t/, $_);
	if ($. == 1){
		foreach my $i($badcolumns..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
			}
		}
	}else{
		next if /^\s*$/;
		my %BC;
		my %total_alleles;
		foreach my $i ($badcolumns..$#a){
			if ($samplepop{$i}){
				if ($iupac_coding eq "TRUE"){
						$a[$i] = $t{$a[$i]};
				}
				unless (($a[$i] eq "NN")or($a[$i] eq "XX")){
					my @bases = split(//, $a[$i]);
					$total_alleles{$bases[0]}++;
					$total_alleles{$bases[1]}++;
					$BC{$samplepop{$i}}{$bases[0]}++;
					$BC{$samplepop{$i}}{$bases[1]}++;
					$BC{$samplepop{$i}}{"Calls"}++;
					$BC{$samplepop{$i}}{"Calls"}++;
				}
			}
		}
		if (keys %total_alleles ==2){
			my @bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
			my $b1 = $bases[1]; #Major
			my $b2 = $bases[0];	#Minor
			foreach my $popname (@poplist){
				if ($BC{$popname}{"Calls"}){
					if ($BC{$popname}{"Calls"}< 10){
						#print "Less than 5 calls\n";
						goto SKIP;
					}
				}else{
					#print "No data for a pop\n";
					goto SKIP;
				}
			}
			my $count = 0;
			my $tmpholder;
			foreach my $popname (@poplist){
				if ($BC{$popname}{$b2}){
					$count++;
					$tmpholder = $popname;
				}
			}
			if ($count == 1){
				$unique_alleles{$tmpholder}++;
			}		

		}
	}
	SKIP:
}

foreach my $popname (@poplist){
	print "$popname\t$unique_alleles{$popname}\n";
}
