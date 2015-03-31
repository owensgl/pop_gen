#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin'; #Location of countbadcolumns.pl
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

my $in = $ARGV[0]; #In SNP table in iupac or not. Will detect the number of columns before data
my $pop = $ARGV[1]; #Population file, (samplename\tpopname\n)
my $poporder = $ARGV[2]; #The file specify which populations represent which samples in the ABBA BABA scheme (popname\t[1-4]\n)

require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;

my %pop;
my %poplist;
my %group;
my %samplegroup;
my @names;

open POP, $pop;
while (<POP>){
        chomp;
        my @a = split (/\t/,$_);
        $pop{$a[0]}=$a[1];
        $poplist{$a[1]}++;
}
close POP;
#Load group information
open GROUP, $poporder;
while (<GROUP>){
        chomp;
        my @a = split (/\t/,$_);
        $group{$a[0]} = $a[1];
}
close GROUP;

my @group1;
my @group2;
my @group3;

open IN, $in;
while (<IN>){
	chomp;
	next if /^\s*$/;
	my @a = split(/\t/,$_);
	if ($. == 1){
		print "chrom\tpos";
		foreach my $i ($badcolumns..$#a){
			if ($pop{$a[$i]}){
				if ($group{$pop{$a[$i]}}){
					$samplegroup{$i} = $group{$pop{$a[$i]}};
					if ($samplegroup{$i} == 1){
						push (@group1, $i);
					}elsif ($samplegroup{$i} == 2){
						push (@group2, $i);
					}elsif ($samplegroup{$i} == 3){
						push (@group3, $i);
					}
				}
			}
		}
		foreach my $pop1sample (@group1){
			foreach my $pop2sample (@group2){
				foreach my $pop3sample (@group3){
					print "\t$a[$pop1sample]-$a[$pop2sample]-$a[$pop3sample]";
				}
			}
		}
				
	}else{
		my $pos = "$a[0]\t$a[1]";
		my %BC;
		my %total_alleles;
		foreach my $i ($badcolumns..$#a){
			if ($samplegroup{$i}){
				if ($iupac_coding eq "TRUE"){
                			$a[$i] = $t{$a[$i]};
                		}
                		unless ($a[$i] eq "NN"){
                			my @bases = split(//, $a[$i]);
                			$BC{$samplegroup{$i}}{$bases[0]}++;
                			$BC{$samplegroup{$i}}{$bases[1]}++;
					$total_alleles{$bases[0]}++;
					$total_alleles{$bases[1]}++;
                		}
			}
		}
		my $skip;
		for my $i (1..3){
			unless(keys %{$BC{$i}}){
				$skip++;	
			}
		}
		if (keys %total_alleles == 1){ $skip++};
		next if ($skip);
		if (keys %{$BC{4}}){
			if (keys %{$BC{4}} == 1){
				print "\n$pos";
				my @outgroupbases = sort { $BC{4}{$b} <=> $BC{4}{$a} } keys %{$BC{4}};
				my $ancestral = $outgroupbases[0];
				foreach my $pop1sample (@group1){
					my $missing1;
					if ($a[$pop1sample] eq "NN"){ $missing1++};
					foreach my $pop2sample (@group2){
						my $missing2;
						if ($a[$pop2sample] eq "NN"){ $missing2++};
						foreach my $pop3sample (@group3){
							my $missing3;
							my $het;
							my %group_alleles;
							if ($a[$pop3sample] eq "NN"){ $missing3++};
							my @allsamples = ($pop1sample, $pop2sample, $pop3sample);
							foreach my $i (@allsamples){
								my @bases = split(//, $a[$i]);
								$group_alleles{$bases[0]}++;
								$group_alleles{$bases[1]}++;
								unless ($bases[0] eq $bases[1]){
									$het++;
								}
							}
							if ((keys %group_alleles == 2) and ($total_alleles{$ancestral})){
								unless (($missing1) or ($missing2) or ($missing3) or ($het)){
									my $state1;
									my $state2;
									my $state3;
									if ($a[$pop1sample] eq "$ancestral$ancestral"){
										$state1 = "A";
									}else{
										$state1 = "B";
									}
									if ($a[$pop2sample] eq "$ancestral$ancestral"){
										$state2 = "A";
									}else{
										$state2 = "B";
									}
									if ($a[$pop3sample] eq "$ancestral$ancestral"){
										$state3 = "A";
									}else{
										$state3 = "B";
									}
									print "\t$state1$state2${state3}A";
								}else{
									print "\t-";
								}
							}else{ 
								print "\t-";
							}
						}
					}
				}
			}
		}
	}
}
