#!/usr/perl
#This script takes in a hapmap file and a population label file, and calculates the genetic distance between each pair of samples for a sliding window.
#Eventually it will call an R script to run a PCoA and extract euclidean distances between pairs.
#Make sure that a sample is not intermediate because it doesn't have any data. 
use warnings;
use strict;
use File::Basename;
use Cwd 'abs_path';
use List::Util qw(shuffle);
use lib '/home/owens/bin/pop_gen/'; #For GObox server

my $path = abs_path($0); #Get bin path for R script calling
$path =~ s/\/SNPtable2distances_SW.pl//g;

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

my %flip; #Flip bases so they're always in alphabetical order
$flip{"TA"} = "AT";
$flip{"GA"} = "AG";
$flip{"CA"} = "AC";
$flip{"GC"} = "CG";
$flip{"TC"} = "CT";
$flip{"TG"} = "GT";



my $min_sites = 100; #Minumum number of sites with data for a pair within a window. 
my $windowsize = 100000; #Number of bases in each sliding window
my $max_null_sites = ($windowsize - $min_sites); #Max number of null sites within a sliding window for a pair
my $max_null_pairs; #max number of excluded pairs for a sliding window. Calculated as half of possible pairs.

my $infile = $ARGV[0]; #A snp table in hmp format. Only include samples you're using.
my $popfile = $ARGV[1]; #A population file.
my $out = $ARGV[2]; #outfile prefix
my $outname = "$out.txt";
open (OUTFILE, "> $outname") or die "Could not open a file\n";

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

require "countbadcolumns.pl"; #Identify the number of columns before genotype data starts
my ($iupac_coding, $badcolumns) = count_bad_columns($infile);
$. = 0;

#Distance hashes
my %nullcounter;
my %sitecount;
my %distance;
my $finalcolumn;
my $pop1;
my $pop2;
my $parenttotal;
my %parenthash;

open IN, $infile;
while (<IN>){
	chomp;
	my @a = split(/\t/,$_);
	
	if ($. == 1){ #Load in sample names associated with column numbers, as well as population
		$finalcolumn = $#a;
		my $number_samples = (($#a - $badcolumns) + 1);
		$max_null_pairs = ($number_samples * ($number_samples -1)) / 2;
		foreach my $i($badcolumns..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				$samplename{$i} = $a[$i];
				if ($samplepop{$i} eq "P1"){
					$pop1++;
					$parenthash{$i} = $samplepop{$i};
				}elsif ($samplepop{$i} eq "P2"){
					$pop2++;
					$parenthash{$i} = $samplepop{$i};
				}
			}
		}
		$parenttotal = $pop1 + $pop2;
	}else{
		my $loc = $a[0];
		foreach my $i($badcolumns..$#a){ #If IUPAC coding used, convert to paired nucleotide.
			if ($iupac_coding eq "TRUE"){
				$a[$i] = $t{$a[$i]};
			}
			else{
				if ($flip{$a[$i]}){
					$a[$i] = $flip{$a[$i]};
				}
			}
		}
		foreach my $i($badcolumns..$#a){
			foreach my $j ($i..$#a){
				if (($a[$i] eq "NN") or ($a[$j] eq "NN")){ #If either is missing, count as missing for that pair.
					$nullcounter{$i}{$j}++;
				}else{
					my @tmp1 = split('',$a[$i]);
					my @tmp2 = split('',$a[$j]); 
					if ($a[$i] eq $a[$j]){ #If they're identical, don't add to the distance but do add to the site count.
						$distance{$i}{$j} += 0;
						$sitecount{$i}{$j}++;
					}elsif (($tmp1[0] eq $tmp2[0]) or ($tmp1[0] eq $tmp2[1]) or ($tmp1[1] eq $tmp2[0]) or ($tmp1[1] eq $tmp2[1])){ #If any have any matches, they are 0.5 distance away.
						$sitecount{$i}{$j}++;
						$distance{$i}{$j} += 0.5;
					}else{
						$sitecount{$i}{$j}++;
						$distance{$i}{$j}++;
					}
				}
			}
		}
	}
}
close IN;
my $paircount = 0;
my $nullpairs = 0;
my $totaldist = 0;
foreach my $i ($badcolumns..$finalcolumn){
	foreach my $j ($i..$finalcolumn){
		if ($nullcounter{$i}{$j} > $max_null_sites){
			$distance{$i}{$j} = "NA"; #Delete pairs without sufficient data within the window.
			$nullpairs++;
		}
		if ($distance{$i}{$j} ne "NA"){
			$distance{$i}{$j} = ($distance{$i}{$j} / $sitecount{$i}{$j});
			if ($i ne $j){
				$paircount++;
				$totaldist += $distance{$i}{$j} #Add up all distances, not including self distances
			}
		}
	}
}
#calculate average distance for window
my $average_dist = ($totaldist / $paircount);
if ($nullpairs < $max_null_pairs){
	foreach my $i ($badcolumns..$finalcolumn){
		foreach my $j ($badcolumns..$finalcolumn){
			unless ($distance{$i}{$j}){
                                $distance{$i}{$j} = $distance{$j}{$i};
			}if ($distance{$i}{$j} eq "NA"){
				$distance{$i}{$j} = $average_dist;
			}if ($j eq $badcolumns){
				print OUTFILE "$distance{$i}{$j}";
			}else{
				print OUTFILE "\t$distance{$i}{$j}";
			}
		}
		print OUTFILE "\n";
	}
}else{
	print OUTFILE "Too much missing data\n";
	exit;
}

#Use cmdscale in R to compress to 2 dimensions, then calculate the euclidean distance between each sample.
my $cmd = "Rscript $path/MDS_pipe.R $outname $out.dist.txt";
system($cmd);
my $i = 0;


my %eucdist;
my $interpopdist;
my $interpopcount;
my $pop1dist;
my $pop2dist;
my $pop1count;
my $pop2count;

open DIST, "$out.dist.txt"; #Load in distance scores from R
while (<DIST>){
	chomp;
	my @a = split(/\t/, $_);
	foreach my $j (0..$#a){
		$eucdist{($i+$badcolumns)}{($j+$badcolumns)} = $a[$j];
	}
	$i++;
}
close DIST;
foreach my $i ($badcolumns..($finalcolumn-1)){
	foreach my $j (($badcolumns+1)..$finalcolumn){
		if ((($samplepop{$i} eq "P1") and ($samplepop{$j} eq "P2")) or (($samplepop{$i} eq "P2") and ($samplepop{$j} eq "P1"))){ #If they're in different parental groups
			$interpopdist += $eucdist{$i}{$j};
			$interpopcount++;
		}elsif(($samplepop{$i} eq "P1") and ($samplepop{$j} eq "P1")){ #if they're within pop1
			$pop1dist += $eucdist{$i}{$j};
			$pop1count++;
		}elsif(($samplepop{$i} eq "P2") and ($samplepop{$j} eq "P2")){ #If they're within pop2
                        $pop2dist += $eucdist{$i}{$j};
                        $pop2count++;
		}
	}
}
my $average_within_dist = (($pop1dist / $pop1count) + ($pop2dist / $pop2count)) / 2;
my $CSS = ($interpopdist / $interpopcount) - $average_within_dist;
print "$CSS\tREAL\n";

foreach (1..100){
	my $loop_interpopdist;
	my $loop_interpopcount;
	my $loop_pop1dist;
	my $loop_pop2dist;
	my $loop_pop1count;
	my $loop_pop2count;
	my @values = shuffle(values %parenthash); #Shuffle hash of population assignment for the parents
	$parenthash{$_} = shift @values for keys %parenthash;
	foreach my $i ($badcolumns..($finalcolumn-1)){
	        foreach my $j (($badcolumns+1)..$finalcolumn){
			if (($parenthash{$i}) and ($parenthash{$j})){
				if ((($parenthash{$i} eq "P1") and ($parenthash{$j} eq "P2")) or (($parenthash{$i} eq "P2") and ($parenthash{$j} eq "P1"))){ #If they're in different parental groups
         	                	$loop_interpopdist += $eucdist{$i}{$j};
                	        	$loop_interpopcount++;
	              		}elsif(($parenthash{$i} eq "P1") and ($parenthash{$j} eq "P1")){ #if they're within pop1
	                        	$loop_pop1dist += $eucdist{$i}{$j};
	                        	$loop_pop1count++;
	                	}elsif(($parenthash{$i} eq "P2") and ($parenthash{$j} eq "P2")){ #If they're within pop2
	                        	$loop_pop2dist += $eucdist{$i}{$j};
	                        	$loop_pop2count++;
				}
			}
		}
	}
	my $loop_average_within_dist = (($loop_pop1dist / $loop_pop1count) + ($loop_pop2dist / $loop_pop2count)) / 2;
	my $loop_CSS = ($loop_interpopdist / $loop_interpopcount) - $loop_average_within_dist;
	print "$loop_CSS\n"; #THESE values are too high. They should go into the negative
}



