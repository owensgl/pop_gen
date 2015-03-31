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


my ( $name, $path, $suffix ) = fileparse( $0, qr{\.[^.]*$} );
#print "NAME=$name\n";
#print "PATH=$path\n";
#print "SFFX=$suffix\n";


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



my $min_sites = 2; #Minumum number of sites with data for a pair within a window. 
my $min_variable_sites = 5; #Number of sites in the dataset in a window. Skips the window if they don't have it. This is before considering missing data.
my $window_size = 1000000; #Number of bases in each sliding window
my $max_null_pairs; #max number of excluded pairs for a sliding window. Calculated as half of possible pairs.
my $permute_number = 200; #Number of permutations.
my $minimum_hybrid_comparisons = 5; #Minumum number of comparisons with each parental for hybrid calculation.
my $infile = $ARGV[0]; #A snp table in hmp format. Only include samples you're using.
my $popfile = $ARGV[1]; #A population file.
my $out = $ARGV[2]; #outfile prefix
open (FINALOUT, "> $out.final.txt") or die "Could not open a file\n"; #open outfile

my %pop;
my %popList;
my %samplepop;
my %samplename;

#Window variables
my $current_chrom;
my $end = $window_size;

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
my %took_average_dist;
my $finalcolumn;
my $pop1;
my $pop2;
my $popH;
my $parenttotal;
my %parenthash;
my %hybridhash;
my @hybridarray;
my $variablesites;

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
				}elsif ($samplepop{$i} eq "H"){
					$popH++;
					$hybridhash{$i} = $samplepop{$i};
					push (@hybridarray, $i);
				}
			}
		}
		$parenttotal = $pop1 + $pop2;
		print FINALOUT "chrom\tstart\tend\tCSS\talpha\tn_sites";
		foreach my $hybrid_number (@hybridarray){
			print FINALOUT "\t$samplename{$hybrid_number}";
		}
	}else{
		my $loc = $a[0];
		my $chrom = $a[2];
		my $pos = $a[3];
		unless ($current_chrom){
			$current_chrom = $chrom;
		}
		if ($. == 2){
			until($end > $pos){ #move the end point of the window until it is after the current marker.
                                $end += $window_size;
                        }
		}
		if (($current_chrom ne $chrom) or ($end < $pos)){ #If it's starting a new window, do all the calculations
			my $start = $end - $window_size;
			if ($variablesites < $min_variable_sites){
				print FINALOUT "\n$current_chrom\t$start\t$end\tNA\tNA\tNA";
				foreach my $hybrid_number (@hybridarray){
					print FINALOUT "\tNA";
				}
				goto RESET;
			}
			print "Calculating distances for $current_chrom, $start to $end\n";
			open (OUTFILE, "> $out.$current_chrom.$end.txt") or die "Could not open a file\n"; #open outfile
			my $paircount = 0;
			my $nullpairs = 0;
			my $totaldist = 0;
			foreach my $i ($badcolumns..($finalcolumn-1)){ #This is to create the average distance between objects for missing data
			        foreach my $j (($i+1)..$finalcolumn){
					unless ($sitecount{$i}{$j}){
						$sitecount{$i}{$j} = 0;
					}
			                if ($sitecount{$i}{$j} <= $min_sites){
                        			$distance{$i}{$j} = "NA"; #Delete pairs without sufficient data within the window.
			                        $nullpairs++;
						#print "Too many null on $i $j\n";
			                }
					else{
						unless($distance{$i}{$j}){
							$distance{$i}{$j} = 0;
						}
					}
			                if ($distance{$i}{$j} ne "NA"){
						#print "For $i $j nullcounts=$nullcounter{$i}{$j}\n";
						#print "For $i $j distance=$distance{$i}{$j}, sitcount=$sitecount{$i}{$j}\n";
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
			#print "$average_dist\t$nullpairs\t$max_null_pairs\n";
			if ($nullpairs < $max_null_pairs){
				foreach my $i ($badcolumns..$finalcolumn){
					foreach my $j ($badcolumns..$finalcolumn){
						if ($i eq $j){
							$distance{$i}{$j} = 0;
						}
						unless (defined $distance{$i}{$j}){
				                	$distance{$i}{$j} = $distance{$j}{$i};
						}if ($distance{$i}{$j} eq "NA"){
							$distance{$i}{$j} = $average_dist;
							$took_average_dist{$i}{$j}++;
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
				close OUTFILE;
				my $start = $end - $window_size;
                                print FINALOUT "\n$current_chrom\t$start\t$end\tNA\tNA\tNA";
				foreach my $hybrid_number (@hybridarray){
					print FINALOUT "\tNA";
				}
				goto RESET;
			}
			close OUTFILE;


			#Use cmdscale in R to compress to 2 dimensions, then calculate the euclidean distance between each sample.
			print "Calculating euclidean distances for $current_chrom, $start to $end\n";
			my $cmd = "Rscript $path/MDS_pipe.R $out.$current_chrom.$end.txt $out.$current_chrom.$end.dist.txt";
			system($cmd);
			my $i = 0;

			my %eucdist;
			my $interpopdist;
			my $interpopcount;
			my $pop1dist;
			my $pop2dist;
			my $pop1count;
			my $pop2count;
			my %hybrid_result;
			open DIST, "$out.$current_chrom.$end.dist.txt"; #Load in distance scores from R
			while (<DIST>){
				chomp;
				my @a = split(/\t/, $_);
				foreach my $j (0..$#a){
					$eucdist{($i+$badcolumns)}{($j+$badcolumns)} = $a[$j];
				}
				$i++;
			}
			close DIST;
			foreach my $hybrid_number (@hybridarray){
				open (HYBFILE, "> $out.$current_chrom.$end.Hyb$hybrid_number.txt") or die "Could not open a file\n"; #open file to put hybrid distances
				my $tmpP1count=0;
				my $tmpP2count=0;
				foreach my $i ($badcolumns..$finalcolumn){
					unless (($took_average_dist{$hybrid_number}{$i}) or ($took_average_dist{$i}{$hybrid_number})){
						if ($samplepop{$i} eq "P1"){
							$tmpP1count++;		
							print HYBFILE "$eucdist{$i}{$hybrid_number}\tP1\n";
						}elsif ($samplepop{$i} eq "P2"){
							$tmpP2count++;
							print HYBFILE "$eucdist{$i}{$hybrid_number}\tP2\n";
						}
					}
				}
				close HYBFILE;
				unless (($tmpP1count > $minimum_hybrid_comparisons) and ($tmpP2count > $minimum_hybrid_comparisons)){
					$hybrid_result{$hybrid_number} = "NA:NA:NA:NA";
					goto SKIPHYB;
				}
				my $cmd2 = "Rscript $path/MDS_pipe2.R $out.$current_chrom.$end.Hyb$hybrid_number.txt $out.$current_chrom.$end.Hyb$hybrid_number.ttest.txt";
				system($cmd2);
				open HYBRESULT, "$out.$current_chrom.$end.Hyb$hybrid_number.ttest.txt";
				while (<HYBRESULT>){
					my $tmp = $_;
					chop($tmp);
					$hybrid_result{$hybrid_number} = $tmp;
				}
				close HYBRESULT;
				SKIPHYB:
				my $remove_hybrid_data = "rm $out.$current_chrom.$end.Hyb$hybrid_number*";
				system($remove_hybrid_data);
				
			}
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
			print "Calculating CSS for $current_chrom, $start to $end. CSS = $CSS\n";
			my @permuted_CSS;
			print "Permuting CSS for $current_chrom, $start to $end\n";
			foreach (1..$permute_number){
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
				push (@permuted_CSS, $loop_CSS);
			}

			my $alpha = 0;
			foreach my $permute (sort {$b <=> $a} @permuted_CSS){
				if ($permute > $CSS){
					$alpha++;
				}elsif ($permute < $CSS){
					my $value = ($alpha / $permute_number);
					my $start = $end - $window_size;
					print FINALOUT "\n$current_chrom\t$start\t$end\t$CSS\t$value\t$variablesites";
					foreach my $hybrid_number (@hybridarray){
						print FINALOUT "\t$hybrid_result{$hybrid_number}";
					}
					goto RESET;
				}
			}
					
			#Reset the variables
			RESET:
			undef %distance;
			undef %sitecount;
			undef $variablesites;
			undef %took_average_dist;
			if ($current_chrom ne $chrom){ #If on the next chromosome, reset the end point for the window
				$current_chrom = $chrom;
				$end = $window_size;
			}
			until($end > $pos){ #move the end point of the window until it is after the current marker.
				$end += $window_size;
			}
		
		}
		$variablesites++;
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
				#nothing happens
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

