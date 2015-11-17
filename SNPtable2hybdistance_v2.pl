#!/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);
#Calculates genetic distance between all species, including standard error using genpofad distance. Must have two individuals per species for each site used.
my $genefile = $ARGV[0];
my $samplefile = $ARGV[1];


my %genehash;
my %genestart;
my %geneend;
my %genechr;
my @genelist;
my %gene_array_name;
my %gene_array_start;
my %gene_array_end;
open GENEFILE, $genefile;
while (<GENEFILE>){
        chomp;
        if ($_ =~/^\>/g){
                my @a = split(/\ /, $_);
                my $begin = $a[2];
                $begin =~ s/begin=//g;
                my $end = $a[3];
                $end =~ s/end=//g;
                my $name = $a[0];
                $name =~ s/>//g;
                my $chr = $a[5];
                $chr =~ s/chr=//g;
                $genestart{$name} = $begin;
                $geneend{$name} = $end;
                $genechr{$name} = $chr;
                push (@{$gene_array_name{$chr}}, $name);
		push(@{$gene_array_start{$chr}}, $begin);
		push(@{$gene_array_end{$chr}}, $end);
        }
}
close GENEFILE;

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
my $current_start_search = 0;
my $current_chrom = "NA";
my $current_gene;
my %totaldistance;
my %totalsites;
my @good_number_list;
my %species_array;
my %distances_hash_array;
#print "gene\tstart\tend\tsample\tspecies\tparent\tdistance\ttotal_sites";
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
		print "gene\tstart\tend\tpair\tdistance\tse\tn";
		next;
	}
	unless ($a[2]){
		next;
	}
  	my $pos = $a[1];
	my $chrom = $a[0];
	if ($current_chrom ne $chrom){
		$current_start_search = 0;
		$current_chrom = $chrom;
	}
	if (($counter % 100000)== 0){
		print STDERR "Processing $chrom $pos...\n";
	}
	my %species_counter;
	my %total_alleles;
	my %geno_hash;
	foreach my $i (2..$#a){
		if (($a[$i] ne "NN") and ($species{$i})){
			$species_counter{$species{$i}}++;
			my @bases = split(//,$a[$i]);
			$total_alleles{$bases[0]}++;
			$total_alleles{$bases[1]}++;
			$geno_hash{$species{$i}}{$bases[0]}++;
			$geno_hash{$species{$i}}{$bases[1]}++;
			$geno_hash{$i}{$bases[0]}++;
			$geno_hash{$i}{$bases[1]}++;
		}
	}
	if (keys %total_alleles == 1){
		#next; #Ignore invariant sites;
	}
	foreach my $species (@species_array){ #Make sure there are atleast two samples genotyped per species
		unless ($species_counter{$species}){
			#print "$species has no data\n";
			goto NEXT_LINE;
		}
		unless($species_counter{$species} >=2){
			#print "$species has less than 2 samples\n";
			goto NEXT_LINE;
		}
	}
				
	REPEAT:
	my $num = @{$gene_array_start{$chrom}} -1;
	unless($current_gene){
		foreach my $j ($current_start_search..$num){
#				print "For $pos I'm trying $gene_array_name{$chrom}[$j] that starts at $gene_array_start{$chrom}[$j]\n";
			if (($gene_array_start{$chrom}[$j] <= $pos) and ($gene_array_end{$chrom}[$j] >= $pos)){
				$current_gene = $gene_array_name{$chrom}[$j];
				$current_start_search = $j;
				goto GOTGENE;
			}elsif ($gene_array_start{$chrom}[$j] > $pos){
				$current_start_search = $j;
#				print "START SEARCH AT $current_start_search\n";
				goto SKIP;
			}
		}
		SKIP:
		next;
  	}else{
    	#it's not in the same gene or in the next chromosome
    		if (($geneend{$current_gene} < $pos) or ($genechr{$current_gene} ne $chrom)){
#				print "Processing $pos which is after $current_gene ($genestart{$current_gene} - $geneend{$current_gene})\n";
			LASTLINE:
        		foreach my $i (0..$#species_array){
                		foreach my $j ($i..$#species_array){
                        		if ($i eq $j){next;}
					my $pair = "$species_array[$i]-$species_array[$j]";
          				my $ave_dist = "NA";
       	   				unless($totaldistance{$pair}){
				        	$totaldistance{$pair} = 0;
          				}
					my $standard_error;
          				if ($totalsites{$pair}){
		       				$ave_dist = $totaldistance{$pair}/$totalsites{$pair};
						my @dist_array = @{$distances_hash_array{$pair}};
						my $ss = 0; #calculate standard error
						foreach my $dist (@dist_array){
							my $value = ($dist - $ave_dist)**2;
							$ss += $value;
						}
						my $divisor = $totalsites{$pair}*($totalsites{$pair} - 1);
						if ($divisor eq "0"){
							$standard_error = "NA";
						}else{
							$standard_error = $ss/ $divisor;
						}
          				}else{
						$totaldistance{$pair} = "NA";
            					$totalsites{$pair} = "NA";
						$standard_error = "NA";
          				}
          				print "\n$current_gene\t$genestart{$current_gene}\t$geneend{$current_gene}\t$pair\t$ave_dist\t$standard_error\t$totalsites{$pair}";
        			}
      			}
      			undef($current_gene);
      			undef(%totalsites);
      			undef(%totaldistance);
			undef(%distances_hash_array);
      			goto REPEAT;
    		}
  	}
  	GOTGENE:
	foreach my $i (0..$#species_array){
        	foreach my $j ($i..$#species_array){
                	if ($i eq $j){next;}
			$totalsites{"$species_array[$i]-$species_array[$j]"}++;
                        my @alleles = keys %{$geno_hash{$species_array[$i]}};
                        my $total_dist;
                        foreach my $n (0..$#alleles){
                        	unless ($geno_hash{$species_array[$j]}{$alleles[$n]}){
                                	$geno_hash{$species_array[$j]}{$alleles[$n]} = 0;
                                }
                                my $dist = ($geno_hash{$species_array[$i]}{$alleles[$n]}/($species_counter{$species_array[$i]}*2)) * (1 - $geno_hash{$species_array[$j]}{$alleles[$n]}/($species_counter{$species_array[$j]}*2)); #p1 * (1-p2)
                                $total_dist += $dist;
                                $totaldistance{"$species_array[$i]-$species_array[$j]"} +=$total_dist;
                                push(@{$distances_hash_array{"$species_array[$i]-$species_array[$j]"}},$total_dist);
			}
     		}
    	}
	if (eof()){
		goto LASTLINE;
	}
	NEXT_LINE:
}

