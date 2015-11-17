#!/bin/perl
use warnings;
use strict;

#Calculates genetic distance between each hybrid sample and it's parents. Also minimum distance. Takes a tab delimited snp table
my $genefile = $ARGV[0];
my $samplefile = $ARGV[1];

my $min_dp = 5;
my $min_qual = 20;

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
open SAMPLEFILE, $samplefile;
while(<SAMPLEFILE>){
	chomp;
  my @a = split(/\t/,$_);
	$printlist{$a[0]}++;
  $specieslist{$a[0]} = $a[1]; #should have P1, P2 and any other name for hybrids
}
close SAMPLEFILE;


my %species;
my %samplelist;
my $counter;
my $current_start_search = 0;
my $current_chrom = "NA";
my $current_gene;
my @parents = qw(P1 P2);
my %totaldistance;
my %totalsites;
my @good_number_list;
print "gene\tstart\tend\tsample\tspecies\tparent\tdistance\ttotal_sites";
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
	if ($current_chrom ne $chrom){
		$current_start_search = 0;
		$current_chrom = $chrom;
	}
  if (($counter % 100000)== 0){
			print STDERR "Processing $chrom $pos...\n";
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
      foreach my $i (@good_number_list){
	if (($specieslist{$samplelist{$i}} ne "P1") and ($specieslist{$samplelist{$i}} ne "P2")){
        foreach my $parent (@parents){
          my $ave_dist = "NA";
          unless($totaldistance{$i}{$parent}){
            $totaldistance{$i}{$parent} = 0;
          }
          if ($totalsites{$i}{$parent}){
            $ave_dist = $totaldistance{$i}{$parent}/$totalsites{$i}{$parent};
          }else{
            $totalsites{$i}{$parent} = "NA";
          }
          print "\n$current_gene\t$genestart{$current_gene}\t$geneend{$current_gene}\t$samplelist{$i}\t$specieslist{$samplelist{$i}}\t$parent\t$ave_dist\t$totalsites{$i}{$parent}";
	}
        }
      }
      undef($current_gene);
      undef(%totalsites);
      undef(%totaldistance);
      goto REPEAT;
    }
  }
  GOTGENE:
  foreach my $i (@good_number_list){
    if (($a[$i] ne "NN") and ($species{$i} ne "P1") and ($species{$i} ne "P2")){ #For each hybrid sample with data
      foreach my $j (@good_number_list){
        if (($a[$j] ne "NN") and (($species{$j} eq "P1") or ($species{$j} eq "P2"))){ #For each parent with data
          my $dist = measure_distance($a[$i], $a[$j]);
          $totalsites{$i}{$species{$j}}++;
          $totaldistance{$i}{$species{$j}}+=$dist; 
          }
        }
      }
    }
  }


sub measure_distance{
  my @baselist = qw(A T C G);
  my $sample1 = shift;
  my $sample2 = shift;
  my @samples1 = split(//,$sample1);
  my @samples2 = split(//,$sample2);
  my %alleles;
  $alleles{$samples1[0]}++;
  $alleles{$samples1[1]}++;
  $alleles{$samples2[0]}++;
  $alleles{$samples2[1]}++;
  my $num_alleles = scalar keys %alleles;
  my $shared_alleles = 0;
  foreach my $base (@baselist){
    if (($sample1 =~ m/$base/) and ($sample2 =~ m/$base/)){
      $shared_alleles++;
    }
  }
  my $distance = 1- $shared_alleles/$num_alleles;
  return($distance);
}
