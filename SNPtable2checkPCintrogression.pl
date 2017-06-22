#!/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);
#This script takes a list of samples with species ID, a list of samples with PC scores, and a list of PC loadings with SNP ID
#It also gets piped in a SNP table that also has those SNPs.
#It divides the samples in the PCA by their score in a PC (+ or -), then for each snp in the PCA, it calculates the allele frequency for + group, - group and then asks if the allele freq in another species is more like one of the groups. It also prints out the PC loading and relative ranking (used for plotting)

my $speciesfile = $ARGV[0];
my $PCsamples_file = $ARGV[1];
my $PCloading_file = $ARGV[2];

#Load the species ID
my %specieslist;
my @species_array;
open POP, $speciesfile;
while (<POP>){
    chomp;
    my @a = split(/\t/, $_);
    $specieslist{$a[0]} = $a[1];
    push(@species_array, $a[1]);
}
my @tmp = uniq @species_array;
@species_array = @tmp;
close POP;

#Load the PCA scores for samples
my %samplePC;
my $number_of_PC; #Not counting zero.
open PCSAMPLE, $PCsamples_file;
while(<PCSAMPLE>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        $number_of_PC = $#a;
        next;
    }
    my $sample = $a[0];
    foreach my $i (1..$#a){
        $samplePC{$i}{$sample} = $a[$i];
    }
}
close PCSAMPLE;

#Load the PC loadings
my %loadingPC;
my %PC_snps;
open PCLOADING, $PCloading_file;
while(<PCLOADING>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        next;
    }
    my $snp_ID = $a[0];
    $PC_snps{$snp_ID}++;
    foreach my $i (1..$number_of_PC){
        $loadingPC{$i}{$snp_ID} = abs($a[$i]);
    }
}
close PCLOADING;

#Rank the PC loadings per snp
my %PC_rankings;

foreach my $i (1..$number_of_PC){
    #my @sorted_snps = sort { %{$loadingPC{$i}{$b}} <=> %{$loadingPC{$i}{$a}} } keys %{$loadingPC{$i}};
    my @sorted_snps = sort { $loadingPC{$i}{$b} <=> $loadingPC{$i}{$a} } keys %{$loadingPC{$i}}; 
    foreach my $j (0..$#sorted_snps){
        $PC_rankings{$i}{$sorted_snps[$j]} = $j;
#	print STDERR "For PC $i, $sorted_snps[$j] = score $loadingPC{$i}{$sorted_snps[$j]}and rank $j\n";
    }
}

my %ID_hash;
#Print header line
print "snp_id\tpc\tabs_pc_loading\tpc_ranking\tspecies\tpos_freq\tneg_freq\tspecies_freq\tpos_dif\tneg_dif\tresult";

#Load in SNP table.
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        foreach my $i (2..$#a){
            $ID_hash{$i} = $a[$i];
        }
    }else{
        my $chr = $a[0];
        my $pos = $a[1];
        my $snp_ID = "${chr}_${pos}";
        unless ($PC_snps{$snp_ID}){ #Skip the site if it isn't used in the PCA
            next;
        }
        #For each PC
        foreach my $pc (1..$number_of_PC){
            my %total_alleles; #For both positive and negative with alleles
            my %group_alleles; #positive and negative separated with alleles
            my $total_count = 0; #counts for all individuals
            my %group_count; #Counts for positive and negatives
            #load in the PCA allele frequencies
            foreach my $i (2..$#a){
                if ($samplePC{$pc}{$ID_hash{$i}}){ #If the sample is in the PCA calculations
                    if ($samplePC{$pc}{$ID_hash{$i}} > 0){ #For positive samples
                        if ($a[$i] eq "NN"){next;} #skip if missing data
                        my @bases = split(//,$a[$i]);
                        $group_alleles{'+'}{$bases[0]}++;
                        $group_alleles{'+'}{$bases[1]}++;
                        $group_count{'+'}+=2;
                        $total_alleles{$bases[0]}++;
                        $total_alleles{$bases[1]}++;
                        $total_count+=2;
                    }elsif ($samplePC{$pc}{$ID_hash{$i}} <= 0){ #For negative samples
                        if ($a[$i] eq "NN"){next;} #skip if missing data
                        my @bases = split(//,$a[$i]);
                        $group_alleles{'-'}{$bases[0]}++;
                        $group_alleles{'-'}{$bases[1]}++;
                        $group_count{'-'}+=2;
                        $total_alleles{$bases[0]}++;
                        $total_alleles{$bases[1]}++;
                        $total_count+=2;
                    }
                }
            }
	    #Check for no data;
	    unless ($group_count{'-'} and $group_count{'+'}){next;}
            #Get allele frequencies for the PCA groups
            my @total_bases = sort { $total_alleles{$a} <=> $total_alleles{$b} } keys %total_alleles ;
            my $major = $total_bases[1];
            my $minor = $total_bases[0];
            my %allele_freq;
            #For the positive group
            if ($group_alleles{'+'}{$major}){
                $allele_freq{'+'} = $group_alleles{'+'}{$major} / $group_count{'+'};
            }else{
                $allele_freq{'+'} = 0;
            }
            #For the negative group
            if ($group_alleles{'-'}{$major}){
                $allele_freq{'-'} = $group_alleles{'-'}{$major} / $group_count{'-'};
            }else{
                $allele_freq{'-'} = 0;
            }
            #Load the data for each species one at a time.
            foreach my $species (@species_array){
                my %species_alleles;
                my $species_count = 0;
                foreach my $i (2..$#a){
                    if ($a[$i] eq "NN"){next;}
                    if ($specieslist{$ID_hash{$i}}){
                        if ($specieslist{$ID_hash{$i}} eq $species){
                            my @bases = split(//,$a[$i]);
                            $species_alleles{$bases[0]}++;
                            $species_alleles{$bases[1]}++;
                            $species_count+=2;
                            #Check if it is still bialleleic
                            foreach my $j (0..1){
                                if (($bases[$j] ne $major) and ($bases[$j] ne $minor)){
                                    goto NEXTSPECIES;
                                }
                            }
                        }
                    }
                }
                if ($species_count == 0){goto NEXTSPECIES;} #If no data move onto next species
                #Get species allele frequency
                my $species_freq;
                if ($species_alleles{$major}){
                    $species_freq = $species_alleles{$major} / $species_count;
                }else{
                    $species_freq = 0;
                }
                #Check if it is closer to positive or negative group
                my $pos_dif = abs($allele_freq{'+'} - $species_freq);
                my $neg_dif = abs($allele_freq{'-'} - $species_freq);
                my $result_closest;
                if ($pos_dif < $neg_dif){
                    $result_closest = "1";
                }elsif($pos_dif > $neg_dif){
                    $result_closest = "-1";
                }else{
                    $result_closest = "NA";
                }
		
		my $tmp_rank = $PC_rankings{$pc}{$snp_ID}+1;
                print "\n$snp_ID\t$pc\t$loadingPC{$pc}{$snp_ID}\t$tmp_rank\t";
                print "$species\t$allele_freq{'+'}\t$allele_freq{'-'}\t$species_freq\t$pos_dif\t$neg_dif\t";
                print "$result_closest";
                NEXTSPECIES:
            }
        }
    }
}
