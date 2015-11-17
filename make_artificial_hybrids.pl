#!/bin/perl
use warnings;
use strict;

my $samplefile = $ARGV[0];
my $junctionfile = $ARGV[1];
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

my %samplelist;
my %species;
my %good_number_hash;
my @good_number_list;
my %P1_number_hash;
my %P2_number_hash;
my %junction_hash;


open JUNCTIONS, $junctionfile;
while(<JUNCTIONS>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){next;}
    $junction_hash{$a[0]}{$a[1]}{$a[2]}++; #Sample, chr, location

}
my $n_samples = 1; #Number of simulated genotypes
my %state_hash;
my $counter = 0;
my $current_chrom;
while(<STDIN>){
        $counter++;
        chomp;
        my $line = "$_";
        my @a = split(/\t/,$line);
        if ($. == 1){
            print "$a[0]";
            foreach my $i (1..11){
                print "\t$a[$i]";
            }
            foreach my $i (12..$#a){
                $samplelist{$i} = $a[$i];
                if ($specieslist{$a[$i]}){
                    $species{$i} = $specieslist{$a[$i]};
                    push(@good_number_list, $i);
                    $good_number_hash{$i}++;
                    if ($species{$i} eq "P1"){
		#	print "\t$a[$i]";
                        $P1_number_hash{$i}++;
                    }elsif ($species{$i} eq "P2"){
		#	print "\t$a[$i]";
                        $P2_number_hash{$i}++;
                    }
                }
            }
            foreach my $sample (1..$n_samples){
                print "\tsample_$sample\tstate_$sample";
            }
            next;
        }
        my $pos = $a[3];
        my $chrom = $a[2];
        my $cm = $a[4];
	if ($cm eq "NA"){next;}
        unless($chrom =~ m/^Ha/){next;}
        if ($chrom =~ m/^Ha0/){next;}
        unless ($current_chrom){
            $current_chrom = $chrom;
        }
        if ($current_chrom ne $chrom){
            undef(%state_hash);
            $current_chrom = $chrom;
        }
        my %P1alleles;
        my %P2alleles;
        my $P1count = 0;
        my $P2count = 0;
        if (($counter % 100000)== 0){
          print STDERR "Processing $chrom $pos...\n";
        }
        foreach my $i (keys %good_number_hash){ #Load up parental alleles
            if ($a[$i] ne "NN"){
                if ($species{$i} eq "P1"){
                    $P1count++;
                }elsif($species{$i} eq "P2"){
                    $P2count++;
                }
            }
        }
        unless(($P1count >=5) and ($P2count >= 5)){next;}
        #Determine the state for each replicate
        foreach my $sample (1..$n_samples){
            unless ($state_hash{$sample}){ #Pick an initial state
                my $rand = rand(2);
                if ($rand > 1){
                    $state_hash{$sample} = "P1";
                }else{
                    $state_hash{$sample} = "P2";
                }
            }
            my @current_junctions = keys %{$junction_hash{$sample}{$chrom}};
            foreach my $junction (@current_junctions){
                if ($junction < $cm){ #We crosses a junction
                    if ($state_hash{$sample} eq "P1"){ #switch your state
                        $state_hash{$sample} = "P2";
                    }else{ $state_hash{$sample} = "P1";}
                    delete($junction_hash{$sample}{$chrom}{$junction}) #Don't cross same junction
                }
            }
        }

        print "\n$a[0]";
        foreach my $i (1..11){
            print "\t$a[$i]";
        }
        foreach my $sample (1..$n_samples){
            if ($state_hash{$sample} eq "P1"){
                foreach my $i (keys %P1_number_hash){
                    if ($a[$i] ne "NN"){
                        print "\t$a[$i]";
                        goto MOVEON;
                    }
                }
            }else{
                foreach my $i (keys %P2_number_hash){
                    if ($a[$i] ne "NN"){
                        print "\t$a[$i]";
                        goto MOVEON;
                    }
                }
            }
            MOVEON:
            print "\t$state_hash{$sample}";
        }
}
