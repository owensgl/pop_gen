#!/bin/perl
use warnings;
use strict;

#This script takes a list of sites (in the initial case it was high Fst sites for ben-lim stickleback), and it asks how many of those sites does each sample have data for.

my $sites = $ARGV[0]; #List of sites in two columns (chrom\tpos)
#Snp table is piped in
open SITES, $sites;
my %sites_hash;
while (<SITES>){
    chomp;
    my @a = split(/\t/,$_);
    $sites_hash{$a[0]}{$a[1]}++;
}
close SITES;

my %samplelist;
my $number_ind;
my %data_counter;
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        foreach my $i (2..$#a){
            $samplelist{$i} = $a[$i];
        }
        $number_ind = $#a;
    }else{
        my $chrom = $a[0];
        my $pos = $a[1];
        if ($sites_hash{$chrom}{$pos}){
            foreach my $i (2..$#a){
                if ($a[$i] ne "NN"){
                    $data_counter{$i}++;
                }
            }
        }

    }
}
print "Sample\tNumber_sites";
foreach my $i (2..$number_ind){
    unless ($data_counter{$i}){
	$data_counter{$i} = 0;
    }
    print "\n$samplelist{$i}\t$data_counter{$i}";
}
