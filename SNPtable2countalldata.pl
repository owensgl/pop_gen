#!/bin/perl
use warnings;
use strict;

#This script counts the number of genotypes with data per sample
#Snp table is piped in

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
            foreach my $i (2..$#a){
                if ($a[$i] ne "NN"){
                    $data_counter{$i}++;
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
