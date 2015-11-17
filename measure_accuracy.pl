#!/bin/perl
use warnings;
use strict;

my $simfile = $ARGV[0];

my %true_hash;
open SIM, $simfile;
while (<SIM>){
    chomp;
    if ($. == 1){next;}
    my @a = split(/\t/,$_);
    my $loc = $a[0];
    my $state = $a[13];
    $true_hash{$loc} = $state;
}
close SIM;
my $n_sites;
my $n_match;
while (<STDIN>){
    chomp;
    if ($. == 1){next;}
    my @a = split(/\t/,$_);
    my $chrom = $a[2];
    my $pos = $a[3];
    my $loc = "${chrom}_${pos}";
    my $state = $a[5];
    $n_sites++;
    if ($state eq $true_hash{$loc}){
        $n_match++;
#	print STDERR "$loc is matched\n";
    }else{
#	print STDERR "$loc is not matched. The measured state is $state, the true state i $true_hash{$loc}\n";
    }
}
my $percent_match = $n_match/$n_sites;
print "$percent_match";
