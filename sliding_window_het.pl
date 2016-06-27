#!/bin/perl
use warnings;
use strict;
#This sums up the counts of heterozygous and potentially heterozygous sites for the hybrid species. Feb 1, 2016. Greg Owens

my $current_chr;
my $window_size = $ARGV[0];
my @samples;
my $start = 0;
my $end = $start+$window_size;
my %species;
my %het;
my %sites;
my $samplefile = "Helianthus.all.forhetexp.specieslist.txt";

open SAMPLE, $samplefile;
while(<SAMPLE>){
    chomp;
    my @a = split(/\t/,$_);
    if(($a[1] ne "P1") and ($a[1] ne "P2")){
        push(@samples,$a[0]);
        $species{$a[0]} = $a[1];
    }
}
close SAMPLE;
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        print "chrom\tstart\tend\tsample\tspecies\tsites\taverage_het";
        next;}
    my $chr = $a[0];
    my $pos = $a[1];
    my $sample = $a[2];
    my $species = $a[3];
    my $H = $a[5];
    unless($current_chr){
        $current_chr = $chr;
    }
    if ($. == 2){
        until($pos < $end){
            $start += $window_size;
            $end = $start+$window_size;
        }
    }
    if (($current_chr ne $chr) or ($pos >= $end)){
		my $current_start = $start;
		my $current_end = $end -1;
        foreach my $name (@samples){
            unless($het{$name}){$het{$name} = 0;}
            unless($sites{$name}){next;}
            my $average_het = $het{$name}/$sites{$name};
            print "\n$current_chr\t$current_start\t$current_end\t$name\t$species{$name}\t$sites{$name}\t$average_het";
        }
        #Reset variables
        undef(%het);
        undef(%sites);
        if ($current_chr ne $chr){
            $current_chr = $chr;
            $start = 0;
            $end =  $start+$window_size;
        }
        until ($pos <$end){
            $start += $window_size;
            $end = $start+$window_size;
        }
    }
    $sites{$sample}++;
    if($H == 1){
        $het{$sample}++;
    }
}
