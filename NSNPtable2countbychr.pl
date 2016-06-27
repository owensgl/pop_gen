#!/bin/perl
use warnings;
use strict;

#This script loads up a tab delimited file for benthic limnetics that is coded such that 0 is benthic and 1 is limnetic. It then goes chrom by chrom and adds up the total percent benthic or limnetic for each population

my $groupfile = $ARGV[0];

open GROUP, $groupfile;

my %poplist;
my %pophash;
while(<GROUP>){
    chomp;
    my @a = split(/\t/,$_);
    $poplist{$a[0]} = $a[1];
}
close GROUP;

my @samplelist;
my %samplehash;
my %chrlist;
my %sitescount;
my %limcount;
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    $chrlist{$chr}++;
    my $pos = $a[1];
    if ($. == 1){
        foreach my $i (2..$#a){
            if ($poplist{$a[$i]}){
                $pophash{$i} = $poplist{$a[$i]};
                $samplehash{$i} = $a[$i];
                push(@samplelist,$a[$i]);
            }
        }
    }else{
        foreach my $i (2..$#a){
		if ($samplehash{$i}){
            		if ($a[$i] ne "NA"){
                		if ($a[$i] eq "00"){
                    			$sitescount{$chr}{$samplehash{$i}}+=2;
                		}elsif ($a[$i] eq "01"){
                    			$sitescount{$chr}{$samplehash{$i}}+=2;
                    			$limcount{$chr}{$samplehash{$i}}++;
                		}elsif ($a[$i] eq "11"){
                    			$sitescount{$chr}{$samplehash{$i}}+=2;
                    			$limcount{$chr}{$samplehash{$i}}+=2;
				}
                	}
            	}
        }
    }
}
print "sample\tpopulation\tchr\tpercent_lim";
my @chrs = sort keys %chrlist;
foreach my $sample(@samplelist){
    foreach my $chr (@chrs){
        if ($sitescount{$chr}{$sample}){
            unless($limcount{$chr}{$sample}){
                $limcount{$chr}{$sample} = 0;
            }
            my $limpercent = $limcount{$chr}{$sample} / $sitescount{$chr}{$sample};
            print "\n$sample\t$poplist{$sample}\t$chr\t$limpercent";
        }
    }
}
