#!/bin/perl
use warnings;
use strict;
#Usage: cat snptable | perl thiscript.pl parentfile > parentcounts.txt
my $parentfile = $ARGV[0]; #This is a just a list of parents. One per line

my %parents;
my @parentlist;
open PARENTS, $parentfile;
while (<PARENTS>){
    chomp;
    $parents{$_}++;
    push(@parentlist,"$_.1");
    push(@parentlist,"$_.2");
}
close PARENTS;

my %sample;
#The STDIN is a tab separated phased file. genotypes should be A|T, or the like.
while(<STDIN>){
    chomp;
    my @a = split(/\t/,$_);
    if ($. == 1){
        foreach my $i (2..$#a){
            $sample{$i} = $a[$i];
        }
        print "$_";
    }else{
        my %data;
        my %count_alelles;
        #Load up information on site for parent samples
        foreach my $i (2..$#a){
            if ($parents{$sample{$i}}){
                if ($a[$i] eq "N|N"){
                    goto NEXTLINE;
                }
                my @alelles = split(/\|/,$a[$i]);
                $count_alelles{$alelles[0]}++;
                $count_alelles{$alelles[1]}++;
                $data{"$sample{$i}.1"} = $alelles[0];
                $data{"$sample{$i}.2"} = $alelles[1];
            }
        }
        if (keys %count_alelles == 0){next;} #Skip sites where the parents have no variation.
        print "\n$a[0]\t$a[1]";
        foreach my $i (2..$#a){
            if ($a[$i] eq "N|N"){
                print "\tNA";
            }else{
                my @alelles = split(/\|/,$a[$i]);
                print "\t";
                foreach my $n (0..1){
                    my $firstprint;
                    foreach my $parent(@parentlist){
                        if ($alelles[$n] eq $data{$parent}){
                            unless($firstprint){
                                print "$parent";
                                $firstprint++;
                            }else{
                                print ",$parent";
                            }
                        }
                    }
                    unless($firstprint){
                        print "*"; #This is when it doesn't match any parental allele
                    }
                    if ($n eq 0){
                        print "|";
                    }
                }

            }
        }
    }
    NEXTLINE:
}
