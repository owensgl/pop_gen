#This script takes a list of fasta files and calculates the genetic distance between each sample.
#Requires all bases to be in one line (not multiline fasta)
#Fasta sequences must all be the same length and not have multi-nucleotide IUPAC bases.
#!/bin/perl
use strict;
use warnings;

my %bases;
my %names;
my $counter = 1;
foreach my $i (0..$#ARGV){
    my $file = $ARGV[$i];
    open FILE, $file;
    while(<FILE>){
        chomp;
        if (/^>/){
            $names{$counter} = $_;
        }else{
            $bases{$counter} = $_;
	    $counter++;
        }
    }
}
foreach my $i (1..($counter-1)){
    foreach my $j (($i+1)..($counter-1)){
        print "$names{$i}|$names{$j}|";
        my @nuc1 = split(//,$bases{$i});
        my @nuc2 = split(//,$bases{$j});
        my $length= $#nuc1;
        my $dist = 0;
        foreach my $n (0..$length){
            if ($nuc1[$n] ne $nuc2[$n]){
                $dist++;
            }
        }
        my $rel_dist = $dist / $length;
        print "$rel_dist\n";
    }
}
