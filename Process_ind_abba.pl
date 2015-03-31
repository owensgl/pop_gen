#!/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];
my $length;
my %ABBA;
my %BABA;
my %total;
my %samplename;
open IN, $in;
while (<IN>){
	chomp;
	my @a = split(/\t/,$_);
	$length = $#a; 
	if ($. == 1){
		for my $i (2..$#a){
			$samplename{$i}= $a[$i];
		}
	}else{
		for my $i (2..$#a){
			if ($a[$i] eq "ABBA"){
				$ABBA{$i}++;
				$total{$i}++;
			}elsif ($a[$i] eq "BABA"){
				$BABA{$i}++;
				$total{$i}++;
			}
		}
	}
}
print "Sample1\tSample2\tSample3\tD\tABBA\tBABA\ttotal";
foreach my $i (2..$length){
	unless ($ABBA{$i}){
		$ABBA{$i} = 0;
	}
	unless ($BABA{$i}){
		$BABA{$i} = 0;
	}
	unless ($total{$i}){
		$total{$i} = 0;
	}
	my @names = split(/-/,$samplename{$i});
	print "\n$names[0]\t$names[1]\t$names[2]";
	my $D;
	if ($total{$i} > 0){
		$D = ($ABBA{$i} - $BABA{$i}) / $total{$i};
	}else{
		$D = "NA";
	}
	print "\t$D\t$ABBA{$i}\t$BABA{$i}\t$total{$i}";
}

