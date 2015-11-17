#!/bin/perl
use warnings;
use strict;
#This merges probabilities for a the forward-backward viterbi algorithm. It assumes the forward and backward result files have been pasted together
while (<STDIN>){
	chomp;
	if ($. == 1){
		print "sample\tspecies\tchrom\tbp\tcm\tstate\tlikelihoodP1\tlikelihoodP2";
		next;
	}
	my @a = split(/\t/,$_);
	my $P1_for = $a[5];
	my $P2_for = $a[6];
	my $P1_rev = $a[12];
	my $P2_rev = $a[13];
	my $P1_both = $P1_for * $P1_rev;
	my $P2_both = $P2_for * $P2_rev;
	my $norm_factor = 1 / ($P1_both + $P2_both);
	$P1_both = $P1_both * $norm_factor;
	$P2_both = $P2_both * $norm_factor;
	my $state;
	if ($P1_both > $P2_both){
		$state = "P1";
	}else{
		$state = "P2";
	}
	print "\n$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$state\t$P1_both\t$P2_both";
}
