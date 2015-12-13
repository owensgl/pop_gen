#!/bin/perl
use warnings;
use strict;
use Math::CDF qw(:all);

my $sum;
foreach my $i (1..100){
	my $value = pbeta(($i)/100,6,6) - pbeta(($i-1)/100, 6, 6);
	$sum+= $value;	
	print "$value\n";
}
print "The total sum is $sum\n";
