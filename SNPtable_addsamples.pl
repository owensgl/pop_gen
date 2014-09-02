#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/owens/bin'; #For GObox server
my %t;
$t{"N"} = "NN";
$t{"A"} = "AA";
$t{"T"} = "TT";
$t{"G"} = "GG";
$t{"C"} = "CC";
$t{"W"} = "TA";
$t{"R"} = "AG";
$t{"M"} = "AC";
$t{"S"} = "CG";
$t{"K"} = "TG";
$t{"Y"} = "CT";

my %h;


my $table1 = $ARGV[0]; #table with snps you want to keep
my $table2 = $ARGV[1]; #larger snp table to be subset


require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;



open TABLE1, $table1;
while (<TABLE1>){
	chomp;
	my @a = split (/\t/,$_);
		if ($. == 1){
        	}
        	else{
        		my $loc = "$a[0]\t$a[1]";
        		$h{$loc}++;
 		}
 	}
}
close TABLE1;

open TABLE2, $table2;
while (<TABLE2>){
	chomp;
		if ($. == 1){
			print "$_\n";
		}else{
			my $loc = "$a[2]\t$a[3]";
			if ($h{$loc}){
				foreach my $i ($badcolumns..$#a){ 
					if ($iupac_coding eq "TRUE"){
						$a[$i] = $t{$a[$i]};
					}
					print "\t$a[$i]";
				}print "\n";
			}else{
				print "\n";
			}
		}
	}
}
	
			
        	

