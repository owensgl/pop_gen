#!/usr/bin/perl

use warnings;
use strict;
use lib './'; #For GObox server
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
my %samplepop;
my %samplename;
my $endcolumn;

my $table1 = $ARGV[0]; #table with snps you want to keep
my $table2 = $ARGV[1]; #larger snp table to be subset


require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($table2);
$. = 0;




open TABLE2, $table2;
while (<TABLE2>){
	chomp;
	my @a = split (/\t/,$_);
	$endcolumn = $#a;
	if ($. == 1){
		foreach my $i ($badcolumns..$#a){ #Get sample names for each column
        		$samplename{$i} = $a[$i];
        	}
	}else{
		my $loc = "$a[2]\t$a[3]";
		foreach my $i ($badcolumns..$#a){ 
			if ($iupac_coding eq "TRUE"){
				$a[$i] = $t{$a[$i]};
			}
			$h{$loc}{$samplename{$i}} = $a[$i];
		}
	}
}



open TABLE1, $table1;
while (<TABLE1>){
	chomp;
	my @a = split (/\t/,$_);
	if ($. == 1){
		print "$_";
		foreach my $i ($badcolumns..$endcolumn){
			print "\t$samplename{$i}";
		}
		print "\n";
       	}
       	else{
       		my $loc = "$a[0]\t$a[1]";
		print "$_";
		foreach my $i ($badcolumns..$endcolumn){
			if ($h{$loc}{$samplename{$i}}){
				print "\t$h{$loc}{$samplename{$i}}";
			}else{
				print "\tNN";
			}
		}
		print "\n";
	}
}

close TABLE1;
			
        	

