#!bin/perl
use strict;
use warnings;
#This script summarizes geno-r2 from vcftools

my $current_pos;
my $site_counter = 0;
my $r2_sum = 0;
print "chr\tpos\tn\tr2";
while (<STDIN>){
    chomp;
    if ($. == 1){next;}
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    my $pos1 =$a[1];
    my $pos2 = $a[2];
    my $r2 = $a[4];
    if ($r2 eq "-nan"){next;}
    unless ($current_pos){
        $current_pos = "$chr\t$pos1";
    }
    if ($current_pos ne "$chr\t$pos1"){
        my $r2_final;
        if ($site_counter){
            $r2_final = $r2_sum / $site_counter;
        }else{
            $r2_final = "NA";
        }
        print "\n$current_pos\t$site_counter\t$r2_final";
        $r2_sum = 0;
        $site_counter = 0;
	$current_pos = "$chr\t$pos1";
    }
    
    $r2_sum+=$r2;
    $site_counter++;
}
