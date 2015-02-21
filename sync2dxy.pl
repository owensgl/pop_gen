#!/bin/perl

my $in = $ARGV[0];

open IN, $in;
while (<IN>){
	chomp;
	my @a = split(/\t/,$_);
	my $coast = $a[3];
	my $inland = $a[4];
	my @c_counts = split(/:/, $coast);
	my @i_counts = split(/:/, $inland);
	my $c_total = 0;
	my $c_alleles  = 0;
	my $i_total = 0;
	my $i_alleles = 0;
	foreach my $i (0..3){
	#	print "$c_counts[$i]\n";
		$c_total += $c_counts[$i];
		if ($c_counts[$i]){
			$c_alleles++;
		}
		$i_total += $i_counts[$i];
                if ($i_counts[$i]){
                        $i_alleles++;
                }
	}
#	print "$i_alleles\t$i_total\n";
	unless (($i_alleles > 2) or ($c_alleles > 2)){
		unless (($i_total < 30) or ($c_total < 30)){
			foreach my $i (0..3){
				if($c_counts[$i] > 0){
					my $p_i = ($i_counts[$i]/$i_total); 
					my $p_c = ($c_counts[$i]/$c_total);
					my $q_i = (1- $p_i); 
					my $q_c = (1- $p_c);
					my $dxy = (($p_i * $q_c) + ($q_i * $p_c));
					print "$a[0]\t$a[1]\t$dxy\n";
					last;
				}
			}
		}
	}
}
