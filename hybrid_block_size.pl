#!/bin/perl

use warnings;
use strict;

my $in = $ARGV[0]; #Summary hybrid windows
open IN, $in;

my %h;
while (<IN>){
	chomp;
	next if $. == 1;
	my @a = split(/\t/,$_);
	my $sample = $a[0];
	my $chrom = $a[1];
	my $start = $a[2];
	my $end = $a[3];
	my $cm_start = $a[4];
	my $cm_end = $a[5];
	my $cm_size = $a[6];
	if ($cm_size eq "NA"){
		next;
	}
	my $status = $a[11];
	$h{$sample}{$chrom}{$start}{'status'} = $status;
	$h{$sample}{$chrom}{$start}{'cm_start'} = $cm_start;
	$h{$sample}{$chrom}{$start}{'cm_end'} = $cm_end;
	$h{$sample}{$chrom}{$start}{'cm_size'} = $cm_size;
	$h{$sample}{$chrom}{$start}{'end'} = $end;
}

foreach my $sample (sort keys %h){
	foreach my $chrom (sort keys %{$h{$sample}}){
		my $ongoingP1;
		my $ongoingP2;
		my $startP1;
		my $endP1;
		my $startP2;
		my $endP2;
		my $emptywindow;
		my @start = (sort {$a<=>$b} keys %{$h{$sample}{$chrom}});
		foreach my $i (0..$#start){
		#	print "$h{$sample}{$chrom}{$start[$i]}{'cm_end'}\n";
			if ($h{$sample}{$chrom}{$start[$i]}{'status'} eq "P1"){
				unless($ongoingP1){
					$startP1 = $start[$i];
				}
				if ($ongoingP2){
					$endP2 = $start[$i];
					my $cm_size = $h{$sample}{$chrom}{$endP2}{'cm_end'} - $h{$sample}{$chrom}{$startP2}{'cm_start'};
					print "\n$sample\t$chrom\t$startP2\t$endP2\t$h{$sample}{$chrom}{$startP2}{'cm_start'}\t$h{$sample}{$chrom}{$endP2}{'cm_end'}\t$cm_size\tP2";
					undef $startP2;
					undef $endP2;
					undef $ongoingP2;
				}
				$ongoingP1++;
				undef $emptywindow;
			}elsif ($h{$sample}{$chrom}{$start[$i]}{'status'} eq "P2"){
				unless($ongoingP2){
                                        $startP2 = $start[$i];
                                }
                                if ($ongoingP1){
                                        $endP1 = $start[$i];
                                        my $cm_size = $h{$sample}{$chrom}{$endP1}{'cm_end'} - $h{$sample}{$chrom}{$startP1}{'cm_start'};
                                        print "\n$sample\t$chrom\t$startP1\t$endP1\t$h{$sample}{$chrom}{$startP1}{'cm_start'}\t$h{$sample}{$chrom}{$endP1}{'cm_end'}\t$cm_size\tP1";
                                        undef $startP1;
                                        undef $endP1;
                                        undef $ongoingP1;
                                }
                                $ongoingP2++;
				undef $emptywindow;
			}elsif ($h{$sample}{$chrom}{$start[$i]}{'status'} eq "NA"){
				if ($emptywindow){
					if($ongoingP1){
						$endP1 = $start[$i];
	                                        my $cm_size = $h{$sample}{$chrom}{$endP1}{'cm_end'} - $h{$sample}{$chrom}{$startP1}{'cm_start'};
        	                                print "\n$sample\t$chrom\t$startP1\t$endP1\t$h{$sample}{$chrom}{$startP1}{'cm_start'}\t$h{$sample}{$chrom}{$endP1}{'cm_end'}\t$cm_size\tP1";
                	                        undef $startP1;
                        	                undef $endP1;
                                	        undef $ongoingP1;
					}elsif($ongoingP2){
	                                        $endP2 = $start[$i];
        	                                my $cm_size = $h{$sample}{$chrom}{$endP2}{'cm_end'} - $h{$sample}{$chrom}{$startP2}{'cm_start'};
                	                        print "\n$sample\t$chrom\t$startP2\t$endP2\t$h{$sample}{$chrom}{$startP2}{'cm_start'}\t$h{$sample}{$chrom}{$endP2}{'cm_end'}\t$cm_size\tP2";
                        	                undef $startP2;
                                	        undef $endP2;
                                       		undef $ongoingP2;
					}
				}
				$emptywindow++;
			}elsif ($h{$sample}{$chrom}{$start[$i]}{'status'} eq "admixed"){
                        	if($ongoingP1){
                                	$endP1 = $start[$i];
                                        my $cm_size = $h{$sample}{$chrom}{$endP1}{'cm_end'} - $h{$sample}{$chrom}{$startP1}{'cm_start'};
                                        print "\n$sample\t$chrom\t$startP1\t$endP1\t$h{$sample}{$chrom}{$startP1}{'cm_start'}\t$h{$sample}{$chrom}{$endP1}{'cm_end'}\t$cm_size\tP1";
                                        undef $startP1;
                                        undef $endP1;
                                        undef $ongoingP1;
                                }elsif($ongoingP2){
                                        $endP2 = $start[$i];
                                        my $cm_size = $h{$sample}{$chrom}{$endP2}{'cm_end'} - $h{$sample}{$chrom}{$startP2}{'cm_start'};
                                        print "\n$sample\t$chrom\t$startP2\t$endP2\t$h{$sample}{$chrom}{$startP2}{'cm_start'}\t$h{$sample}{$chrom}{$endP2}{'cm_end'}\t$cm_size\tP2";
                                        undef $startP2;
                                        undef $endP2;
                                        undef $ongoingP2;
                                }
				undef $emptywindow;			
			}
		}
                if($ongoingP1){
                        $endP1 = $start[$#start];
                        my $cm_size = $h{$sample}{$chrom}{$endP1}{'cm_end'} - $h{$sample}{$chrom}{$startP1}{'cm_start'};
                        print "\n$sample\t$chrom\t$startP1\t$endP1\t$h{$sample}{$chrom}{$startP1}{'cm_start'}\t$h{$sample}{$chrom}{$endP1}{'cm_end'}\t$cm_size\tP1";
               	}elsif($ongoingP2){
               		$endP2 = $start[$#start];
                        my $cm_size = $h{$sample}{$chrom}{$endP2}{'cm_end'} - $h{$sample}{$chrom}{$startP2}{'cm_start'};
                       	print "\n$sample\t$chrom\t$startP2\t$endP2\t$h{$sample}{$chrom}{$startP2}{'cm_start'}\t$h{$sample}{$chrom}{$endP2}{'cm_end'}\t$cm_size\tP2";
		}
	}
}
