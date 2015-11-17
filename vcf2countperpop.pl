#!/bin/perl
use warnings;
use strict;
use File::Basename;
my $min_gt_qual = 20;
my $min_mq = 20;
my $min_qual = 20;
my $min_dp = 5;
my $max_dp =100000;


my $in = $ARGV[0];
my %pophash;
my %poplist;
open IN, $in;
while (<IN>){
	chomp;
	my @a = split(/\t/,$_);
	$pophash{$a[0]} = $a[1];
	$poplist{$a[1]}++;
}
my @poplist = keys %poplist;
my %samplelist;
#old VCF2VERTICAL had format from mpileup: GT:PL:DP:GQ
#THIS VERSION:    has format from GATK-UG: GT:AD:DP:GQ:PL
#GLO VERSION Sept2014: includes ./.:
while(<STDIN>){
	if(eof()){
		#print "\n";
	}
	else{
		my %populations_sampled;
		my $line = "$_";
		chomp $line;
		my @fields = split /\t/,$line;
	    	if($line=~m/^##/){
			next;
		}
		elsif($fields[7]=~m/^NCC/) {
			next;
		}
		else{
			my $chrome = shift @fields;
			my $pos =    shift @fields;
			my $id =     shift @fields;
			my $ref =    shift @fields;
			my $alt =    shift @fields;
			my $qual =   shift @fields;
			my $filter = shift @fields;
			my $info =   shift @fields;
			my $format = shift @fields;
			my $mq = "NA";
			if($info=~m/MQ=(\d+)/){
				$mq = "$1";
			}
			my $meta = "$chrome\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";
			if($line=~m/^#/){
				foreach my $i (0..$#fields){
					$samplelist{$i} = $pophash{$fields[$i]};
				}
				print "chrom\tpos";
				foreach my $population (@poplist){
					print "\t$population";
				}
				next;
			}
            if ((length($ref) > 1) or (length($alt) > 1)){ #If its an indel, skip the line
                next;
            }
			unless ($info =~m/DP/){
				next;
			}
			if ($alt eq '.'){ #If there are no alternate allele
				print "\n$chrome\t$pos";
				foreach my $i (0..$#fields){
					my @tmp = split(/:/,$fields[$i]);
					my $depth = $tmp[1];
					if ($depth eq "."){
						$depth = 0;
					}
					if ($depth >= $min_dp){
						if ($samplelist{$i}){
							$populations_sampled{$samplelist{$i}}++;
						}
					}
				}
				foreach my $population (@poplist){
					if ($populations_sampled{$population}){
						print "\t$populations_sampled{$population}";
					}else{
						print "\t0";
					}
				}
			}
			else{
				print "\n$chrome\t$pos";
				foreach my $i (0..$#fields){
					my @tmp = split(/:/,$fields[$i]);
					my $depth = $tmp[2];
					if ($depth eq "."){
                                                $depth = 0;
                                        }
					if ($depth >= $min_dp){
						if($samplelist{$i}){
							$populations_sampled{$samplelist{$i}}++;
						}
					}
				}
				foreach my $population (@poplist){
					if ($populations_sampled{$population}){
						print "\t$populations_sampled{$population}";
					}else{
						print "\t0";
					}
				}
			}
		}
	}
}
