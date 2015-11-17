#!/bin/perl
use warnings;
use strict;

my %state_hash;
my %start_bp_hash;
my %start_cm_hash;
my $current_chrom;
my %samplelist;
my %specieslist;
my $last_position_bp;
my $last_position_cm;
while (<STDIN>){
	if ($. == 1){
		print "sample\tspecies\tchrom\tstate\tbp_start\tbp_end\tbp_size\tcm_start\tcm_end\tcm_size";
		next;
	}
	chomp;
	my @a = split(/\t/,$_);
	my $sample = $a[0];
	my $species = $a[1];
	my $chrom = $a[2];
	my $bp = $a[3];
	my $cm = $a[4];
	my $state = $a[5];
	unless ($samplelist{$sample}){
		$samplelist{$sample}++; #Make a list of all samples in a hash
		$specieslist{$sample} = $species;
	}
	unless($current_chrom){ #if it's a the start of the script
		$current_chrom = $chrom;
	}
	if ($current_chrom ne $chrom){ #If its a new chromosome
		foreach my $name (sort keys %samplelist){
			my $cm_size = $last_position_cm - $start_cm_hash{$name};
			my $bp_size = $last_position_bp - $start_bp_hash{$name};
			print "\n$name\t$specieslist{$name}\t$current_chrom\t$state_hash{$name}";
			print "\t$start_bp_hash{$name}\t$last_position_bp\t$bp_size";
			print "\t$start_cm_hash{$name}\t$last_position_cm\t$cm_size";
			delete($start_cm_hash{$name});
			delete($start_bp_hash{$name});
			delete($state_hash{$name});
		}
		$current_chrom = $chrom;
	}
        unless($state_hash{$sample}){ #if its the start of the script
                $start_cm_hash{$sample} = $cm;
                $start_bp_hash{$sample} = $bp;
                $state_hash{$sample} = $state;
        }
	if ($state_hash{$sample} ne $state){ #It's switched state
		my $cm_size = $cm - $start_cm_hash{$sample};
		my $bp_size = $bp - $start_bp_hash{$sample};
		print "\n$sample\t$specieslist{$sample}\t$chrom\t$state_hash{$sample}";
		print "\t$start_bp_hash{$sample}\t$bp\t$bp_size";
		print "\t$start_cm_hash{$sample}\t$cm\t$cm_size";
		#Start a new hash of infos
		$state_hash{$sample} = $state;
		$start_bp_hash{$sample} = $bp;
		$start_cm_hash{$sample} = $cm;
	}
	#save last positions for end of the chrom printing.
	$last_position_bp = $bp;
	$last_position_cm = $cm;
}

foreach my $name (sort keys %samplelist){ #if its the end of the file, print last windows                       
	my $cm_size = $last_position_cm - $start_cm_hash{$name};
       	my $bp_size = $last_position_bp - $start_bp_hash{$name};
	print "\n$name\t$specieslist{$name}\t$current_chrom\t$state_hash{$name}";
	print "\t$start_bp_hash{$name}\t$last_position_bp\t$bp_size";
	print "\t$start_cm_hash{$name}\t$last_position_cm\t$cm_size";
}
