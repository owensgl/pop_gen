#!/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);
my $samplefile = $ARGV[0];
#This script jackknifes the Neis D estimates by dropping each sample and calculating D for the species without it.
my %printlist;
my $snpfile = "germplasm.freebayes.full.tab";
my $nei_script = "/home/owens/bin/pop_gen/SNPtable2NeiD.pl";
my $outfile = "germplasm.freebayes.full.NeiD.jackknife.tab";
my %specieslist;
my @species_array;
my @sample_list;
system("rm $outfile");
open SAMPLEFILE, $samplefile;
while(<SAMPLEFILE>){
	chomp;
	my @a = split(/\t/,$_);
	$printlist{$a[0]}++;
	$specieslist{$a[0]} = $a[1]; #should have real species names, don't include species you don't want compared.
	push(@species_array,$a[1]);
	push(@sample_list,$a[0]);
}
close SAMPLEFILE;

my @tmp = uniq @species_array;
@species_array = @tmp;

foreach my $focal_sample (@sample_list){
	my $current_species = $specieslist{$focal_sample};
	open (my $file, '>', "tmp_samplelist.txt"); 
	foreach my $sample (@sample_list){
		if (($sample ne $focal_sample) and ($specieslist{$sample} eq $current_species)){
			print $file "$sample\t$specieslist{$sample}\n";
		}
	}
	close $file;
	system("cat $snpfile | perl $nei_script tmp_samplelist.txt >> $outfile");	
#	exit;
}
