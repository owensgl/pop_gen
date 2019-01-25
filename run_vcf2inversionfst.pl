#!/usr/bin/perl
#This runs multiple instances of vcf2inversionfst.pl for each mds outlier selected. 
use warnings;
use strict;
use Parallel::ForkManager;

my $pm = new Parallel::ForkManager(4);
my $genotype_file = $ARGV[0]; #This is the Ha412HO_inv.dec11.inversions.genotypes.v1.txt file 
my $vcf_script = "perl /home/owens/bin/pop_gen/vcf2inversionfst.pl";
my $prefix = "/media/owens/Copper/wild_gwas_2018/";
my @outputname = split(/\//,$genotype_file);
my $outputname = $outputname[$#outputname];
$outputname =~ s/.inversions.genotypes.v.*.txt//g;

my $file_start = "tranche90.snp.";
my $file_end = ".90.bi.remappedHa412HO.vcf.gz";

open GENO, $genotype_file;

my %mds;
my %geno;
my %chr;
my %species;

my %fullspecies;
$fullspecies{"annuus"} = "Annuus";
$fullspecies{"petpet"} = "Petiolaris";
$fullspecies{"petfal"} = "Petiolaris";
$fullspecies{"argophyllus"} = "Argophyllus";

while(<GENO>){
	chomp;
	if ($. == 1){next;}
	my @a = split(/\t/,$_);
	my $mds = "$a[4].$a[5].$a[6]";
	my $sample = $a[0];
	my $genotype = $a[3];
	my $chr = $a[5];
	my $species = $a[4];
	$mds{$mds}++;
	$geno{$mds}{$sample}=$genotype;
	$chr{$mds} = $chr;
	$species{$mds} = $species;
}

foreach my $mds (sort keys %mds){
	my $filename = "TMP.$mds.genolist.txt";
	open (my $fh1,'>', $filename);
	my $first_line;
	foreach my $sample (sort keys %{$geno{$mds}}){
		unless($first_line){
			print $fh1 "$sample\t$geno{$mds}{$sample}";
			$first_line++;
		}else{
			print $fh1 "\n$sample\t$geno{$mds}{$sample}";
		}
	}
	close $fh1;
}
foreach my $mds (sort keys %mds){
	$pm->start and next;
	#dataset name in the middle of the vcf file name needs to be selected; Custom hacky script.
	my $middle = $species{$mds}; 
	if ($species{$mds} eq "annuus"){
		$middle = "env";
	}elsif($species{$mds} eq "argophyllus"){
		$middle = "gwas";
	}
	my $gzvcf_file = "${prefix}$species{$mds}/$fullspecies{$species{$mds}}.$file_start$middle$file_end";
	my $command = "bcftools view -r $chr{$mds} -O v $gzvcf_file | $vcf_script TMP.$mds.genolist.txt |gzip  > $prefix$species{$mds}/$outputname.$mds.fst.txt.gz";
	print("$command\n");
	system($command);
        system("rm TMP.$mds.genolist.txt");
	$pm->finish;
}
$pm->wait_all_children;




