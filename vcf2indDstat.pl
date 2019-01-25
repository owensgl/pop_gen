#!/bin/perl
use warnings;
use strict;
#This is running an abba baba test on individuals. It looks for sites where the outgroup samples are monomorphic. For each individual in group 2 it randomly chooses one allele. Then it selects a single allele from group 1 and group 3, and outputs the pattern. If there are not only two alleles, it skips the site for the individual.
my %h;
$h{"AAAA"} = 1;
$h{"BAAA"} = 2;
$h{"ABAA"} = 3;
$h{"AABA"} = 4;
$h{"BBAA"} = 5;
$h{"BABA"} = 6;
$h{"ABBA"} = 7;
$h{"BBBA"} = 8;

my $info_file = $ARGV[0];
my $group_file = $ARGV[1];

open INFO, $info_file;
my $min_in_group = 4; #minimum number of indviduals genotyped in group 1 and 3. 
my $min_out_group = 4; #minimum number of individuals genotyped in group 4.
my %pop;
my %species;
my %species_id;
my %bam;
while(<INFO>){
  chomp;
  my @a = split(/\t/,$_);
  my $bam = $a[0];
  my $sample = $a[1];
  my $pop = $a[3];
  my $species = $a[2];
  my $type = $a[4];
  if ($type ne "wild"){next;} #Use this to avoid domestic samples;
  $species_id{$sample} = $species;
  $pop{$sample} = $pop;
}
close INFO;
open GROUPS, $group_file;
my %group;
while(<GROUPS>){
  chomp;
  my @a = split(/\t/,$_);
  $group{$a[0]} = $a[1];
}
close GROUPS;

my %sample;
my $header;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^##/){
    next;
  }
  if ($_ =~ m/^#/){
    my @a = split(/\t/,$_);
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    unless($header){
      print "chr\tpos";
      foreach my $i (9..$#a){
        unless($species_id{$sample{$i}}){next;}
	unless($group{$species_id{$sample{$i}}}){next;}
        if ($group{$species_id{$sample{$i}}} == 1){
          print "\t$sample{$i}";
        }
      }
      $header++;
    }
  }else{
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    my $bp = $a[1];
    my %data;
    my %group_data;
    my %group_freq;
    my %group_count;
    foreach my $i (9..$#a){
      unless($species_id{$sample{$i}}){next;}
      unless($group{$species_id{$sample{$i}}}){next;} #Skip samples without a group
      my @fields = split(/:/,$a[$i]);
      if (($fields[0] eq '.') or ($fields[0] eq './.')){
        next;
      }
      my @genotypes = split('/',$fields[0]);
      if ($group{$species_id{$sample{$i}}} == 1){
        $data{$sample{$i}}{0} = $genotypes[0];
        $data{$sample{$i}}{1} = $genotypes[1];
      }
      push (@{$group_data{$group{$species_id{$sample{$i}}}}}, $genotypes[0]);
      push (@{$group_data{$group{$species_id{$sample{$i}}}}}, $genotypes[1]);
      $group_freq{$group{$species_id{$sample{$i}}}}{$genotypes[0]}++;
      $group_freq{$group{$species_id{$sample{$i}}}}{$genotypes[1]}++;
      $group_count{$group{$species_id{$sample{$i}}}}++;
    }
    #Check to make sure it is monomorphic in group 4 and has enough genotypes.
    unless(($group_count{4}) and ($group_count{1}) and ($group_count{3})){next;}
    if ($group_count{4} < $min_out_group){next;}
    my @outgroup_alleles = keys %{$group_freq{4}};
    if (scalar(@outgroup_alleles) != 1){next;}

    #Check to make sure group 1 and group 3 have enough genotypes.
    if (($group_count{1} < $min_in_group) or ($group_count{3} < $min_in_group)){next;}
    #For each individual in group 1, output a randomly chosen state.
    print "\n$chr\t$bp";
    foreach my $i (9..$#a){
      unless($species_id{$sample{$i}}){next;}
      unless($group{$species_id{$sample{$i}}}){next;}
      if ($group{$species_id{$sample{$i}}} == 1){
        unless (exists $data{$sample{$i}}{0}){
          print "\t0";
          next;
        }
        my $ancestor = $outgroup_alleles[0];
        my $group_2_genotype = $data{$sample{$i}}{int(rand(2))};
        my @data = keys %data;
        my $group_1_ind = $sample{$i};
        until (($group_1_ind ne $sample{$i})and (exists $data{$group_1_ind}{0})){
          $group_1_ind = $data[int rand( @data)]; 
        }
        my $group_1_genotype = $data{$group_1_ind}{int(rand(2))};
        my $tmp_length = @{ $group_data{3} };
	my $group_3_genotype  = $group_data{3}[rand $tmp_length];
        my %chosen_alleles;
        $chosen_alleles{$ancestor}++;
        $chosen_alleles{$group_1_genotype}++;
        $chosen_alleles{$group_2_genotype}++;
        $chosen_alleles{$group_3_genotype}++;
        my $n_alleles = keys %chosen_alleles;
        if ($n_alleles == 1){
          print "\t1"; #AAAA
          next;
        }elsif ($n_alleles > 2){
          print "\t0"; #More than two alleles
          next;
        }
        my $pattern;
        if ($group_1_genotype eq $ancestor){
          $pattern.="A";
        }else{
          $pattern.="B";
        }
        if ($group_2_genotype eq $ancestor){
          $pattern.="A";
        }else{
          $pattern.="B";
        }
        if ($group_3_genotype eq $ancestor){
          $pattern.="A";
        }else{
          $pattern.="B";
        }
        $pattern.="A";
        print "\t$h{$pattern}";
      }
    }
  }
}
