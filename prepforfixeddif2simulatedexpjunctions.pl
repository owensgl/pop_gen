#!/bin/perl
use warnings;
use strict;
use POSIX;
use Statistics::Distribution::Generator qw( :all );
#This script takes the output from "/home/owens/bin/reformat/ctab2prepforparentblocks.pl" and simulates samples with different junction densities, with an exponential distribution. It uses the missing data patterns of my real samples and also simulates wrongly called identity using a bayesian strategy.

my $min_junctions_theta = 250;
my $max_junctions_theta = 20000;
my $increment = 250;
my $rep = 20; #Reps per species;
my $error_multiplier = $ARGV[0];

my @species_list = ("Ano","Des","Par");
my %error;
$error{5}=0.005893446;
$error{6}=0.003928964;
$error{7}=0.002750274;
$error{8}=0.002000199;
$error{9}=0.001500149;
$error{10}=0.001153961;
$error{11}=0.0009066832;

foreach my $key (sort keys %error) {
    $error{$key} = $error{$key}* $error_multiplier;
}

my %species;
$species{"Des1484"}="Des";
$species{"des2458"}="Des";
$species{"Sample_DES1476"}="Des";
$species{"Ano1495"}="Ano";
$species{"Sample_Ano1506"}="Ano";
$species{"Sample_des1486"}="Des";
$species{"Sample_Des2463"}="Des";
$species{"Sample_desA2"}="Des";
$species{"Sample_desc"}="Des";
$species{"king141B"}="Par";
$species{"king145B"}="Par";
$species{"king147A"}="Par";
$species{"King151"}="Par";
$species{"king152"}="Par";
$species{"King156B"}="Par";
$species{"Sample_king1443"}="Par";
$species{"Sample_king159B"}="Par";

my %name;
my %data;
my %site;
my %location;
my $counter;
my %err_hash;
my %chromosomes;
#Load in site data including where missing data is.
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (4..$#a){
      $name{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $cm = $a[2];
    my $current_error = $error{$a[3]};
    if ($chr =~ m/Chr00/){next;}
    $counter++;
    $err_hash{$counter} = $current_error;
    $site{$counter} = $chr;
    $location{$counter} = $cm;
    $chromosomes{$chr}++;
    foreach my $i (4..$#a){
      if ($a[$i] eq "N"){
        $data{$name{$i}}{$counter} = 0;
      }else{
        $data{$name{$i}}{$counter} = 1;
      }
    }  
  }
}
my $max_chr_length = 104;
foreach my $current_species (@species_list){
  for (my $density = $min_junctions_theta; $density <=$max_junctions_theta; $density += $increment){
    foreach my $template (sort keys %species){
      if ($species{$template} ne $current_species){
        next;
      }
      for (my $j = 1; $j <= $rep; $j+=1){
	#Define the boundaries of the simulated blocks.
	my %sim_junctions;
#	my $start_run = time();
	foreach my $chr (sort keys %chromosomes){
 	  my $current_spot = 0;
          my $current_state = int(rand(2)) * 2;
	  until($current_spot == $max_chr_length){
	    my $increment = exponential($density);
	    my $end = $current_spot+$increment;
	    my $index = ceil($end);
	    $sim_junctions{$chr}{$index}{$end}=$current_state;
	    if ($current_state == 2){
	      $current_state = 0;
	    }else{
	      $current_state = 2;
	    }
            $current_spot = $end;
            if ($current_spot > $max_chr_length){
              $current_spot = $max_chr_length;
	    }
          }
        }
#	my $end_run = time();
#	my $run_time = $end_run - $start_run;
#	print STDERR "Making exponential junctions took $run_time seconds\n";
        my $junc_counter;
        until($species{$template} eq $current_species){
          $template = (keys %species)[rand keys %species];
        }
        #Keep track of parentage using even and odd divisions of the increment;
        my $current_state;
        my $current_chr;
#	$start_run = time();
        foreach my $n (1..$counter){
          my $chr = $site{$n};
          unless($current_chr){
            $current_chr = $chr;
          }
          if ($current_chr ne $chr){
            undef($current_state);
            $current_chr = $chr;
          }
          
          if ($data{$template}{$n}){ #Only continue if it's got data in the template
            my $cm = $location{$n};
            my $true_state;
	    my $index = ceil($cm);
	    foreach my $spot (sort {$a <=> $b} keys %{$sim_junctions{$chr}{$index}}){
	      if ($spot >= $cm){
                $true_state = $sim_junctions{$chr}{$index}{$spot};
		goto NEXTSPOT;
	      }
            }
	    foreach my $plus (1..50){
	      my $new_index = $index + $plus;
              foreach my $spot (sort {$a <=> $b} keys %{$sim_junctions{$chr}{$new_index}}){
                if ($spot >= $cm){
                  $true_state = $sim_junctions{$chr}{$new_index}{$spot};
                  goto NEXTSPOT;
                }
	      }
            }
            NEXTSPOT:
unless(defined $true_state){print STDERR "No true state: cm = $cm, chr = $chr"};
            #print STDERR "\n$cm\t$window\t$true_state";
#print "\n$cm\tTrue state is $true_state";
            my $accurate = 0; #Check to see if marker is randomly wrong.
    	    my $rand = rand(1);
            if ($rand > $err_hash{$n}){
              $accurate = 1;
            }
#print "\taccuracy = $accurate";
            my $viewed_state;
            if ($accurate){
              $viewed_state = $true_state;
            }else{
              if ($true_state == 0){
                $viewed_state = 2;
              }else{
                $viewed_state = 0;
              }
            }
#print "\tViewed_state = $viewed_state";
            if (defined $current_state){
  	      my $added_junc = abs($viewed_state - $current_state)/2;
              $junc_counter += $added_junc;
#print "\t+$added_junc";
              $current_state = $viewed_state;
            }else{
              $current_state = $viewed_state;
            }
          }
        }
#	$end_run = time();
#        $run_time = $end_run - $start_run;
#	print STDERR "Counting junctions took $run_time seconds\n";
        print "\n$current_species\t$template\t$density\t$junc_counter";
      }
    }
  }
}

