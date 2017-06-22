#!/bin/perl
use warnings;
use strict;
use POSIX;

#This script takes the output from "/home/owens/bin/reformat/ctab2prepforparentblocks.pl" and simulates samples with different junction densities. It uses the missing data patterns of my real samples and also simulates wrongly called identity using a bayesian strategy.
#It outputs per X cM window. It also runs simulations until they pass the 90% confidence interval of the real junction number. 
my $window_size = 10; #in cM
my $min_junctions_per_cm = 1;
my $max_junctions_per_cm = 10000;
my $increment = 10;
my $rep = 100; #Reps per species;
my $total_cm = 1394.91; #For XRQ annuus genome.
my $error_multiplier = 1;

my @species_list = ("Ano","Des","Par");
my %error;
$error{5}=0.007567428 * 2;
$error{6}=0.005486178 * 2;
$error{7}=0.004158341 * 2;
$error{8}=0.003259651 * 2;
$error{9}=0.002623258 * 2;
$error{10}=0.002156185 * 2;
$error{11}=0.00180328 * 2;

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
my %current_state;
my $current_chr;
my %junctions;
my $window_end = $window_size;
my $total_columns;
my %marker_counter;
#Load in site data including where missing data is.
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (4..$#a){
      $name{$i} = $a[$i];
    }
    $total_columns = $#a;
    print "type\tspecies\tsample\tchr\twindow_end\tmarkers\tdensity\tjunctions";
  }else{
    my $chr = $a[0];
    my $cm = $a[2];
    my $current_error = $error{$a[3]};
    if ($chr =~ m/Chr00/){next;}
    unless($current_chr){
      $current_chr = $chr;
    }
    if ($current_chr ne $chr){
      &print_junction_counts();
      &simulate_junctions();
      undef(%current_state);
      undef($counter);
      undef(%data);
      undef(%site);
      undef(%err_hash);
      undef(%location);
      undef(%junctions);
      undef(%marker_counter);
      $current_chr = $chr;
      $window_end = $window_size;
      until($cm < $window_end){
        $window_end+=$window_size;
      }
    }
    if ($cm > $window_end){
      &print_junction_counts();
      &simulate_junctions();
      undef(%current_state);
      undef($counter);
      undef(%data);
      undef(%site);
      undef(%err_hash);
      undef(%location);
      undef(%junctions);
      undef(%marker_counter);
      until($cm < $window_end){
        $window_end+=$window_size;
      }
    }
    $counter++;
    $err_hash{$counter} = $current_error;
    $site{$counter} = $chr;
    $location{$counter} = $cm;
    foreach my $i (4..$#a){
      if ($a[$i] eq "N"){
        $data{$name{$i}}{$counter} = 0;
      }else{
        $marker_counter{$name{$i}}++;
        $data{$name{$i}}{$counter} = 1;
	if(defined $current_state{$name{$i}}){ #If it has a previous state.
          my $state_dif = abs($current_state{$name{$i}} - $a[$i])/2;
          $junctions{$name{$i}}+=$state_dif;
          $current_state{$name{$i}} = $a[$i];

        }else{
          $current_state{$name{$i}} = $a[$i];
        }
      }
    }  
  }
}
#Print real junction counts for the window;
sub print_junction_counts{
  foreach my $i (4..$total_columns){
    unless (defined $junctions{$name{$i}}){
      $junctions{$name{$i}} = 0;
    }
    unless (defined $marker_counter{$name{$i}}){
      $marker_counter{$name{$i}} = 0;
    }
    print "\nmeasured\t$species{$name{$i}}\t$name{$i}\t$current_chr\t$window_end\t$marker_counter{$name{$i}}\tNA\t$junctions{$name{$i}}";
   }
}

#Simulate increasing junction density. Make a note when the 90% percentile of the simulations overlap withe real value and stop simulating after then 10% percentile is above the real value;
sub simulate_junctions{
  foreach my $i (4..$total_columns){
    my $template = $name{$i};
    my $current_species = $species{$name{$i}};
    my $finished_sims;
    my $overlap_empirical;
    my $density = $min_junctions_per_cm;
    until(defined $finished_sims){
      my $junc_dist = 1 / $density;
      my @sim_junctions;
      for (my $j = 1; $j <= $rep; $j+=1){
        my $junc_counter = 0;
        my $start = rand($junc_dist);
        #Keep track of parentage using even and odd divisions of the increment;
        my $current_state;
        foreach my $n (1..$counter){
          
          if ($data{$template}{$n}){ #Only continue if it's got data in the template
            my $cm = $location{$n};
            my $true_state;
            if ($cm < $start){
              $true_state = 0;
            }else{
              my $window = floor(($cm - $start)/$junc_dist);
              if ($window % 2 == 0){
                $true_state = 2;
              }else{
                $true_state = 0;
              }
#print STDERR "\n$cm\t$window\t$true_state";
            }
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
  	push(@sim_junctions,$junc_counter);	
        print "\nsim\t$current_species\t$template\t$current_chr\t$window_end\t$marker_counter{$template}\t$density\t$junc_counter";
      }
      #Calculate top and bottom 10% values of junction counts
      my @sorted_sim_junctions = sort { $a <=> $b } @sim_junctions;
      my $top_spot = ceil($#sorted_sim_junctions * 0.9);
      my $bottom_spot = floor($#sorted_sim_junctions * 0.1);
      #if the 90% percentile overlaps with the real value, print that out for plotting.
      unless($overlap_empirical){
        if ($sorted_sim_junctions[$top_spot] > $junctions{$template}){
          $overlap_empirical++;
          print "\nbottom_limit\t$current_species\t$template\t$current_chr\t$window_end\t$marker_counter{$template}\t$density\tNA";
        }
      }
      if ($sorted_sim_junctions[$bottom_spot] > $junctions{$template}){
        print "\nupper_limit\t$current_species\t$template\t$current_chr\t$window_end\t$marker_counter{$template}\t$density\tNA";
        $finished_sims++;
      }
      if ($density >= $max_junctions_per_cm){
        $finished_sims++;
	print "\nupper_limit\t$current_species\t$template\t$current_chr\t$window_end\t$marker_counter{$template}\t$density\tNA";
      }
      $density+=$increment;
    }
  }
}

