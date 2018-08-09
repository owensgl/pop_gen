#!/bin/perl
use warnings;
use strict;

#This script takes two csvr files for two haplotypes in repulsion. It combines the information for both haplotypes and then interpolates between known points to try to maximize the amount of data. For use in Texanus BC1 mapping with extreme missing data.

my $csvr1 = $ARGV[0];
my $csvr2 = $ARGV[1];

my $output_csvr1 = $csvr1;
$output_csvr1 =~ s/.csvr/.imputed.csvr/g;
my $output_csvr2 = $csvr2;
$output_csvr2 =~ s/.csvr/.imputed.csvr/g;
my $output_conflicted = $csvr1;
$output_conflicted =~ s/.csvr/.conflicted.csvr/g;

my %reverse;
$reverse{"H"} = "A";
$reverse{"A"} = "H";
$reverse{"-"} = "-";

my %header;
my $headercounter = 0;
my %data;
my $nsamples;
my %loci_label;
#Load in data for the first haplotype, not reversed.
open CSVR1, $csvr1;
while(<CSVR1>){
  chomp;
  unless ($_ =~ m/^HanXRQ/){
    $header{$headercounter} = $_;
    $headercounter++;
  }else{
    my @a = split(/,/,$_);
    my $locus = $a[0];
    my $chr = $a[1];
    my $cm = $a[2];
    $nsamples = $#a;
    $loci_label{$chr}{$cm} = $locus;
    foreach my $i (3..$#a){
      if ($a[$i] eq "-"){next;}
      $data{$chr}{$cm}{$i} = $a[$i];
    }
  }
}
close CSVR1;

#Load in data for the second haplotype. This one reverse the genotypes
open CSVR2, $csvr2;
while(<CSVR2>){
  chomp;
  unless ($_ =~ m/^HanXRQ/){
  }else{
    my @a = split(/,/,$_);
    my $locus = $a[0];
    my $chr = $a[1];
    my $cm = $a[2];
    $loci_label{$chr}{$cm} = $locus;
    foreach my $i (3..$#a){
      if ($a[$i] eq "-"){next;}
      $data{$chr}{$cm}{$i} = $reverse{$a[$i]};
    }
  }
}
close CSVR2;
my %imputed_data;
my %conflicted_data;
foreach my $chr (sort keys %data){
  my @loci = (sort keys %{$data{$chr}});
  foreach my $sample (3..$nsamples){
    foreach my $n (1..($#loci -1)){ #not trying the edges because they can never be imputed.
      if ($data{$chr}{$loci[$n]}{$sample}){
    #    $imputed_data{$chr}{$loci[$n]}{$sample} = $data{$chr}{$loci[$n]}{$sample};
        next;
      }
      my $before_state;
      my $after_state;
      #Look backwards and get status
      foreach my $j (reverse 0..$n){
        if ($data{$chr}{$loci[$j]}{$sample}){
          $before_state = $data{$chr}{$loci[$j]}{$sample};
          goto NEXTLOOK_1;
        }
      }
      NEXTLOOK_1:
      foreach my $j ($n..$#loci){
        if ($data{$chr}{$loci[$j]}{$sample}){
          $after_state = $data{$chr}{$loci[$j]}{$sample};
          goto NEXTLOOK_2;
        }
      }
      NEXTLOOK_2:
      if (($after_state) and ($before_state)){
        if ($after_state eq $before_state){
          $imputed_data{$chr}{$loci[$n]}{$sample} = $after_state;
        }else{
          $conflicted_data{$chr}{$loci[$n]}{$sample}++;

        }
      }
    }
  }
}
#Print out imputed file 1.
open(my $out1, '>', $output_csvr1);
#print header
foreach my $n (0..($headercounter-2)){
  print $out1 "$header{$n}\n";
}
print $out1 "$header{($headercounter-1)}\n";

foreach my $chr (sort keys %imputed_data){
  my @loci = (sort keys %{$imputed_data{$chr}});
  foreach my $loci (@loci){
    print $out1 "\n$loci_label{$chr}{$loci},$chr,$loci";
    foreach my $sample (3..$nsamples){
      if ($imputed_data{$chr}{$loci}{$sample}){
        print $out1 ",$imputed_data{$chr}{$loci}{$sample}";
      }else{
        print $out1 ",-";
      }
    }
  }
}

#Print out imputed file 2.
open(my $out2, '>', $output_csvr2);
#print header
foreach my $n (0..($headercounter-2)){
  print $out2 "$header{$n}\n";
}
print $out2 "$header{($headercounter-1)}\n";

foreach my $chr (sort keys %imputed_data){
  my @loci = (sort keys %{$imputed_data{$chr}});
  foreach my $loci (@loci){
    print $out2 "\n$loci_label{$chr}{$loci},$chr,$loci";
    foreach my $sample (3..$nsamples){
      if ($imputed_data{$chr}{$loci}{$sample}){
        print $out2 ",$reverse{$imputed_data{$chr}{$loci}{$sample}}";
      }else{
        print $out2 ",-";
      }
    }
  }
}

#Print out imputed file 2.
open(my $out3, '>', $output_conflicted);
#print header
foreach my $n (0..($headercounter-2)){
  print $out3 "$header{$n}\n";
}
print $out3 "$header{($headercounter-1)}\n";

foreach my $chr (sort keys %imputed_data){
  my @loci = (sort keys %{$imputed_data{$chr}});
  foreach my $loci (@loci){
    print $out3 "\n$loci_label{$chr}{$loci},$chr,$loci";
    foreach my $sample (3..$nsamples){
      if ($conflicted_data{$chr}{$loci}{$sample}){
        print $out3 ",$conflicted_data{$chr}{$loci}{$sample}";
      }else{
        print $out3 ",-";
      }
    }
  }
}
