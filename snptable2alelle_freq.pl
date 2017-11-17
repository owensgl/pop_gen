#!/bin/perl
use warnings;
use strict;

##This takes a tab delimited file and calculates allele frequency for populations. It also takes a second file that has "groups", which are made up of one or more populations. It calculates the allele frequency for both populations and any groups
#usage: cat SNPTABLE.tab | perl this_script.pl sampleinfo.txt groupinfo.txt > output.txt
#sampleinfo is two columns, first column sample name, second column population name
#groupinfo is two columns, first column population name, second column group name

my $min_size = 4; #minimum number of samples needed to calculate an allele frequency. 

my $samplefile = $ARGV[0];

my %pop;
my %poplist;
my %group;
my %grouplist;
open SAMPLE, $samplefile;
while(<SAMPLE>){
  chomp;
  my @a = split(' ',$_);
  my $name = $a[0];
  my $pop = $a[1];
  $pop{$name} = $pop;
  $poplist{$pop}++;
}
close SAMPLE;
if ($ARGV[1]){
  my $metagroupfile = $ARGV[1];
  open GROUP, $metagroupfile;
  while(<GROUP>){
    chomp;
    my @a = split(' ',$_);
    my $pop = $a[0];
    my $group = $a[1];
    $group{$pop} = $group;
    $grouplist{$group}++;
  }
}
my %name;

print "chr\tpos";
foreach my $pop (sort keys %poplist){
  print "\t$pop";
}
if ($ARGV[1]){
  foreach my $group (sort keys %grouplist){
    print "\t$group";
  }
}
while(<STDIN>){
  chomp;
  my @a = split(' ',$_);
  if ($. == 1){
    foreach my $i (2..$#a){
      $name{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $pos = $a[1];
    my %bases;
    foreach my $i (2..$#a){
      if ($a[$i] eq "NN"){next;}
      my @alleles = split(//,$a[$i]);
      foreach my $j (0..1){
        $bases{$alleles[$j]}++;
      }
    }
    my @ordered_bases = sort { $bases{$b} <=> $bases{$a} } keys %bases;
    my $n_alleles = scalar @ordered_bases;
    if ($n_alleles != 2){next;} #Skip not biallelic sites.
    my $major = $ordered_bases[0];
    my %major_count;
    my %total_count;
    foreach my $i (2..$#a){
      if ($a[$i] eq "NN"){next;}
      if ($pop{$name{$i}}){
        my @alleles = split(//,$a[$i]);
        foreach my $j (0..1){
          if ($alleles[$j] eq $major){
            $major_count{$pop{$name{$i}}}++;
            $total_count{$pop{$name{$i}}}++;
          }else{
            $total_count{$pop{$name{$i}}}++;
          }
        }
      }
      if ($group{$pop{$name{$i}}}){
        my @alleles = split(//,$a[$i]);
        foreach my $j (0..1){
          if ($alleles[$j] eq $major){
            $major_count{$group{$pop{$name{$i}}}}++;
            $total_count{$group{$pop{$name{$i}}}}++;
          }else{
            $total_count{$group{$pop{$name{$i}}}}++;
          }
        }
      }  
    }
    print "\n$chr\t$pos";
    foreach my $pop (sort keys %poplist){
      unless ($total_count{$pop}){
        print "\tNA";
      }elsif ($total_count{$pop} < $min_size){
        print "\tNA";
      }else{
        unless ($major_count{$pop}){
          $major_count{$pop} = 0;
        }
        my $freq = $major_count{$pop} / $total_count{$pop};
        print "\t$freq";
      }
    }
    if ($ARGV[1]){
      foreach my $group (sort keys %grouplist){
	unless ($total_count{$group}){
          print "\tNA";
        }elsif ($total_count{$group} < $min_size){
          print "\tNA";
        }else{
          unless ($major_count{$group}){
            $major_count{$group} = 0;
          }
          my $freq = $major_count{$group} / $total_count{$group};
          print "\t$freq";
        }
      }
    }
  }
}
