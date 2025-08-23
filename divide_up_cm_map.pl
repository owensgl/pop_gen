#!/bin/perl
#This script adds spaces to genetic maps with identical marker locations.

my $previous_cm = "NA";
my $previous_line;
my @stored_lines;
while (<STDIN>){
    chomp;
    if ($. == 1){
        print "$_";
    }else{
        my @a = split(/\t/,$_);
        my $cm = $a[2];
        my $line = "$a[0]\t$a[1]";
        if ($cm != $previous_cm){ #If this is different than the previous cm. Print stored lines, clear stored lines, and then store this one.
            my $divisions = $#stored_lines + 2;
            my $floor = $previous_cm - 0.005;
            my $unit = 0.01 / $divisions;
            foreach my $i (0..$#stored_lines){
                my $new_cm = $floor + ($unit * ($i+1));
                print "\n$stored_lines[$i]\t$new_cm";
            }
            undef(@stored_lines);
            push(@stored_lines,$line);
            $previous_cm = $cm;
        }else{
            #Store lines
            push(@stored_lines,$line);

        }
    }
}
