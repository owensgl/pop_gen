#!/usr/bin/perl
#Sliding window for Fst and dxy 
#V3 introduces variable window size
use List::MoreUtils qw(uniq);
use warnings;

my %h;
my %poslist;
my @chrlist;
my $Start = 0;
my $Size = $ARGV[0];
my @statlist = ("Hexp");
my @colname;
my %count;

while (<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	my $chr;
	my $pos;
	my %stat_hash;
	if ($. == 1){
		foreach my $i (0..$#a){
			$colname[$i] = $a[$i];
		}
	}else{ #Load in each variable as well as a hash of positions and an array of chromosomes.
		foreach my $i (0..$#a){
			if ($colname[$i] eq "CHROM"){
				$chr = $a[$i];
			}elsif ($colname[$i] eq "POS"){
				$pos = $a[$i];
			}elsif ($colname[$i] eq "Dxy"){
				$stat_hash{'dxy'} = $a[$i];
			}elsif ($colname[$i] eq "FstNum"){
				$stat_hash{'FstNum'} = $a[$i];
			}elsif ($colname[$i] eq "FstDenom"){
				$stat_hash{'FstDenom'} = $a[$i];
			}elsif ($colname[$i] eq "Hexp"){
				$stat_hash{'Hexp'} = $a[$i];
			}elsif ($colname[$i] eq "N1"){
				$stat_hash{'N1'} = $a[$i];
			}
		}	#Take each stat and put it into a hash
		next if ($stat_hash{'N1'} eq "0"); 
		$h{'usablecount'}{$chr}{$pos}++;
		push (@chrlist, $chr);
		$poslist{$chr}{$pos}++;
		if ($stat_hash{'Hexp'} ne "0"){
			$h{'variablecount'}{$chr}{$pos}++;
		}
		foreach my $stat (@statlist){
			$h{$stat}{$chr}{$pos} = $stat_hash{$stat};
			$count{$stat}{$chr}{$pos}++;
		}
	}
}	

my @uniq_chr = uniq(@chrlist);
print "Chr\tStartPos\tEndPos\tMidPos\tTotalSites\tVariableSites\tHexp\n";
foreach my $chr (@uniq_chr){ #For each chromosome
	my $s = $Start; #Start position of scan
	my $e = $Size; #Size of each window
	my %total;
	my $count;
	foreach my $pos (sort {$a<=>$b} keys %{$poslist{$chr}}){
		if (($pos >= $s)and($pos<$e)){
			$count++;
			foreach my $stat(@statlist){
				if($h{$stat}{$chr}{$pos}){
					if ($h{$stat}{$chr}{$pos} ne "NA"){
						if ($total{$stat}){
							$total{$stat} += $h{$stat}{$chr}{$pos};
						}else{
							$total{$stat} = $h{$stat}{$chr}{$pos};
						}
					}
				}
			}
			if ($h{'usablecount'}{$chr}{$pos}){
				if ($total{'usablecount'}){
					$total{'usablecount'} += $h{'usablecount'}{$chr}{$pos};
				}else{
					$total{'usablecount'} = $h{'usablecount'}{$chr}{$pos};
				}
			}
                        if ($h{'variablecount'}{$chr}{$pos}){
                                if ($total{'variablecount'}){
                                        $total{'variablecount'} += $h{'variablecount'}{$chr}{$pos};
                                }else{
                                        $total{'variablecount'} = $h{'variablecount'}{$chr}{$pos};
                                }
                        }
		}elsif($pos>$e){
                      #then the last one is done and the next one starts
			if ($total{'variablecount'}){ #If there are variable sites in this window.
				foreach my $stat (@statlist){
					unless($total{$stat}){
						$total{$stat} = "0";
					}
				}
				unless($total{'usablecount'}){
					$total{'usablecount'} = "0";
				}
               	        	unless($total{'variablecount'}){
                                	$total{'variablecount'} = "0";
                        	}
				
				if (($total{"FstNum"}) and ($total{"FstDenom"})){
					$total{"Fst"} = ($total{"FstNum"}/$total{"FstDenom"});
				}else{
					$total{"Fst"} = "0";
				}
				my %ave;
				foreach my $stat (@statlist){
					if ($total{$stat}){
						$ave{$stat} = ($total{$stat} / $total{"usablecount"});	
					}else{
						$ave{$stat} = "0";
					}
				}
				print "$chr\t$s\t$e\t".($s+($Size/2))."\t$total{'usablecount'}\t$total{'variablecount'}\t$ave{'Hexp'}\n";
			}else{ #If there are usable sites but no variable sites
				unless ($total{'usablecount'}){
					$total{'usablecount'} = "0";
				}
				print "$chr\t$s\t$e\t".($s+($Size/2))."\t$total{'usablecount'}\t0\t0\n";
			}
			# reset these
			$count = 1;
 			foreach my $stat (@statlist){
				if ($h{$stat}{$chr}{$pos}){
					if ($h{$stat}{$chr}{$pos} ne "NA"){
						$total{$stat} = $h{$stat}{$chr}{$pos};
					}else{
						$total{$stat} = 0;
					}
				}else{
					$total{$stat} = 0;
				}
			}
			if ($h{'usablecount'}{$chr}{$pos}){
				$total{'usablecount'} = $h{'usablecount'}{$chr}{$pos};
			}else{
				$total{'usablecount'} = 0;
			}
                        if ($h{'variablecount'}{$chr}{$pos}){
                                $total{'variablecount'} = $h{'variablecount'}{$chr}{$pos};
                        }else{
                                $total{'variablecount'} = 0;
                        }          
                        until ($pos<$e){
                                $s += $Size;
                                $e += $Size;
				if ($pos>$e){
					print "$chr\t$s\t$e\t".($s+($Size/2))."\t0\t0\tNA\n";
				}
                        }    
		}
	}
	if ($total{'variablecount'}){
	        foreach my $stat (@statlist){
        		unless($total{$stat}){
                		$total{$stat} = "0";
               		}
        	}
        	unless($total{'usablecount'}){
        		$total{'usablecount'} = "0";
        	}
        	unless($total{'variablecount'}){
                	$total{'variablecount'} = "0";
        	}

        	if (($total{"FstNum"}) and ($total{"FstDenom"})){
        		$total{"Fst"} = ($total{"FstNum"}/$total{"FstDenom"});
        	}else{
        		$total{"Fst"} = "0";
        	}
        	my %ave;
        	foreach my $stat (@statlist){
        		if ($total{$stat}){
	                	$ave{$stat} = ($total{$stat} / $total{"usablecount"});
	              	}else{
        	        	$ave{$stat} = "0";
	                }
		}

		print "$chr\t$s\t$e\t".($s+($Size/2))."\t$total{'usablecount'}\t$total{'variablecount'}\t$ave{'Hexp'}\n";
	}else{
		print "$chr\t$s\t$e\t".($s+($Size/2))."\t$total{'usablecount'}\t0\t0\t0\n";
	}
}

