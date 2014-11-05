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
my @statlist = ("dxy","FstNum","FstDenom","Hexp1","Hexp2");
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
	}else{
		foreach my $i (0..$#a){
			if ($colname[$i] eq "CHROM"){
				$chr = $a[$i];
				push (@chrlist, $chr);
			}elsif ($colname[$i] eq "POS"){
				$pos = $a[$i];
				$poslist{$chr}{$pos}++;
			}elsif ($colname[$i] eq "Dxy"){
				$stat_hash{'dxy'} = $a[$i];
			}elsif ($colname[$i] eq "FstNum"){
				$stat_hash{'FstNum'} = $a[$i];
			}elsif ($colname[$i] eq "FstDenom"){
				$stat_hash{'FstDenom'} = $a[$i];
			}elsif ($colname[$i] eq "Hexp1"){
				$stat_hash{'Hexp1'} = $a[$i];
			}elsif ($colname[$i] eq "Hexp2"){
				$stat_hash{'Hexp2'} = $a[$i];
			}
		}
		foreach my $stat (@statlist){
			if ($stat_hash{$stat} ne "NA"){
				$h{$stat}{$chr}{$pos} = $stat_hash{$stat};
				$count{$stat}{$chr}{$pos}++;
				if ($stat eq "dxy"){
					$h{'usablecount'}{$chr}{$pos}++;
				}
			}
		}
	}
}	

my @uniq_chr = uniq(@chrlist);
print "Chr\tStartPos\tEndPos\tMidPos\tUsableSites\tDxy\tFst\tHexp1\tHexp2\n";
foreach my $chr (@uniq_chr){
	my $s = $Start;
	my $e = $Size;
	my %total;
	my $count;
	foreach my $pos (sort {$a<=>$b} keys %{$poslist{$chr}}){
		if (($pos >= $s)and($pos<$e)){
			$count++;
			foreach my $stat(@statlist){
				if ($h{$stat}{$chr}{$pos}){
					if ($total{$stat}){
						$total{$stat} += $h{$stat}{$chr}{$pos};
						$totalCount{$stat} += $count{$stat}{$chr}{$pos};
					}else{
						$total{$stat} = $h{$stat}{$chr}{$pos};
						$totalCount{$stat} = $count{$stat}{$chr}{$pos};
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
		}elsif($pos>$e){
                      #then the last one is done and the next one starts
			foreach my $stat (@statlist){
				unless($total{$stat}){
					$total{$stat} = "NA";
				}
				unless($totalCount{$stat}){
					$totalCount{$stat} = "0";
				}
			}
			unless($total{'usablecount'}){
				$total{'usablecount'} = "0";
			}
				
			if (($total{"FstNum"} ne "NA") and ($total{"FstDenom"} ne "NA")){
				$total{"Fst"} = ($total{"FstNum"}/$total{"FstDenom"});
			}else{
				$total{"Fst"} = "NA";
			}
			my %ave;
			foreach my $stat (@statlist){
				if (($totalCount{$stat}) and ($total{$stat} ne "NA")){
					$ave{$stat} = ($total{$stat} / $totalCount{$stat});	
				}else{
					$ave{$stat} = "NA";
				}
			}

				
			
			print "$chr\t$s\t$e\t".($s+($Size/2))."\t$totalCount{'dxy'}\t$ave{'dxy'}\t$total{'Fst'}\t$ave{'Hexp1'}\t$ave{'Hexp2'}\n";
			# reset these
			$count = 1;
 			foreach my $stat (@statlist){
				if ($h{$stat}{$chr}{$pos}){
					$total{$stat} = $h{$stat}{$chr}{$pos};
					$totalCount{$stat} = $count{$stat}{$chr}{$pos};
				}else{
					$total{$stat} = 0;
					$totalCount{$stat} = 0;
				}
			}
			if ($h{'usablecount'}{$chr}{$pos}){
				$total{'usablecount'} = $h{'usablecount'}{$chr}{$pos};
			}else{
				$total{'usablecount'} = 0;
			}          
                        until ($pos<$e){
                                $s += $Size;
                                $e += $Size;
				if ($pos>$e){
					print "$chr\t$s\t$e\t".($s+($Size/2))."\t0\tNA\tNA\tNA\tNA\n";
				}
                        }    
		}
	}
        foreach my $stat (@statlist){
        	unless($total{$stat}){
                	$total{$stat} = "NA";
                }
	        unless($totalCount{$stat}){
                	$totalCount{$stat} = "0";
                }
        }
        unless($total{'usablecount'}){
        	$total{'usablecount'} = "0";
        }

        if (($total{"FstNum"} ne "NA") and ($total{"FstDenom"} ne "NA")){
        	$total{"Fst"} = ($total{"FstNum"}/$total{"FstDenom"});
        }else{
        	$total{"Fst"} = "NA";
        }
        my %ave;
        foreach my $stat (@statlist){
        	if (($totalCount{$stat}) and ($total{$stat} ne "NA")){
                	$ave{$stat} = ($total{$stat} / $totalCount{$stat});
              	}else{
                	$ave{$stat} = "NA";
                }
	}

        if ($total{"Fst"} ne "NA"){
        	$ave{"Fst"} = $total{"Fst"}/$totalCount{"FstDenom"};
        }	
	else{
		$ave{"Fst"} = "NA";
	}

	print "$chr\t$s\t$e\t".($s+($Size/2))."\t$totalCount{'dxy'}\t$ave{'dxy'}\t$total{'Fst'}\t$ave{'Hexp1'}\t$ave{'Hexp2'}\n";
}

