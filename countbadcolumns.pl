sub count_bad_columns{
	open IN, $_[0];
	while (<IN>){
		my @a = split (/\t/,$_);
		if ($. == 2){
			foreach my $i (0..$#a){
				if ($a[$i] =~ /^(AA|AT|AC|AG|TT|TA|TC|TG|CC|CA|CT|CG|GG|GC|GA|GT|NN)$/){
					return ("FALSE", "$i");
				}
				elsif ($a[$i] =~ /^(A|T|G|C|W|R|M|S|K|Y|N)$/){
					return ("TRUE", "$i");
				}
			}
		}
	}
	close IN;
}
1;					
