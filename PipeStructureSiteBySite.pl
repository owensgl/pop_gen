#!/usr/bin/perl
use warnings;
use strict;

my $in = $ARGV[0];
my $des = $ARGV[1];
my $name = $in;
$name =~ s/snp_tables\///;
my $table = $in;
my $home = "/home/owens/RNAseq";

foreach my $lg (1..17){
	#pull out that LG
	if ($lg < 10){
		$lg = "0".$lg;
	}
	my $tmp;
	my $tmp_name = $name."_tmp";
	my $tmp_header = $name."_header";
	my $numCols;
	my $snpTable = $table."_LG".$lg;
	my $dir = "linkage_structure/$name"."LG".$lg;
	my $structureFile = "$dir/$name".$lg.".structure";

	unless (-d "$dir") {
		#system ("mkdir $dir");
		system ("mkdir $dir");
	}

	my $lg_formatted = "HanXRQChr".$lg;
	$tmp = 'awk \'$1'." == \"".$lg_formatted."\"\' $table > $tmp_name\n";
	print $tmp."\n";
	system ("$tmp");
	# get line count
	my $lineCount = `wc -l $tmp_name`; 
	if($lineCount =~m/(\d+)/){
		$lineCount = "$1";
	}
	else{
		die;
	}
	#stick the header on
	#system "cat header str_tmp > $snpTable\n";
	system "head -n1 $table > $tmp_header\n";
	system "cat $tmp_header $tmp_name > $snpTable\n";
	#make a directory
	###system "mkdir $dir\n";
	system "cp $home/linkage_structure/extraparams $dir/extraparams\n";
        system "cp $home/linkage_structure/mainparams $dir/mainparams\n";
	
	#format for structure
	system "perl /home/owens/bin/reformat/SNPtable2structure_sitebysite.pl $snpTable $des $structureFile\n";
	#print "perl bin/format4structure_special.pl $snpTable $des $structureFile\n";

	#get the number of inds
	my $numInds;
	$numInds = `wc -l $des`;
	chomp($numInds);
	#run structure
	#print "cd /home/greg/wb/DomTech/$dir/\n";
	
	chdir("$dir") or die "$!";
	print "/home/owens/bin/console/structure -L $lineCount -N $numInds -i $home/$structureFile -o $home/$dir/results\n";
	system "/home/owens/bin/console/structure -L $lineCount -N $numInds -i $home/$structureFile -o $home/$dir/results";
	chdir("../../") or die "$!";
	system "rm $tmp_name";
	system "rm $tmp_header";
}


















