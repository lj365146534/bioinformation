#! /usr/bin/perl -w

######################################################################################
# This script takes a subset from fastafiles
# Author: Lei Yang
# Date: 2013-10-8
# Annotation: 
# Usage: perl sub_dataset.pl inputFastaFile inputNamelist
######################################################################################
use strict;
use constant DE_BUG => 0;

my $usage = "Usage: perl sub_dataset.pl inputFastaFile inputNamelist matchpatten matchoperate\n
			matchpatten : eq or lk\tequal or like\n
			matchoperate : in or ex\tinclude or exclude\n";

my $infile = shift or die("1" . $usage);	# input fasta file or a directory 
					     # containing fasta files
my $namelistfile = shift or die("2" . $usage);	# a namelist seperate by line
my ($m,$matchpatten,$matchoperate) = (1,undef,undef);
$matchpatten = shift or $matchpatten = "eq";
$matchoperate = shift or $m = 1;

if($matchoperate eq 'in'){$m = 1;}
if($matchoperate eq 'ex'){$m = 0;}

open (IN, $namelistfile) or die "3. Couldn't open $namelistfile: $!\n";
my @namelist;
while(<IN>){
	chomp;
	push @namelist,$_;
}

#my @namelist = <IN>;
close IN;
print "@namelist";

my @fasfiles;

if(-d $infile){
	opendir(DIRHANDER,$infile) || die "4 can't opendir $infile: $!";
	@fasfiles = grep { /.*\.fas.*/ && -f "$infile/$_" } readdir(DIRHANDER);	
	for (my $i=0;$i<=$#fasfiles;$i++)
	{$fasfiles[$i] = "$infile/$fasfiles[$i]";}
	closedir(DIRHANDER);	
}else{
	open (IN, $infile) or die "5 Couldn't open $infile: $!\n";
	close IN;
	@fasfiles = ($infile);
}

foreach my $file (@fasfiles){
	#print "$file\n";
	subset("$file",$matchpatten,$m);
}

print "All done!\n";

exit 0;

sub subset{
	my ($filein,$matchpatten,$matchoperate) = @_;
	my ($path,$filename);
	open (IN, $filein) or die "6 Couldn't open $filein: $!\n";
	$filein =~ /(.*)\/(.+?)\.fas/;
	unless($1){$filein =~ /(.*)\.fas/;$filename=$1;}
	$path = ($1) ? "$1\/" : "";
	$filename = ($filename) ? $filename : $2;
	unless (-e $path . "subdata"){mkdir $path . "subdata";}
	open (OUT, ">$path" . "subdata\/$filename.fas") || die "7 Cant open $path" . "subdata\/$filename.fas: $!\n";
	my $is_print = 0;
	while (my $line = <IN>){
		chomp $line;
		next if $line =~ /^\s*$/;	#yang# pass the blank row
		if ($line =~ /^>(\S+)/) {	#yang# find a seq begin with '>'
			my $find = 0;
			my $text = $1;
			foreach(@namelist){  #find name?
				chomp;
				if(DE_BUG){print "$_\t$text\n";}
				if((($matchpatten eq "eq") && (lc($text) eq lc($_))) || (($matchpatten eq "lk") && ($text =~ /\Q$_/gi))){$find = 1;}
			}
			if($find == 1){$is_print = 1;}else{$is_print = 0;}
			#set print label
		}
		if ($is_print == $matchoperate){print OUT "$line\n";}
	}

	close IN;
	close OUT;
}
