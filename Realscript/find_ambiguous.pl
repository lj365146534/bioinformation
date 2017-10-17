#! /usr/bin/perl -w

######################################################################################
# This script find ambiguous sites from fastafiles for each seq
# Author: Lei Yang
# Date: 2015-8-6
# Annotation: 
# Usage: perl find_ambiguous.pl inputFastaFile > out_csv
######################################################################################

use strict;
use constant DE_BUG => 0;

my $usage = "Usage: find_ambiguous.pl inputFastaFile > out.csv\n";

my $infile = shift or die("1" . $usage);	# input fasta file or a directory 
					     # containing fasta files
my @fasfiles;
my ($path,$filename)=("","");
#if(-d xxxx)  xxxx为目录
if(-d $infile){
	opendir(DIRHANDER,$infile) || die "4 can't opendir $infile: $!";
	@fasfiles = grep { /.*\.fas.*/ && -f "$infile/$_" } readdir(DIRHANDER);	
	closedir(DIRHANDER);
	$path = "$infile\/";
}else{
	open (IN, $infile) or die "5 Couldn't open $infile: $!\n";
	close IN;
	$infile =~ /(.*)\/(.*)/;
	if($1){$path = "$1\/";$filename=$2;}
	@fasfiles = ($infile);	
}

unless (-e $path . "subdata"){mkdir $path . "subdata";}


foreach my $file (@fasfiles){
	#print "$file\n";
	find($file);
}

print "All done!\n";

exit 0;

#print ambiguous("actgrglstea");

sub ambiguous{
	my ($seq) = @_;
	$seq =~ s/([actg]|-)//ig;#g全局，i忽略大小写，s为替换actg为空
	return $seq;
}

sub find{
	my $filein = shift;
	open (IN, $filein) or die "6 Couldn't open $filein: $!\n";
	if($filename){$filein=$filename;}
	open (OUT, ">$path" . "subdata\/$filein") || die "7 Cant open $path" . "subdata\/$filein: $!\n";
	
	while (my $line = <IN>){
		chomp $line;
		next if $line =~ /^\s*$/;	#yang# pass the blank row
		if($line =~ /^>(\S+)/) {	#yang# find a seq begin with '>'
			print OUT "\n$line\t";			
		}else{
			$line =ambiguous($line);
			print OUT "$line";
		}
	}
	close IN;
	close OUT;
}
