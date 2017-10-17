#! /usr/bin/perl -w

######################################################################################
# This script takes alignment sequence fasta file and converts it to phylip file
# Author: Lei Yang
# Date: 2013-10-29
# Annotation: 
# Usage: perl newick_tree_nodes.pl inputFastaFile outputPhilipFile
######################################################################################
use strict;

my $usage = "perl newick_tree_nodes.pl inputtextFile namelistFile\n";
my $infile = shift or die($usage);	# input nexus file
my $nameFile = shift or die($usage);	# output phylip file

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
	Replacename($file,$nameFile);
}


print "All done!\n";

exit 0;


######################################################################################

sub Replacename {
	my ($infile, $nameFile) = @_;
	open (IN, $infile) or die "Couldn't open $infile: $!\n";
	open (NAME_IN, $nameFile) or die "Couldn't open $nameFile: $!\n";
	my @buffer = <NAME_IN>;
	close NAME_IN;
	my $lines;
	while (<IN>){
		$lines .= $_;
	}
	close IN;
	my ($key,$value);
	foreach my $line (@buffer){
		chomp $line;
		$line =~ s/\r//g;
		next if($line =~ /^\s*$/);
		($key,$value) = split /\t/,$line;
			#print $key,$value;
			#while(<STDIN>){last;}
		$lines =~ s/\Q$key/$value/gi;
	}
	open (OUT,">$infile") or die "Couldn't open $infile: $!\n";
	print OUT $lines;
}
	
