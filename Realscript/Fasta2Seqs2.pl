#! /usr/bin/perl -w

######################################################################################
# This script takes alignment sequence fasta file and converts it to phylip file
# Author: Wenjie Deng
# Date: 2007-01-29
# Annotation: Lei Yang 2014-12-12
# Usage: perl Fasta2Seqs.pl inputFastaFile foldernamelist
######################################################################################
use strict;

my $usage = "Usage: perl Fasta2Seqs.pl inputFastaFile foldernamelist\n";
my $infile = shift or die($usage);	# input nexus file
my $foldernamelist = shift;
my $unixFile = "tmp.unix";

my @fasfiles;

if(-d $infile){
	opendir(DIRHANDER,$infile) || die "4 can't opendir $infile: $!";
	@fasfiles = grep { /.*\.(fa|fas|fasta)$/ && -f "$infile/$_" } readdir(DIRHANDER);	
	for (my $i=0;$i<=$#fasfiles;$i++)
	{$fasfiles[$i] = "$infile/$fasfiles[$i]";}
	closedir(DIRHANDER);	
}else{
	open (IN, $infile) or die "5 Couldn't open $infile: $!\n";
	close IN;
	@fasfiles = ($infile);
}

foreach my $file (@fasfiles){
	ConvertToUnix ($file, $unixFile);	#yang# replace the \r or \r\n to \n
	ChangetoSeqs($unixFile, $foldernamelist);
	unlink ($unixFile);
}
print "All done!\n";

exit 0;


######################################################################################
sub ConvertToUnix {
	my ($infile, $unixFile) = @_;
	open (IN, $infile) or die "Couldn't open $infile: $!\n";
	open (OUT, ">$unixFile") or die "Couldn't open $unixFile: $!\n";
	my @buffer = <IN>;
	close IN;
	my $line = "";			#yang# store the input Array to Scalary
	foreach my $element (@buffer) {
		$line .= $element;
	}
	if ($line =~ /\r\n/) {
		$line =~ s/\r//g;	#yang# replace \r\n to \n
	}elsif ($line =~ /\r/) {
		$line =~ s/\r/\n/g;	#yang# replace \r to \n
	}
	print OUT $line;	
	close OUT;	
}


###############################################################################################################
sub ChangetoSeqs {
	my ($unixFile,$foldernamelist) = @_;
	my $texon;
	my $isprint;
	#my @foldernamelist;
	open(IN, $foldernamelist) || die "Can't open $foldernamelist\n";
	while(<IN>) {
		chomp;
		mkdir $_;
	}
	close IN;
	#print @foldernamelist;
	
	
	open(IN, $unixFile) || die "Can't open $unixFile\n";
	#open(OUT, ">test");
	while(my $line = <IN>) {
		chomp $line;	
		next if($line =~ /^\s*$/);	
		if($line =~ /^>(\S+)/) {	#yang# find a seq begin and store the name to $1 using (\S+)
			$isprint = 0;
			$texon = substr($1,0,11);
			if (-e $texon && -d $texon){
				$isprint = 1;
				close OUT;
				open (OUT, ">$texon/$1.seq") || die "11 Couldn't open File $texon/$1.seq: $!\n";
				#print OUT "\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n";
				print OUT "\n\n\n\n\n\n\n\n";
			}
		}
		
		
		
		else {
			if($isprint){print OUT $line;}		#yang# seq store in one line
		}
	}
	close IN;
	close OUT;
	
}

