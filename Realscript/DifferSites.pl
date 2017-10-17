#!use#! /usr/bin/perl -w
use strict;
#use List::Util qw(first max maxstr min minstr reduce shuffle sum);

######################################################################################
# This script takes alignment sequence fasta file and export the differsites table
# Author: Lei Yang                                      <<<<<<<<< YANG >>>>>>>>>>>
# Date: 2010-05-24                                      <<<<<<<<< LEI  >>>>>>>>>>>
# Annotation: Lei Yang 2010-05-11
# LastModify: Lei Yang 2011-07-18
# Usage: perl Fasta2Phylip.pl -i inputFastaFile -o outputFile -d differNO. -e exportNO.
######################################################################################

#declare vars
my $help = "
		This script use the alignment fasta file
		Usage: perl DifferSites [options] [value]
		options	value
		-i	inputfile must be enter
		-o	outputfile must be enter
		-r	start-end|site[,start-end|site] to search
		-d	number of differs
		-e	export sites 
		-g	except gap '-'
		-hide	hide the same vs concens use [.]
		-f	use first seq as conseq
		-s	[startnumber] all the sites output will add it
use like this: perl DifferSites.pl -i aa714.fas -o aa714.sites -f -hide -d 2
		";

my ($infile,$outfile,$Range,$Differ,$unixFile,$Gap,$Hide,$First,$Start);

for(my $i=0;$i<=$#ARGV;)
{
	if($ARGV[$i] eq '-i'){ $infile = $ARGV[++$i];}
	if($ARGV[$i] eq '-o'){ $outfile = $ARGV[++$i];}
	if($ARGV[$i] eq '-r'){ $Range = $ARGV[++$i];}
	if($ARGV[$i] eq '-d'){ $Differ = $ARGV[++$i];}
	if($ARGV[$i] eq '-g'){ $Gap = 1;}
	if($ARGV[$i] eq '-hide'){ $Hide = '.';}
	if($ARGV[$i] eq '-f'){ $First = 1;}
	if($ARGV[$i] eq '-s'){ $Start = $ARGV[++$i];}
	$i++;
}

if(! $infile){ die $help;}
unless(defined $First){$First = 0;}
$unixFile = $infile.".unix";

open(IN,$infile) or die "Couldn't open $infile: $!\n";
open(OUT,">$outfile") or die "Couldn't open $outfile: $!\n";
close IN;
close OUT;

ConvertToUnix ($infile, $unixFile);	#yang# replace the \r or \r\n to \n
print "filein:" , $infile , "\n" if $infile;
print "fileout:", $outfile , "\n" if $outfile;
print "hide:" , $Hide , "\n" if $Hide;
print "first seq:" , $First , "\n" if defined $First;
print "differ than " , $Differ , "\n" if $Differ;
print "find within range:" , $Range , "\n" if $Range;
print "ignore gap?" ,$Gap , "\n" if defined $Gap;
Report($unixFile,$outfile,$Range,$Differ,$Gap,$Hide,$First,$Start);

unlink ($unixFile);
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

######################################################################################
sub Report{
	my ($unixFile, $outfile,$Range,$Differ,$Gap,$Hide,$First,$Start) = @_;
	my (@names,@seqs,@differsites);
	my $Consenus;

	my @Range = GetRanges($Range) if ($Range);
	my ($seqCount,$seqLen) = ReadFastaSeqs($unixFile,\@names,\@seqs);

	@Range = (1..$seqLen) unless(@Range);
	unless($Differ){ $Differ = 1;}
	#unless($Gap){$Gap = 1;}
	
	foreach my $site (@Range)
	{
		#unless($site){die "zero found in @Range"; }
		my %arpha;
		foreach (@seqs) {$arpha{substr($_,$site-1,1)} += 1;}
		my ($sum,$max,$maxkey) = (0,0,0);
		while ( my ($key, $value) = each %arpha )
		{
			unless($Gap && $key eq '-'){
  				$sum += $value;
				if ($max < $value)
				{
					$max = $value;
					$maxkey = $key;
				}
			}
		}
		push @differsites,$site if($sum - $max >= $Differ);
		$Consenus .= $maxkey;
	}
	$Consenus = $seqs[0] if ($First);
	print "Sites are found!\n";
	open(OUT,">$outfile") or die "Couldn't open $outfile: $!\n";
	print OUT "There are $seqCount sequnces in $unixFile file with $seqLen sites.\n";
	print OUT "$#differsites sites found differ more than $Differ seqs.\n";
	print OUT "SitesNumber\t" . TabArray($Start,\@differsites) . "\tN-Glyc\n";
	print OUT (($First) ? $names[0] : "Consenus") . "\t" . GetSites($Consenus,@differsites) . "\n";
	#2016-08-29 ADD  N-Glycosylation sites printing progress###
	for(my $i=$First; $i<$seqCount; $i++){
		if($Hide){print OUT $names[$i] . "\t" . GetHideSites($Consenus,$seqs[$i],@differsites);}
		else{print OUT $names[$i] . "\t" . GetSites($seqs[$i],@differsites);}
		$_ = $seqs[$i];
		$_ =~ s/-//g;
		my @N_Glyc = /N.[S|T]/g;
		print OUT "\t". (my $N_Glyc = @N_Glyc) . "\n";
	}
	
	#if($Hide){
	#	for(my $i=$First; $i<$seqCount; $i++){
	#		print OUT $names[$i] . "\t" . GetHideSites($Consenus,$seqs[$i],@differsites) . "\n";
	#	}
	#}else{
	#	for(my $i=$First; $i<$seqCount; $i++){
	#		print OUT $names[$i] . "\t" . GetSites($seqs[$i],@differsites) . "\n";
	#	}
	#}
	close OUT;
}
######################################################################################
sub GetRanges{

	my $Range = shift;
	my @Range; 
	foreach (split /,/,$Range)
	{
		if(/^(\d+)$/){push @Range,$1 if ($1>0);}
		elsif(/(\d+)-(\d+)/){ 
			if($1<=$2){push @Range,$1..$2 if ($1>0);}
			elsif($1>$2){ die "the range $1 is bigger than $2\n";}
		}
		else{ die "Range Option should like start-end|site[,start-end|site]\n";}		
	}
	@Range;
}
######################################################################################
sub ReadFastaSeqs{

	my ($unixFile,$names,$seqs) = @_;

	my $seqCount = my $seqLen = 0;
	my $seq = my $seqName = "";

	die "No array to store names\n" if (! $names);
	die "No array to store seqs\n" if (! $seqs);
	open(IN, $unixFile) || die "Can't open $unixFile\n";

	while(my $line = <IN>) {
		chomp $line;	
		next if($line =~ /^\s*$/);	#yang# pass the blank row
	
		if($line =~ /^>(\S+)/) {	#yang# find a seq begin and store the name to $1 using (\S+)
			if ($seqCount == 1) 
				{$seqLen = length $seq;}
			if ($seqCount) {	#yang# if it a new seq but the first one, we need to store the latest
				my $len = length $seq;
				if ($len == $seqLen) {
					push @{$names},$seqName;
					push @{$seqs},uc $seq;		#yang# export one Seq
					$seq = $seqName = "";
				}else {					#yang# find Seqlength Error
					unlink $unixFile;
					die "Error: the sequence length of $seqName is not same as others.\n";
				}
			}	
			$seqName = $1;
			$seqCount++;
		}else {
			$seq .= $line;		#yang# seq store in one line
		}		
	}
	close IN;
	# check the length of last sequence
	my $len = length $seq;
	if ($len == $seqLen) {
		push @{$names},$seqName;
		push @{$seqs},$seq;		#yang# export one Seq
		$seq = $seqName = "";
	}else {
		unlink $unixFile;
		die "Error: the sequence length of $seqName is not same as others.\n";
	}
	
	print "Reading finished!\n";
	($seqCount,$seqLen);

}
######################################################################################
sub GetSites{
	my $seq = shift;
	my $sites;
	foreach (@_){
		$sites .= substr($seq,$_-1,1);
		$sites .= "\t";
	}
	$sites;
}
######################################################################################
sub GetHideSites{
	my $con = shift;	
	my $seq = shift;
	my $sites;
	foreach (@_){
		if(substr($seq,$_-1,1) eq substr($con,$_-1,1)){$sites .= "\.";}
		else{$sites .= substr($seq,$_-1,1);}	
		$sites .= "\t";
	}
	$sites;
}
######################################################################################
sub TabArray{
	my $start = shift;
	my $sites = shift;
	my $s;
	foreach(@{$sites}){
		$s .= (($start) ? $_+$start : $_) . "\t";
	}
	$s;
}







