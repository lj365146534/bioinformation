#! /usr/bin/perl -w
use strict;

my $usage = "perl seperate_fasta.pl infile.fas\n";

my $infile = shift or die($usage);	# input nexus file
my $unixFile = $infile.".unix";
ConvertToUnix($infile, $unixFile);	#yang# replace the \r or \r\n to \n
my %fas;
my (@names,@seqs);
my ($seqCount,$seqLen,$isaligned) = ReadFastaSeqs($unixFile,\@names,\@seqs);
	
	###set fasta seqs to hash
	for(my $i=0; $i<$seqCount; $i++){
		$fas{$names[$i]} = $seqs[$i];
	}

my @sort_names = sort{my $da = get_tags($a);my $db = get_tags($b);$da cmp $db;} @names;
my $tag = "";
open MYOUT , ">temp";
foreach my $name (@sort_names)
{
	unless(get_tags($name) eq $tag){
		$tag = get_tags($name); 
		close MYOUT;
		open(MYOUT,">$tag.fas") or die "Couldn't open $tag: $!\n";
	}
	print MYOUT ">$name\n";
	print MYOUT $fas{$name};
	print MYOUT "\n";
}
close MYOUT;
unlink ($unixFile);
unlink ("temp");
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
sub ReadFastaSeqs{

	my ($unixFile,$names,$seqs) = @_;

	my ($seqCount,$seqLen,$isaligned,$len) = (0,0,1,0);
	my $seq = my $seqName = "";

	die "No array to store names\n" if (! $names);
	die "No array to store seqs\n" if (! $seqs);
	open(IN, $unixFile) || die "Can't open $unixFile\n";

	while(my $line = <IN>) {
		chomp $line;	
		next if($line =~ /^\s*$/);	#yang# pass the blank row
	
		if($line =~ /^>(\S+)/) {	#yang# find a seq begin and store the name to $1 using (\S+)
			if ($seqCount == 1)		#initial $seqCount as the first seq
				{$seqLen = length $seq;}
			if ($seqCount) {	#yang# if it a new seq but the first one, we need to store the latest
				$len = length $seq;
				push @{$names},$seqName;
				push @{$seqs},uc $seq;		#yang# export one Seq
				$seq = $seqName = "";
				if ($len != $seqLen){$isaligned=0;$seqLen = ($len > $seqLen) ? $len : $seqLen;}
			}
			$seqName = $1;
			$seqCount++;
		}else{
			$seq .= $line;		#yang# seq store in one line
		}
	}
	close IN;
	# check the length of last sequence
	$len = length $seq;
	if ($len != $seqLen){$isaligned=0;}
	push @{$names},$seqName;
	push @{$seqs},$seq;		#yang# export one Seq
	$seq = $seqName = "";
	
	print "Reading finished!\n";
	($seqCount,$seqLen,$isaligned);

}

######################################################################################
sub get_tags{
	my $da = shift;
	$da =~ s/.+[-_]//;
	return $da;
}