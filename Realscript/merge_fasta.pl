#! /usr/bin/perl -w
use strict;
use constant DE_BUG => 1;

my $file_extension = "fa|fas|fasta";
my $usage = "perl merge_fasta.pl \n";
my $infile = ".";
my @fasfiles;

opendir(DIRHANDER,$infile) || die "4 can't opendir $infile: $!";
#@fasfiles = grep { /.*\.fas(?!.*\..*)/ && -f "$infile/$_" } readdir(DIRHANDER);
@fasfiles = grep { /.*\.($file_extension)$/ && -f "$infile/$_" } readdir(DIRHANDER);
closedir(DIRHANDER);

print join("\n",@fasfiles) if DE_BUG;

open OUT, ">temp.fas";
foreach my $file (@fasfiles){
	open (IN, $file) or die "6 Couldn't open $file: $!\n";
	$file =~ /(.+)\.($file_extension)$/;
	my $name = $1;
	while (my $line = <IN>){
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			$line =~ s/_<.*>//gi;
			$line =~ s/>//gi;
			$line = join "",">",$name,"-",$line;
		}
		print OUT "$line\n";
	}
	close IN;
}
print "All done!\n";
close OUT;

exit 0;
