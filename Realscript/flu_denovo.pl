#! /usr/bin/perl -w
use strict;
use warnings;
use grep_flu_fregment;

my $outdir = "/home/lei/Desktop/flu_toolkit/flu_out";
my $indir  = "/home/lei/Desktop/flu_toolkit/flu_in";
my $threads = 10;

chdir $indir;
my @filelist = split (/\s/ , `ls`);
print "@filelist\n";
#print $filelist[1];
foreach my $eachfile (@filelist){
	chomp $eachfile;
	my $basename = $eachfile;  
	$basename =~ s/\Q.fastq//;
	$basename =~ s/\Q.fq//;
	if ($basename=~m/_R1$/){
		$basename=~s/_R1$//;
	}else{
		next;
	}
	my $base_1 = "$basename"."_R1";
	my $base_2 = "$basename"."_R2";
	print "\n$indir/$base_2.fastq\n";

	die 1 unless (-e "$indir/$base_1.fastq");
	die 2 unless (-e "$indir/$base_2.fastq");
	#print $basename;
	mkdir "$outdir/$basename";



	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq $outdir/$basename/trim_1.fastq $outdir/$basename/unpair_1.fastq $outdir/$basename/trim_2.fastq $outdir/$basename/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	system ("fastqc -o $outdir/$basename $outdir/$basename/trim_1.fastq $outdir/$basename/trim_2.fastq");
	system ("perl /home/lei/Desktop/Chen/velvet/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl $outdir/$basename/trim_1.fastq $outdir/$basename/trim_2.fastq $outdir/$basename/trim_reads.fastq");
	system ("velveth $outdir/$basename/velvet 33 -short -fastq $outdir/$basename/trim_reads.fastq");
	system ("velvetg $outdir/$basename/velvet -max_branch_length 100 -max_divergence 0.2 -max_gap_count 3 -cov_cutoff 55  -scaffolding yes ");#原来为(250,0.4,10)

}

print "PRESS ENTER";
