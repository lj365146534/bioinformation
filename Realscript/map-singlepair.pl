#! /usr/bin/perl -w
use strict;
use grep_flu_fregment;
my $outdir = "/home/lei/Desktop/flu_toolkit/flu_out";
my $indir  = "/home/lei/Desktop/flu_toolkit/flu_in";
my $threads = 10;
chdir $indir;
my @filelist = split (/\s/ , `ls`);
print @filelist;
#print $filelist[1];
foreach my $eachfile (@filelist){
	chomp $eachfile;
	my $basename = $eachfile;  
	$basename =~ s/\Q.fastq//;
	$basename =~ s/\Q.fq//;
	#print $basename;
	mkdir "$outdir/$basename";
#	system ("fastqc -o $outdir/$basename $indir/$eachfile");
=trim
	open (BE,"<$indir/$eachfile") or die;
	open (AF,">$outdir/$basename/trim.fastq") or  die;
	while(<BE>){
		my $a = <BE>;
		my $b = <BE>;
		my $c = <BE>;
		print AF "$_$a$b$c" if length $a > 100;	
	}
	close AF;
	close BE;
=cut	
	
=pod	system ("fastx_trimmer -Q33 -t 20 -m 100 -i $indir/$eachfile -o $outdir/$basename/trim.fastq");
	system ("fastqc -o $outdir/$basename $outdir/$basename/trim.fastq");
#	system ("tail -n 30000 $indir/$eachfile > $outdir/$basename/cat_file.fa");
	system ("head -50000 $outdir/$basename/trim.fastq > $outdir/$basename/cat_file.fq");
#	system ("velveth $outdir/$basename/velvet 33 -short -fastq $outdir/$basename/cat_file.fa");
#	system ("velvetg $outdir/$basename/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes   ");
#	system ("velvetg $outdir/$basename/velvet  ");
	system("runAssembly -cpu 10  -o $outdir/$basename/velvet $outdir/$basename/cat_file.fq");
	system ("blastn -query $outdir/$basename/velvet/454AllContigs.fna -out $outdir/$basename/velvet/contigs.blast -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10");

###
	system ("tail -n 100000 $outdir/$basename/trim.fastq > $outdir/$basename/cat_file2.fa");
	system ("velveth $outdir/$basename/velvet 33 -short -fastq $outdir/$basename/cat_file2.fa");
	system ("velvetg $outdir/$basename/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes ");
#	system ("velvetg $outdir/$basename/velvet  ");
	system ("blastn -query $outdir/$basename/velvet/contigs.fa  -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10 >> $outdir/$basename/velvet/contigs.blast");
###


	my $aa = "$outdir/$basename/velvet";
	grep_flu_fregment::grep_flu ( $aa, "contigs.blast");
	mkdir  "$outdir/$basename/first_map";
	system ("blastdbcmd -db $outdir/../fasta.fa -entry_batch $aa/blast_list > $outdir/$basename/blast.fa");
	system ("fasta_formatter -i $outdir/$basename/blast.fa -o $outdir/$basename/first_map/blast_w.fa -w 0");
	
	open (IN,"<$outdir/$basename/first_map/blast_w.fa") or die;
	open (OUT,">$outdir/$basename/first_map/$basename._first_ref.fa") or die;
	while (my $line = <IN>){
		if($line =~ /^>/){
			$line =~ s/.*\(/>/;
			$line =~ s/\Q) /_/;			
		}
		else {
			chomp $line;
			$line = $line."\n";
		}
		print OUT $line;
	}
	close IN;
	close OUT;

#	
	system ("bowtie2-build $outdir/$basename/first_map/$basename._first_ref.fa $outdir/$basename/first_map/$basename._first_ref.fa");
	system ("bowtie2 --very-sensitive-local  -x $outdir/$basename/first_map/$basename._first_ref.fa -q $outdir/$basename/trim.fastq -S $outdir/$basename/first_map/$basename._first.sam -p 10");
	system ("samtools view -bS $outdir/$basename/first_map/$basename._first.sam > $outdir/$basename/first_map/$basename._first.bam");
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/first_map/$basename._first.bam -o $outdir/$basename/first_map/$basename._first_sort.bam");
	system ("samtools index $outdir/$basename/first_map/$basename._first_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	print "\n\n\n";
	mkdir  "$outdir/$basename/second_map"; 
	system ("samtools mpileup -Q 0 $outdir/$basename/first_map/$basename._first_sort.bam > $outdir/$basename/first_map/bbb");
	system ("perl $outdir/../pileup2fasta.pl -S 0 -d 1 $outdir/$basename/first_map/bbb > $outdir/$basename/second_map/first_consensus.fa");
=cut	#system ("samtools mpileup -vuf  $outdir/$basename/first_map/$basename._first_ref.fa  $outdir/$basename/first_map/$basename._first_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/second_map/first_consensus.fa");
	system ("bowtie2-build $outdir/$basename/second_map/first_consensus.fa $outdir/$basename/second_map/first_consensus.fa");
	system ("bowtie2 --very-sensitive-local -x $outdir/$basename/second_map/first_consensus.fa -q $outdir/$basename/trim.fastq --un $outdir/$basename/second_map/unmap.fq -S $outdir/$basename/second_map/second.sam -p 10");
	system ("samtools view -bS $outdir/$basename/second_map/second.sam > $outdir/$basename/second_map/$basename._second.bam");
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/second_map/$basename._second.bam -o $outdir/$basename/second_map/$basename._second_sort.bam");
	system ("samtools index $outdir/$basename/second_map/$basename._second_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	print "\n\n\n";
	system ("samtools mpileup -Q 0 $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/second_pileup");
	system ("samtools mpileup -Q 0 -f $outdir/$basename/second_map/first_consensus.fa $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/ref_dep");
	system ("perl $outdir/../draw_depth.pl $outdir/$basename/second_map/ref_dep > $outdir/$basename/second_map/coverage");
	system ("Rscript $outdir/../DR.r $outdir/$basename/");
	system ("perl $outdir/../pileup2fasta.pl -S 0.3 -d 3 $outdir/$basename/second_map/second_pileup >  $outdir/$basename/$basename.fa");
	#	system ("samtools mpileup -vuf  $outdir/$basename/second_map/first_consensus.fa  $outdir/$basename/second_map/$basename._second_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/consensus.fa");
#	system ("java -jar $outdir/../FluAnn.jar $outdir/$basename/$basename.fa $outdir/$basename/$basename.html");
	# (fastq_to_fasta -Q33 -i te.fq -o te.fa) also works;
	open (IN,"<$outdir/$basename/second_map/unmap.fq") or die;
	open (OUT,">$outdir/$basename/second_map/unmap.fa") or die;
	my $i = 1;
	while(<IN>){
		print OUT ">$_\n" if ($i%4 ==1);
		print OUT "$_\n"  if ($i%4 ==2);
		$i++;
	}
	close IN;
	close OUT;
	#
	#system ("blastn -query $outdir/$basename/second_map/unmap.fa -out $outdir/$basename/unmap.blast -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10");
=edit
	system ("cp $indir/$eachfile $outdir/$basename/second_map");
	system ("cp $outdir/$basename/second_map/first_consensus.fa $outdir/$basename");
	system ("cp $outdir/$basename/second_map/second.sam $outdir/$basename");
	
	system ("rm -rf $outdir/$basename/first_map");
	system ("rm -rf  $outdir/$basename/velvet");
	system ("rm $outdir/$basename/blast.fa $outdir/$basename/cat_file.fa $outdir/$basename/trim_fastqc.zip $outdir/$basename/trim.fastq");
=cut
}
