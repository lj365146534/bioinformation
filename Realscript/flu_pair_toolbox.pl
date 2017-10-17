#! /usr/bin/perl -w
use strict;
use warnings;
use grep_flu_fregment;

my $outdir = "/home/lei/Desktop/storage/flu_toolkit/flu_out";
my $indir  = "/home/lei/Desktop/storage/flu_toolkit/flu_in";
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
#	system ("fastqc -o $outdir/$basename $indir/$eachfile");
	
#	system ("fastx_trimmer -Q33 -t 3  -i $indir/$base_1.fastq -o $outdir/$basename/trim_1.fastq");
#	system ("fastx_trimmer -Q33 -t 3  -i $indir/$base_2.fastq -o $outdir/$basename/trim_2.fastq");
	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq $outdir/$basename/trim_1.fastq $outdir/$basename/unpair_1.fastq $outdir/$basename/trim_2.fastq $outdir/$basename/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	
	system ("fastqc -o $outdir/$basename $outdir/$basename/trim_1.fastq $outdir/$basename/trim_2.fastq");
#	system ("tail -n 30000 $indir/$eachfile > $outdir/$basename/cat_file.fa");
	
	system ("head -25000 $outdir/$basename/trim_1.fastq > $outdir/$basename/cat_file.fq");
#	system ("velveth $outdir/$basename/velvet 33 -short -fastq $outdir/$basename/cat_file.fa");
#	system ("velvetg $outdir/$basename/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes   ");
#	system ("velvetg $outdir/$basename/velvet  ");
	system("runAssembly -cpu $threads  -o $outdir/$basename/velvet $outdir/$basename/cat_file.fq");
	system ("blastn -query $outdir/$basename/velvet/454AllContigs.fna -out $outdir/$basename/velvet/contigs.blast -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads $threads");

###
	system ("tail -n 50000 $outdir/$basename/trim_2.fastq > $outdir/$basename/cat_file2.fa");
	system ("velveth $outdir/$basename/velvet 33 -short -fastq $outdir/$basename/cat_file2.fa");
	system ("velvetg $outdir/$basename/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes ");
#	system ("velvetg $outdir/$basename/velvet  ");
	system ("blastn -query $outdir/$basename/velvet/contigs.fa  -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads $threads >> $outdir/$basename/velvet/contigs.blast");
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
#		else {
#			chomp $line;
#			$line = "NNNNNNNNNNNNN".$line."NNNNNNNNNNNNN"."\n";
#		}
		print OUT $line;
	}
	close IN;
	close OUT;

#	
	system ("bowtie2-build $outdir/$basename/first_map/$basename._first_ref.fa $outdir/$basename/first_map/$basename._first_ref.fa");
	system ("bowtie2 --very-sensitive-local  -x $outdir/$basename/first_map/$basename._first_ref.fa -1 $outdir/$basename/trim_1.fastq -2 $outdir/$basename/trim_2.fastq -S $outdir/$basename/first_map/$basename._first.sam -p $threads");
	system ("samtools view -bS $outdir/$basename/first_map/$basename._first.sam > $outdir/$basename/first_map/$basename._first.bam");#samtool view命令把sam文件转换为bam文件
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/first_map/$basename._first.bam -o $outdir/$basename/first_map/$basename._first_sort.bam");#samtool sort对bam排序，此处应该out为prefix文件
	system ("samtools index $outdir/$basename/first_map/$basename._first_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai；bai用于快速检索reads，对bam文件进行默认排序后才能进行index
	print "\n\n\n";
	mkdir  "$outdir/$basename/second_map"; 
	system ("samtools mpileup -Q 0 $outdir/$basename/first_map/$basename._first_sort.bam > $outdir/$basename/first_map/bbb");
	system ("perl $outdir/../pileup2fasta.pl -S 0 -d 1 $outdir/$basename/first_map/bbb > $outdir/$basename/second_map/first_consensus.fa");
	#system ("samtools mpileup -vuf  $outdir/$basename/first_map/$basename._first_ref.fa  $outdir/$basename/first_map/$basename._first_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/second_map/first_consensus.fa");
	system ("bowtie2-build $outdir/$basename/second_map/first_consensus.fa $outdir/$basename/second_map/first_consensus.fa");
	system ("bowtie2 --very-sensitive-local -x $outdir/$basename/second_map/first_consensus.fa -1 $outdir/$basename/trim_1.fastq -2 $outdir/$basename/trim_2.fastq  --un-conc $outdir/$basename/second_map/unmap.fq -S $outdir/$basename/second_map/second.sam -p $threads");
	system ("samtools view -bS $outdir/$basename/second_map/second.sam > $outdir/$basename/second_map/$basename._second.bam");
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/second_map/$basename._second.bam -o $outdir/$basename/second_map/$basename._second_sort.bam");
	system ("samtools index $outdir/$basename/second_map/$basename._second_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	print "\n\n\n";
	system ("samtools mpileup -Q 0 $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/second_pileup");
	system ("samtools mpileup -Q 0 -f $outdir/$basename/second_map/first_consensus.fa $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/ref_dep");
	system ("perl $outdir/../draw_depth.pl $outdir/$basename/second_map/ref_dep > $outdir/$basename/second_map/coverage");
	system ("Rscript $outdir/../DR.r $outdir/$basename/");
	system ("perl $outdir/../pileup2fasta.pl -d 3 -S 0.3 $outdir/$basename/second_map/second_pileup >  $outdir/$basename/$basename.fa");
	#	system ("samtools mpileup -vuf  $outdir/$basename/second_map/first_consensus.fa  $outdir/$basename/second_map/$basename._second_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/consensus.fa");
	system ("java -jar $outdir/../FluAnn.jar $outdir/$basename/$basename.fa $outdir/$basename/$basename.html");
	# (fastq_to_fasta -Q33 -i te.fq -o te.fa) also works;
=ann	
	open (IN,"<$outdir/$basename/second_map/unmap.1.fq") or die;
	open (OUT,">$outdir/$basename/second_map/unmap.fa") or die;
	my $i = 1;
	while(<IN>){
		print OUT ">$_\n" if ($i%4 ==1);
		print OUT "$_\n"  if ($i%4 ==2);
		$i++;
	}
	close IN;
	
	open (IN,"<$outdir/$basename/second_map/unmap.2.fq") or die;
	 $i = 1;
	while(<IN>){
		print OUT ">$_\n" if ($i%4 ==1);
		print OUT "$_\n"  if ($i%4 ==2);
		$i++;
	}
	close IN;
	close OUT;
	#
	system ("blastn -query $outdir/$basename/second_map/unmap.fa -out $outdir/$basename/unmap.blast -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads $threads");
	
	
=edit
	system ("cp $indir/$eachfile $outdir/$basename/second_map");
	system ("cp $outdir/$basename/second_map/first_consensus.fa $outdir/$basename");
	system ("cp $outdir/$basename/second_map/second.sam $outdir/$basename");
	
	system ("rm -rf $outdir/$basename/first_map");
	system ("rm -rf  $outdir/$basename/velvet");
	system ("rm $outdir/$basename/blast.fa $outdir/$basename/cat_file.fa $outdir/$basename/trim_fastqc.zip $outdir/$basename/trim.fastq");
=cut
}

print "PRESS ENTER";
