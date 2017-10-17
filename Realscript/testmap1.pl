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

=pod	system ("bowtie2 --very-sensitive-local -x /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -1 $outdir/$basename/trim_1.fastq -2 $outdir/$basename/trim_2.fastq  --un-conc $outdir/$basename/second_map/unmap.fq -S $outdir/$basename/second_map/second.sam -p $threads");
	#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
	#	-x <bt2-idx> 由bowtie2-build所生成的索引文件的前缀。首先 在当前目录搜寻，然后在环境变量 BOWTIE2_INDEXES 中制定的文件夹中搜寻。
	#-1 <m1> 双末端测寻对应的文件1。可以为多个文件，并用逗号分开；多个文件必须和 -2 <m2> 中制定的文件一一对应。比如:"-1 flyA_1.fq,flyB_1.fq -2 flyA_2.fq,flyB_2.fq". 测序文件中的reads的长度可以不一样。-2 <m2> 双末端测寻对应的文件2.
	#-U <r> 非双末端测寻对应的文件。可以为多个文件，并用逗号分开。测序文件中的reads的长度可以不一样。
	#-S <hit> 所生成的SAM格式的文件前缀。默认是输入到标准输出。
	#--very-sensitive-local Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
	#--un-conc <path> 将不能和谐比对的paired-end reads写入<path>
	system ("samtools view -bS $outdir/$basename/second_map/second.sam > $outdir/$basename/second_map/$basename._second.bam");
	#samtool view -bS命令把sam文件变成bam文件
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/second_map/$basename._second.bam -o $outdir/$basename/second_map/$basename._second_sort.bam");
	#sort排序（不明白有啥用）
	system ("samtools index $outdir/$basename/second_map/$basename._second_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	print "\n\n\n";
	system ("samtools mpileup -Q 0 $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/second_pileup");
	#samtools mpileup -Q [INT] Minimum base quality for a base to be considered [13]
	#-f FILE  The faidx-indexed reference file in the FASTA format. The file can be optionally compressed by razip. [null] 
    #-l FILE  BED or position list file containing a list of regions or sites where pileup or BCF should be generated [null] 

	#system ("copy /home/lei/Desktop/flu_toolkit/zhuwenfei/ASichuan12009.fasta /home/lei/Desktop/flu_toolkit/$outdir/$basename/second_map/");

	#system ("bowtie2-build $outdir/$basename/second_map/reference.fa $outdir/$basename/second_map/reference.fa");

	system ("samtools mpileup -Q 0 -f /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/ref_dep");
	system ("perl $outdir/../draw_depth.pl $outdir/$basename/second_map/ref_dep > $outdir/$basename/second_map/coverage");
	system ("Rscript $outdir/../DR.r $outdir/$basename/");
=cut
	system ("perl $outdir/../pileup2fasta.pl -d 3 -S 0.3 $outdir/$basename/second_map/second_pileup >  $outdir/$basename/$basename.30.fa");
	#	system ("samtools mpileup -vuf  $outdir/$basename/second_map/first_consensus.fa  $outdir/$basename/second_map/$basename._second_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/consensus.fa");
#	system ("java -jar $outdir/../FluAnn.jar $outdir/$basename/$basename.fa $outdir/$basename/$basename.html");
}	
	print "DONE!\nPREES ENTER!"
