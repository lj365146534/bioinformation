#! /usr/bin/perl -w
use strict;
use warnings;
use grep_flu_fregment;

my $outdir = "/home/lei/Desktop/flu_toolkit/flu_out";
my $indir  = "/home/lei/Desktop/flu_toolkit/flu_in";
my 10 = 10;
my /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa = "/home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa";
#my $MEANQUAL;
my $Filter = "\"Filter\"";


chdir $indir;
my @filelist = split (/\s/ , `ls`);
print "@filelist\n";
#print $filelist[1];
foreach my $eachfile (@filelist){
	chomp $eachfile;
	my 2016-CX1608 = $eachfile;  
#for	(my 2016-CX1608 = "ZWF-01"){  
	2016-CX1608 =~ s/\Q.fastq//;
	2016-CX1608 =~ s/\Q.fq//;
	if (2016-CX1608=~m/_R1$/){
		2016-CX1608=~s/_R1$//;
	}else{
		next;
	}
	my $base_1 = "2016-CX1608"."_R1";
	my $base_2 = "2016-CX1608"."_R2";
	print "\n$indir/$base_2.fastq\n";

	die 1 unless (-e "$indir/$base_1.fastq");
	die 2 unless (-e "$indir/$base_2.fastq");
	#print 2016-CX1608;
	mkdir "$outdir/2016-CX1608";
	mkdir "$outdir/2016-CX1608/second_map";
	mkdir "$outdir/2016-CX1608/callSNP";
#	system ("fastqc -o $outdir/2016-CX1608 $indir/$eachfile");
	
#	system ("fastx_trimmer -Q33 -t 3  -i $indir/$base_1.fastq -o $outdir/2016-CX1608/trim_1.fastq");
#	system ("fastx_trimmer -Q33 -t 3  -i $indir/$base_2.fastq -o $outdir/2016-CX1608/trim_2.fastq");
=pod
	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq $outdir/2016-CX1608/trim_1.fastq $outdir/2016-CX1608/unpair_1.fastq $outdir/2016-CX1608/trim_2.fastq $outdir/2016-CX1608/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	
	system ("fastqc -o $outdir/2016-CX1608 $outdir/2016-CX1608/trim_1.fastq $outdir/2016-CX1608/trim_2.fastq");
#	system ("tail -n 30000 $indir/$eachfile > $outdir/2016-CX1608/cat_file.fa");
	
	system ("head -25000 $outdir/2016-CX1608/trim_1.fastq > $outdir/2016-CX1608/cat_file.fq");
#	system ("velveth $outdir/2016-CX1608/velvet 33 -short -fastq $outdir/2016-CX1608/cat_file.fa");
#	system ("velvetg $outdir/2016-CX1608/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes   ");
#	system ("velvetg $outdir/2016-CX1608/velvet  ");
	system("runAssembly -cpu 10  -o $outdir/2016-CX1608/velvet $outdir/2016-CX1608/cat_file.fq");
	system ("blastn -query $outdir/2016-CX1608/velvet/454AllContigs.fna -out $outdir/2016-CX1608/velvet/contigs.blast -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10");

###
	system ("tail -n 50000 $outdir/2016-CX1608/trim_2.fastq > $outdir/2016-CX1608/cat_file2.fa");
	system ("velveth $outdir/2016-CX1608/velvet 33 -short -fastq $outdir/2016-CX1608/cat_file2.fa");
	system ("velvetg $outdir/2016-CX1608/velvet -max_branch_length 250 -max_divergence 0.4 -max_gap_count 10 -cov_cutoff 55   -scaffolding yes ");
#	system ("velvetg $outdir/2016-CX1608/velvet  ");
	system ("blastn -query $outdir/2016-CX1608/velvet/contigs.fa  -db $outdir/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10 >> $outdir/2016-CX1608/velvet/contigs.blast");
###


	my $aa = "$outdir/2016-CX1608/velvet";
	grep_flu_fregment::grep_flu ( $aa, "contigs.blast");
	mkdir  "$outdir/2016-CX1608/first_map";
	system ("blastdbcmd -db $outdir/../fasta.fa -entry_batch $aa/blast_list > $outdir/2016-CX1608/blast.fa");
	system ("fasta_formatter -i $outdir/2016-CX1608/blast.fa -o $outdir/2016-CX1608/first_map/blast_w.fa -w 0");
	
	open (IN,"<$outdir/2016-CX1608/first_map/blast_w.fa") or die;
	open (OUT,">$outdir/2016-CX1608/first_map/2016-CX1608._first_ref.fa") or die;
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
	system ("bowtie2-build $outdir/2016-CX1608/first_map/2016-CX1608._first_ref.fa $outdir/2016-CX1608/first_map/2016-CX1608._first_ref.fa");

	system ("bowtie2 --very-sensitive-local  -x $outdir/2016-CX1608/first_map/2016-CX1608._first_ref.fa -1 $outdir/2016-CX1608/trim_1.fastq -2 $outdir/2016-CX1608/trim_2.fastq -S $outdir/2016-CX1608/first_map/2016-CX1608._first.sam -p 10");
	system ("samtools view -bS $outdir/2016-CX1608/first_map/2016-CX1608._first.sam > $outdir/2016-CX1608/first_map/2016-CX1608._first.bam");
	print "\n\n\n";
	system ("samtools sort $outdir/2016-CX1608/first_map/2016-CX1608._first.bam -o $outdir/2016-CX1608/first_map/2016-CX1608._first_sort.bam");
	system ("samtools index $outdir/2016-CX1608/first_map/2016-CX1608._first_sort.bam ");#$outdir/2016-CX1608/first_map/2016-CX1608._first_sort.bam.bai
	print "\n\n\n";
	mkdir  "$outdir/2016-CX1608/second_map"; 
	system ("samtools mpileup -Q 0 $outdir/2016-CX1608/first_map/2016-CX1608._first_sort.bam > $outdir/2016-CX1608/first_map/bbb");
	system ("perl $outdir/../pileup2fasta.pl -S 0 -d 1 $outdir/2016-CX1608/first_map/bbb > $outdir/2016-CX1608/second_map/first_consensus.fa");
	#system ("samtools mpileup -vuf  $outdir/2016-CX1608/first_map/2016-CX1608._first_ref.fa  $outdir/2016-CX1608/first_map/2016-CX1608._first_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/2016-CX1608/second_map/first_consensus.fa");
	#system ("bowtie2-build $outdir/2016-CX1608/second_map/first_consensus.fa $outdir/2016-CX1608/second_map/first_consensus.fa");

	
	
	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_1.fastq $outdir/2016-CX1608/unpair_1.fastq $outdir/2016-CX1608/trim_2.fastq $outdir/2016-CX1608/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	
	system ("fastqc -o $outdir/2016-CX1608 $outdir/2016-CX1608/trim_1.fastq $outdir/2016-CX1608/trim_2.fastq");
	
	
	
	system ("bowtie2 --very-sensitive-local -x /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -1 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_1.fastq -2 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_2.fastq  --un-conc /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.fq -S /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second.sam -p 10");
	#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
	#	-x <bt2-idx> 由bowtie2-build所生成的索引文件的前缀。首先 在当前目录搜寻，然后在环境变量 BOWTIE2_INDEXES 中制定的文件夹中搜寻。
	#-1 <m1> 双末端测寻对应的文件1。可以为多个文件，并用逗号分开；多个文件必须和 -2 <m2> 中制定的文件一一对应。比如:"-1 flyA_1.fq,flyB_1.fq -2 flyA_2.fq,flyB_2.fq". 测序文件中的reads的长度可以不一样。-2 <m2> 双末端测寻对应的文件2.
	#-U <r> 非双末端测寻对应的文件。可以为多个文件，并用逗号分开。测序文件中的reads的长度可以不一样。
	#-S <hit> 所生成的SAM格式的文件前缀。默认是输入到标准输出。
	#--very-sensitive-local Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
	#--un-conc <path> 将不能和谐比对的paired-end reads写入<path>
	system ("samtools view -bS $outdir/2016-CX1608/second_map/second.sam > $outdir/2016-CX1608/second_map/2016-CX1608._second.bam");
	#samtool view -bS命令把sam文件变成bam文件
	print "\n\n\n";
=cut
	system ("samtools view -bhS /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second.sam > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second.bam");
	system ("samtools sort /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam");
	#sort排序（不明白有啥用）
	#system ("samtools index /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam ");#/home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/first_map/2016-CX1608._first_sort.bam.bai
	#print "\n\n\n";
	#system("bowtie2 -rg-id sample -rg"PL:ILLYMINA" -rg "SM:sample" -x /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -1 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_1.fastq -2 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_2.fastq --un-conc /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.fq -S /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second.sam -p 10");
	#system("java -Xmx10g -jar /home/lei/Desktop/Chen/picard_tools/SortSam.jar INPUT=/home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second.sam OUTPUT= /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam SORT_ORDER=coordinate");
	##system ("samtools mpileup -Q 0 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second_pileup");
	#system ("java -Xmx10g -jar /home/lei/Desktop/Chen/picard_tools/MarkDuplicates.jar INPUT=2016-CX1608/second_map/2016-CX1608._second_sort.bam  OUTPUT=/home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam METRICS_FILE=/home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._dup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000");
	#标记PCR duplicates
	#samtools mpileup 
	system ("java -Xmx10g -jar /home/lei/Desktop/Chen/picard_tools/CreateSequenceDictionary.jar R=/home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa O=/home/lei/Desktop/flu_toolkit/zhuwenfei/reference.dict");
	system ("samtools index /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam");
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T RealignerTargetCreator -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.intervals");
	system ("java −Xmx10g−jarGATKHome/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T IndelRealigner -targetIntervals /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.intervals -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.bam");
	##以上4个system()为建立索引并重新比对
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T UnifiedGenotyper  -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._gatk.raw1.vcf --read_filter BadCigar -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 0");
	system ("samtools mpileup -DSugf /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.bam | bcftools view -Ncvg -> /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._samtools.raw1.vcf");
	#以上2个system()为分别生成GATK和samtools的vcf文件(SNP和INDEL calling)
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T SelectVariants --variant /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._gatk.raw1.vcf --concordance /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._samtools.raw1.vcf -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.raw1.vcf");
	system ("\$MEANQUAL = `awk '{ if (\$1! ~ /#/) { total += \$6; count++ }} END { print total/count }' /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.raw1.vcf`");
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T VariantFiltration --filterExpression 'QD < 20.0 || ReadPosRankSum < -8.0 ||  FS > 10.0 || QUAL < \$MEANQUAL' --filterName LowQualFilter --variant /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.raw1.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing  --logging_level ERROR -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.flt1.vcf");
	system ("grep -v "Filter" /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.flt1.vcf > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.filter1.vcf");
	#以上4个system()为对两个版本的vcf结果进行综合和过滤
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T BaseRecalibrator -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._recal_data.grp -knownSites /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.filter1.vcf");
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T PrintReads -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._realn.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._recal.bam -BQSR /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._recal_data.grp");
	#以上2个system()为对realign的bam文件进行校正，生成较准确的-knowSite
	#(以上6个system为第一遍calling)
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T UnifiedGenotyper -I /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._recal.bam -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._gatk.raw2.vcf  --read_filter BadCigar -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 0");
	system ("samtools mpileup -DSugf /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._recal.bam |  bcftools view -Ncvg -> /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._samtools.raw2.vcf");
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar  -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T SelectVariants  --variant /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._gatk.raw2.vcf --concordance /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._samtools.raw2.vcf -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.raw2.vcf");
	system ("java −Xmx10g −jar /home/lei/Desktop/flu_toolkit/flu_out/../../Chen/GATK/GenomeAnalysisTK.jar -R /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -T VariantFiltration --filterExpression 'QD < 10.0 || ReadPosRankSum < -8.0 ||  FS > 10.0 || QUAL < 30'  --filterName LowQualFilter --variant /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.raw2.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing  --logging_level ERROR -o /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.flt2.vcf");
	system ("grep -v "Filter" /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.flt2.vcf >  /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.filter2.vcf");
	system ("cp /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._concordance.filter2.vcf /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/callSNP/2016-CX1608._final.vcf");
	#以上6个system为第二遍callSNP，call出最终结果
	
	
=pod1	
	system ("samtools mpileup -Q 0 -f /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/ref_dep");
	system ("perl /home/lei/Desktop/flu_toolkit/flu_out/../draw_depth.pl /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/ref_dep > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/coverage");
	system ("Rscript /home/lei/Desktop/flu_toolkit/flu_out/../DR.r /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/");
	system ("perl /home/lei/Desktop/flu_toolkit/flu_out/../pileup2fasta.pl -d 3 -S 0.1 /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second_pileup >  /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/2016-CX1608.10.fa");
	#	system ("samtools mpileup -vuf  /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/first_consensus.fa  /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/2016-CX1608._second_sort.bam | bcftools view -cgS - | vcfutils.pl vcf2fq  > /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/consensus.fa");
	system ("java -jar /home/lei/Desktop/flu_toolkit/flu_out/../FluAnn.jar /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/2016-CX1608.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/2016-CX1608.html");
	# (fastq_to_fasta -Q33 -i te.fq -o te.fa) also works;
	
	open (IN,"</home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.1.fq") or die;
	open (OUT,">/home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.fa") or die;
	my $i = 1;
	while(<IN>){
		print OUT ">$_\n" if ($i%4 ==1);
		print OUT "$_\n"  if ($i%4 ==2);
		$i++;
	}
	close IN;
	
	open (IN,"</home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.2.fq") or die;
	 $i = 1;
	while(<IN>){
		print OUT ">$_\n" if ($i%4 ==1);
		print OUT "$_\n"  if ($i%4 ==2);
		$i++;
	}
	close IN;
	close OUT;
	#
	system ("blastn -query /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/unmap.fa -out /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/unmap.blast -db /home/lei/Desktop/flu_toolkit/flu_out/../fasta.fa -outfmt '6 stitle sseqid' -max_target_seqs 1 -num_threads 10");
	
	
=edit
	system ("cp $indir/$eachfile /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map");
	system ("cp /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/first_consensus.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608");
	system ("cp /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/second_map/second.sam /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608");
	
	system ("rm -rf /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/first_map");
	system ("rm -rf  /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/velvet");
	system ("rm /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/blast.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/cat_file.fa /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim_fastqc.zip /home/lei/Desktop/flu_toolkit/flu_out/2016-CX1608/trim.fastq");
=cut
}

print "PRESS ENTER";
