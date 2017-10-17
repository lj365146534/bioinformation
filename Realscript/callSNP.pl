#! /usr/bin/perl -w
use strict;
use warnings;
use grep_flu_fregment;

my $outdir = "/home/lei/Desktop/flu_toolkit/flu_out";
my $indir  = "/home/lei/Desktop/flu_toolkit/flu_in";
my $threads = 10;

#my $MEANQUAL;
my $Filter = "\"Filter\"";
my $PL = "\"PL:ILLUMINA\"";
my $SM = "\"SM:sample\"";
my $QD = "\"QD < 20.0 \|\| ReadPosRankSum < -8.0 \|\|  FS > 10.0 \|\| QUAL < \$MEANQUAL\"";


chdir $indir;
my @filelist = split (/\s/ , `ls`);
print "@filelist\n";
#print $filelist[1];
foreach my $eachfile (@filelist){
	chomp $eachfile;
	my $basename = $eachfile;  
#for	(my $basename = "ZWF-01"){  
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
	mkdir "$outdir/$basename/second_map";
	mkdir "$outdir/$basename/callSNP";
	chdir "$outdir/$basename/callSNP/";
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
	system ("samtools view -bS $outdir/$basename/first_map/$basename._first.sam > $outdir/$basename/first_map/$basename._first.bam");
	print "\n\n\n";
	system ("samtools sort $outdir/$basename/first_map/$basename._first.bam -o $outdir/$basename/first_map/$basename._first_sort.bam");
	system ("samtools index $outdir/$basename/first_map/$basename._first_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	print "\n\n\n";
	mkdir  "$outdir/$basename/second_map"; 
	system ("samtools mpileup -Q 0 $outdir/$basename/first_map/$basename._first_sort.bam > $outdir/$basename/first_map/bbb");
	system ("perl $outdir/../pileup2fasta.pl -S 0 -d 1 $outdir/$basename/first_map/bbb > $outdir/$basename/second_map/first_consensus.fa");
	#system ("samtools mpileup -vuf  $outdir/$basename/first_map/$basename._first_ref.fa  $outdir/$basename/first_map/$basename._first_sort.bam | /home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools/bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/second_map/first_consensus.fa");
	#system ("bowtie2-build $outdir/$basename/second_map/first_consensus.fa $outdir/$basename/second_map/first_consensus.fa");

	
	
	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq $outdir/$basename/trim_1.fastq $outdir/$basename/unpair_1.fastq $outdir/$basename/trim_2.fastq $outdir/$basename/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	
	system ("fastqc -o $outdir/$basename $outdir/$basename/trim_1.fastq $outdir/$basename/trim_2.fastq");
	
	
	
#	system ("bowtie2 --very-sensitive-local -x /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa -1 $outdir/$basename/trim_1.fastq -2 $outdir/$basename/trim_2.fastq  --un-conc $outdir/$basename/second_map/unmap.fq -S second.sam -p $threads");
	#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
	#	-x <bt2-idx> 由bowtie2-build所生成的索引文件的前缀。首先 在当前目录搜寻，然后在环境变量 BOWTIE2_INDEXES 中制定的文件夹中搜寻。
	#-1 <m1> 双末端测寻对应的文件1。可以为多个文件，并用逗号分开；多个文件必须和 -2 <m2> 中制定的文件一一对应。比如:"-1 flyA_1.fq,flyB_1.fq -2 flyA_2.fq,flyB_2.fq". 测序文件中的reads的长度可以不一样。-2 <m2> 双末端测寻对应的文件2.
	#-U <r> 非双末端测寻对应的文件。可以为多个文件，并用逗号分开。测序文件中的reads的长度可以不一样。
	#-S <hit> 所生成的SAM格式的文件前缀。默认是输入到标准输出。
	#--very-sensitive-local Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
	#--un-conc <path> 将不能和谐比对的paired-end reads写入<path>
	system ("samtools view -bS -h second.sam > $outdir/$basename/second_map/$basename._second.bam");#-h输出头文件+map文件
	#samtool view -bS命令把sam文件变成bam文件
	print "\n\n\n";

	system ("java -jar /home/lei/Desktop/Chen/trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 10 $indir/$base_1.fastq $indir/$base_2.fastq $outdir/$basename/trim_1.fastq $outdir/$basename/unpair_1.fastq $outdir/$basename/trim_2.fastq $outdir/$basename/unpair_2.fastq ILLUMINACLIP:/home/lei/Desktop/Chen/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	
	system ("fastqc -o $outdir/$basename $outdir/$basename/trim_1.fastq $outdir/$basename/trim_2.fastq");
	
	system ("cp $outdir/$basename/$basename.fa $outdir/$basename/callSNP/$basename.fa");
	system ("bowtie2-build $basename.fa $basename");
	system ("bowtie2 -p $threads --rg-id sample --rg $PL --rg $SM -x $basename -1 $outdir/$basename/trim_1.fastq -2 $outdir/$basename/trim_2.fastq -S second.sam");
	#system ("samtools view -bhS second.sam > $outdir/$basename/second_map/$basename._second.bam");#-h输出头文件+map文件
	system ("java -Xmx10g -jar $outdir/../../Chen/picard_tools/SortSam.jar INPUT=second.sam OUTPUT=$basename.second_sort.bam SORT_ORDER=coordinate");
	#sort排序（不明白有啥用）
	#system ("samtools index $outdir/$basename/second_map/$basename._second_sort.bam ");#$outdir/$basename/first_map/$basename._first_sort.bam.bai
	#print "\n\n\n";
	##system ("samtools mpileup -Q 0 $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/second_pileup");
	system ("java -Xmx10g -jar $outdir/../../Chen/picard_tools/MarkDuplicates.jar INPUT=$basename.second_sort.bam  OUTPUT=$basename.dup.bam METRICS_FILE=$basename.dup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000");
	#标记PCR duplicates
	#samtools mpileup 
	system ("samtools faidx $outdir/$basename/callSNP/$basename.fa");
	system ("java -Xmx10g -jar $outdir/../../Chen/picard_tools/CreateSequenceDictionary.jar R=$outdir/$basename/callSNP/$basename.fa O=$outdir/$basename/callSNP/$basename.dict");
	system ("samtools index $outdir/$basename/callSNP/$basename.dup.bam");
	
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T RealignerTargetCreator -I $basename.dup.bam -o $outdir/$basename/callSNP/$basename.realn.intervals");
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T IndelRealigner -targetIntervals $basename.realn.intervals -I $basename.dup.bam -o $basename.realn.bam");
	##以上4个system()为建立索引并重新比对
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T UnifiedGenotyper  -I $basename.realn.bam -o $basename.gatk.raw1.vcf --read_filter BadCigar -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 0");

	system ("samtools mpileup -t DP -t SP -gf $basename.fa $basename.realn.bam > $basename.samtools.raw1.bcf.gz");#有.GZ
	system ("/home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools call -cvO z $basename.samtools.raw1.bcf.gz -o $basename.samtools.raw1.vcf");
	#system ("samtools mpileup -t DP -t SP -gf $basename.fa $basename.realn.bam | /home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools/bcftools view -Go -> $basename.samtools.raw1.vcf");
	#system ("samtools mpileup -t DP -t SP -guf $basename.fa $basename.realn.bam | /home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools/bcftools call -vmO v -> $basename.samtools.raw1.vcf");
	system ("/home/lei/Desktop/liujia/bcftools-1.3/vcfutils.pl varFilter $basename.samtools.raw1.vcf > final_variants.vcf");
	#以上2个system()为分别生成GATK和samtools的vcf文件(SNP和INDEL calling)
=pod
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T SelectVariants -V $basename.gatk.raw1.vcf -conc $basename.samtools.raw1.vcf -o $outdir/$basename/callSNP/$basename.concordance.raw1.vcf");
	#system ("\$MEANQUAL=`awk '{ if (\$1! ~ /#/) { total += \$6; count++ }} END { print total/count }' $basename.concordance.raw1.vcf`");
	#system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T VariantFiltration --filterExpression $QD --filterName LowQualFilter --variant $basename.concordance.raw1.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing  --logging_level ERROR -o $basename.concordance.flt1.vcf");
	#system ("grep -v $Filter $basename.concordance.flt1.vcf > $basename.concordance.filter1.vcf");
	#以上4个system()为对两个版本的vcf结果进行综合和过滤
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T BaseRecalibrator -I $basename.realn.bam -o $basename.recal_data.grp -knownSites $basename.concordance.filter1.vcf");
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T PrintReads -I $basename.realn.bam -o $basename.recal.bam -BQSR $basename.recal_data.grp");
	#以上2个system()为对realign的bam文件进行校正，生成较准确的-knowSite
	#(以上6个system为第一遍calling)
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T UnifiedGenotyper -I $basename.recal.bam -o $basename.gatk.raw2.vcf  --read_filter BadCigar -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 0");
	system ("samtools mpileup -DSugf $basename.recal.bam |  /home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools/bcftools call -cvO z -> $basename.samtools.raw2.vcf.gz");
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar  -R $basename.fa -T SelectVariants  --variant $basename.gatk.raw2.vcf --concordance $basename.samtools.raw2.vcf -o $basename.concordance.raw2.vcf");
	system ("java -Xmx10g -jar $outdir/../../Chen/GATK/GenomeAnalysisTK.jar -R $basename.fa -T VariantFiltration --filterExpression 'QD < 10.0 || ReadPosRankSum < -8.0 ||  FS > 10.0 || QUAL < 30'  --filterName LowQualFilter --variant $basename.concordance.raw2.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing  --logging_level ERROR -o $basename.concordance.flt2.vcf");
	system ("grep -v $Filter $basename.concordance.flt2.vcf >  $basename.concordance.filter2.vcf");
	system ("cp $basename.concordance.filter2.vcf $basename.final.vcf");
	#以上6个system为第二遍callSNP，call出最终结果
=cut	
	
=pod1	
	system ("samtools mpileup -Q 0 -f /home/lei/Desktop/flu_toolkit/zhuwenfei/reference.fa $outdir/$basename/second_map/$basename._second_sort.bam > $outdir/$basename/second_map/ref_dep");
	system ("perl $outdir/../draw_depth.pl $outdir/$basename/second_map/ref_dep > $outdir/$basename/second_map/coverage");
	system ("Rscript $outdir/../DR.r $outdir/$basename/");
	system ("perl $outdir/../pileup2fasta.pl -d 3 -S 0.1 $outdir/$basename/second_map/second_pileup >  $outdir/$basename/$basename.10.fa");
	#	system ("samtools mpileup -vuf  $outdir/$basename/second_map/first_consensus.fa  $outdir/$basename/second_map/$basename._second_sort.bam | /home/lei/Desktop/liujia/bcftools-1.3/bin/bcftools/bcftools view -cgS - | vcfutils.pl vcf2fq  > $outdir/$basename/consensus.fa");
	system ("java -jar $outdir/../FluAnn.jar $basename.fa $outdir/$basename/$basename.html");
	# (fastq_to_fasta -Q33 -i te.fq -o te.fa) also works;
	
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
	system ("cp second.sam $outdir/$basename");
	
	system ("rm -rf $outdir/$basename/first_map");
	system ("rm -rf  $outdir/$basename/velvet");
	system ("rm $outdir/$basename/blast.fa $outdir/$basename/cat_file.fa $outdir/$basename/trim_fastqc.zip $outdir/$basename/trim.fastq");
=cut
}

print "PRESS ENTER";
