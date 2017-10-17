#!/usr/bin/perl -w

##########################################################################################
# This script find ambiguous sites from fasta files for each seq
# Author: Jia Liu
# Date: 2016-11-14
# Annotation: 仅供刘佳及其朋友学习交流使用（致敬陈连福），用前推荐用mega对齐···也可不对齐
# 2016-12-15修改，增加提取氨基酸位置和三联码
# Usage: perl find_ambiguous_pro.pl [FILE]
##########################################################################################



#my@display_name = "";
#my@desc = "";
#my@seq = "";


my $display_name;
my $seq;
my $position;
my $seq_length;
my $infile = $ARGV[0];
my $nt_pos;
my $triplet_codon;
my $AA_pos;
open (BABABABA,$infile)||die("can't open file BABABABA, $!");         # 打开ecorho.fasta文件
open (CCCCCCCC, ">$infile.txt")||die("can't open file CCCCCCCC, $!");
#foreach my$numb(0..$#index)
print CCCCCCCC "display_name\tposition\tambigous\tlength\ttriplet_codon\tAA_pos\n";
foreach(<BABABABA>){
			$/ = "(\n>|>)";
			#上一句把系统默认分隔符\n改变为>或者\n>;
		
			chomp;
			#print $_;
			#@ccc = $_;
			#print CCCCCCCC @ccc;
			if(/^>/){   
				($display_name,$seq)=/>(.*?)\n(.*)$/;
				#上一句把毒株名字和序列分别赋值给$display_name和$seq;
				#也可以把($display_name,$seq)=/>(.*?)\n(.*)$/换成($display_name,$seq)=/>(.*?)\n(.*\n)$/，这样可以不把$/ = "(\n>|>)"，只设成">"就可以
				
			}
			else{
        $seq = $_;
		#print $seq;
			}
			#print $display_name;
			while($seq =~ m/(r|y|m|k|s|w|h|b|v|d|n)/ig){
			#匹配简并碱基；
			$seq_length = length($seq)-1;#序列全长；
			$position = pos($seq);#简并碱基的位置；
			my $ambigous = substr($seq,$position-1,1);#截取简并碱基；
			$nt_pos = $position % 3;
				if($nt_pos == 1){
					$triplet_codon = substr($seq, $position-1,3);
					$AA_pos = ($position+2)/3;
					}
				elsif($nt_pos == 2){
					 $triplet_codon = substr($seq, $position-2,3);
					 $AA_pos = ($position+1)/3;
					}
				elsif($nt_pos == 0){
					 $triplet_codon = substr($seq, $position-3,3);
					 $AA_pos = ($position)/3;
					}
			
			print CCCCCCCC "$display_name\t$position\t$ambigous\t$seq_length\t$triplet_codon\t$AA_pos\n";
			}
			print "$display_name\n";
		}
=cut
			print $;
			my ( $pos, $now ) = ( 0, -1 );
			until ( $pos == -1 ) {
						$pos = index( $seq, $key, $now + 1 );
						$now = $pos;
						say $pos unless $pos < 0;
=cut
						
			
    #  如果这一行的开头是>，就说明是描述的一行，可以提取序列的名称和描述
=cut			if($ccc =~ /^>/){    
	# 从>开始到第一个空白出现为“名称”，之后的内容为“描述”
			#($display_name[$numb],$seq[$numb])=/^>(.*?)\n(^actg)$/i;
			#print $display_name[$numb];
			#print $seq[$numb];
			($display_name,$seq)=/^>(.*?)\n(.*?)$/;
			print $display_name;
			print $seq;

=cut   #  如果这一行的开头不是>，则就是序列行。由于序列可以分为好几行，所以要把每一行的序列都连接起来。


			$/ = "\n";

#  现在“三要素”已经提取出来，可以进一步分析了。先来计算序列的长度。
			#$seq_length =length($seq[$numb]);
#  我们可以判断序列是DNA还是蛋白质
			#if ($seq[$numb]=~/[^atgc]/i){
			#$seq_type="protein";
#}			else{
#			$seq_type="dna";
#}

#print "sequence name: $display_name\n";
#print "sequence :$seq\n";
#print "sequence description: $desc[$numb]";
#print "sequence type: $seq_type\t";
#print "sequence lengt: $seq_length";



close BABABABA;
close CCCCCCCC;
print ("DONE!!!PRESS 666 TO CONTINUE");