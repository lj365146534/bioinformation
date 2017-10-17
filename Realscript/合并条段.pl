#!/usr/bin/perl -w
#usage:perl perl.pl INFILE
#缺点：此脚本对第一个的>判断不为分隔符···，所以会出现第一个序列的名字变为“>>seqname”，因此使用前先把各个片段copy *.fas infile.fa整合成1个fa,再删除infile第一个序列的>就可以了。
#刘佳写于2017-03-07
###################################################################################################################################################################
######################################################################################################################################################
use strict;
my %hash;
my $infile = $ARGV[0];
$/ = ">";
open (DATA, "$infile")|| die;
open (OUTFILE, ">mm.fas")||die;
$/ = "\n";
while (<DATA>){
my $name = $1 if (/(\S+)\n/);
$/ = ">";
my $seq = <DATA>;
$seq =~ s/>$|\n|\r//g;
push @{$hash{$name}},$seq;
$/ = "\n";
}
foreach my $k (sort keys %hash){
print OUTFILE ">$k\n";
print OUTFILE join "",@{$hash{$k}};
print OUTFILE "\n";
}
close DATA;
close OUTFILE;
