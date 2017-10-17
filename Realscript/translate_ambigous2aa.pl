#! /usr/bin/perl

##########################################################################################
# This script translate seq with ambiguous from fasta files 
# Author: Jia Liu
# Date: 2017-01-09
# Annotation: 仅供刘佳及其朋友学习交流使用,用前推荐用mega对齐···也可不对齐
# Usage: perl translate_ambiguous2aa.pl [FILE.fasta]			  生成FILE.fasta.txt文件
##########################################################################################


use strict;    
use warnings;   
#my $dna      ='';    
my $protein;    
my($display_name,$seq); 
my $infile = $ARGV[0];





open (BABABABA,$infile)||die("can't open file BABABABA, $!");         # 打开ecorho.fasta文件
open (CCCCCCCC, ">$infile.prot.txt")||die("can't open file CCCCCCCC, $!");

foreach(<BABABABA>){
			$/ = "(\n>|>)";
			#上一句把系统默认分隔符\n改变为>或者\n>;
		
			chomp;
			#print $_;
			#@ccc = $_;
			#print CCCCCCCC @ccc;
			if(/^>/){   
				($display_name,$seq)=/>(.*?)\n(.*)$/;
				#($seq) = /^(.*)$/;
				#上一句把毒株名字和序列分别赋值给$display_name和$seq;
				#也可以把($display_name,$seq)=/>(.*?)\n(.*)$/换成($display_name,$seq)=/>(.*?)\n(.*\n)$/,这样可以不把$/ = "(\n>|>)",只设成">"就可以
				
			}
			else{
				$seq = $_;
			}
		
		#print $seq;
		#$sequence = "$display_name\n$seq\n";
		#print $sequence;
		if($seq =~ /[ATCG]+/){		
		for(my $i = 0;$i<(length($seq)-3);$i+=3)    {    
		
         $protein .= codon2aa(substr($seq,$i,3)); 
    }
	print CCCCCCCC 	"\>$display_name\n$protein\n";
	$protein = "";			#这个弄了半天，每次都变长，我的天，最后归零了才搞定。
	#print CCCCCCCC  "$seq\n"
	}
	$/ = "\n";	
}	
#$/ = "\n";
#print CCCCCCCC 	"$display_name\n$protein\n";
close BABABABA;
close CCCCCCCC;


sub codon2aa       
{       
    
    #第三种方法      
    #也就是运用哈希      
    #我们将所有的密码子作为hash的key,然后将代表的氨基酸作为hash的value      
    #然后进行匹配      
    # codon2aa       
    # A subroutine to translate a DNA 3-character codon to an amino acid       
    # Version 3, using hash lookup       
    my($codon) = @_;       
       
    $codon = uc $codon;#uc=uppercase;lc=lowercase      
                   #也就是大小写转换,uc表示将所有的小写 转换为大写      
               #lc将所有的大写转换为小写      
        
    my(%genetic_code) = (       
           
    'TCA' => 'S',    # Serine       
    'TCC' => 'S',    # Serine       
    'TCG' => 'S',    # Serine       
    'TCT' => 'S',    # Serine       
    'TTC' => 'F',    # Phenylalanine       
    'TTT' => 'F',    # Phenylalanine       
    'TTA' => 'L',    # Leucine       
    'TTG' => 'L',    # Leucine       
    'TAC' => 'Y',    # Tyrosine        
    'TAT' => 'Y',    # Tyrosine       
    'TAA' => '*',    # Stop       
    'TAG' => '*',    # Stop       
    'TGC' => 'C',    # Cysteine       
    'TGT' => 'C',    # Cysteine       
    'TGA' => '*',    # Stop       
    'TGG' => 'W',    # Tryptophan       
    'CTA' => 'L',    # Leucine       
    'CTC' => 'L',    # Leucine       
    'CTG' => 'L',    # Leucine       
    'CTT' => 'L',    # Leucine       
    'CCA' => 'P',    # Proline       
    'CCC' => 'P',    # Proline       
    'CCG' => 'P',    # Proline       
    'CCT' => 'P',    # Proline       
    'CAC' => 'H',    # Histidine       
    'CAT' => 'H',    # Histidine       
    'CAA' => 'Q',    # Glutamine       
    'CAG' => 'Q',    # Glutamine       
    'CGA' => 'R',    # Arginine       
    'CGC' => 'R',    # Arginine       
    'CGG' => 'R',    # Arginine       
    'CGT' => 'R',    # Arginine       
    'ATA' => 'I',    # Isoleucine       
    'ATC' => 'I',    # Isoleucine       
    'ATT' => 'I',    # Isoleucine       
    'ATG' => 'M',    # Methionine       
    'ACA' => 'T',    # Threonine       
    'ACC' => 'T',    # Threonine       
    'ACG' => 'T',    # Threonine       
    'ACT' => 'T',    # Threonine       
    'AAC' => 'N',    # Asparagine       
    'AAT' => 'N',    # Asparagine       
    'AAA' => 'K',    # Lysine       
    'AAG' => 'K',    # Lysine       
    'AGC' => 'S',    # Serine       
    'AGT' => 'S',    # Serine       
    'AGA' => 'R',    # Arginine       
    'AGG' => 'R',    # Arginine       
    'GTA' => 'V',    # Valine       
    'GTC' => 'V',    # Valine       
    'GTG' => 'V',    # Valine       
    'GTT' => 'V',    # Valine       
    'GCA' => 'A',    # Alanine       
    'GCC' => 'A',    # Alanine       
    'GCG' => 'A',    # Alanine       
    'GCT' => 'A',    # Alanine           
    'GAC' => 'D',    # Aspartic Acid       
    'GAT' => 'D',    # Aspartic Acid       
    'GAA' => 'E',    # Glutamic Acid       
    'GAG' => 'E',    # Glutamic Acid       
    'GGA' => 'G',    # Glycine       
    'GGC' => 'G',    # Glycine       
    'GGG' => 'G',    # Glycine       
    'GGT' => 'G',    # Glycine       
	#############################################
	#############################################
	#	以下为带简并碱基且氨基酸不变的组合		#
	#											#
	#############################################
	#############################################
	
	'TGY' => 'C',  
	'GAY' => 'D',  
	'GAR' => 'E',  
	'TTY' => 'F',  
	'CAY' => 'H',  
	'ATH' => 'I',  
	'ATM' => 'I',  
	'ATW' => 'I',  
	'ATY' => 'I',  
	'AAR' => 'K',  
	'TTR' => 'L',  
	'YTA' => 'L',  
	'YTG' => 'L',  
	'YTR' => 'L',  
	'AAY' => 'N',  
	'CAR' => 'Q',  
	'AGR' => 'R',  
	'MGA' => 'R',  
	'MGG' => 'R',  
	'MGR' => 'R',  
	'AGY' => 'S',  
	'TAY' => 'Y',  
	'TAR' => '*',  
	'GCB' => 'A',  
	'GCF' => 'A',  
	'GCH' => 'A',  
	'GCK' => 'A',  
	'GCM' => 'A',  
	'GCN' => 'A',  
	'GCR' => 'A',  
	'GCS' => 'A',  
	'GCV' => 'A',  
	'GCW' => 'A',  
	'GCY' => 'A',  
	'GGB' => 'G',  
	'GGF' => 'G',  
	'GGH' => 'G',  
	'GGK' => 'G',  
	'GGM' => 'G',  
	'GGN' => 'G',  
	'GGR' => 'G',  
	'GGS' => 'G',  
	'GGV' => 'G',  
	'GGW' => 'G',  
	'GGY' => 'G',  
	'CTB' => 'L',  
	'CTF' => 'L',  
	'CTH' => 'L',  
	'CTK' => 'L',  
	'CTM' => 'L',  
	'CTN' => 'L',  
	'CTR' => 'L',  
	'CTS' => 'L',  
	'CTV' => 'L',  
	'CTW' => 'L',  
	'CTY' => 'L',  
	'CCB' => 'P',  
	'CCF' => 'P',  
	'CCH' => 'P',  
	'CCK' => 'P',  
	'CCM' => 'P',  
	'CCN' => 'P',  
	'CCR' => 'P',  
	'CCS' => 'P',  
	'CCV' => 'P',  
	'CCW' => 'P',  
	'CCY' => 'P',  
	'CGB' => 'R',  
	'CGF' => 'R',  
	'CGH' => 'R',  
	'CGK' => 'R',  
	'CGM' => 'R',  
	'CGN' => 'R',  
	'CGR' => 'R',  
	'CGS' => 'R',  
	'CGV' => 'R',  
	'CGW' => 'R',  
	'CGY' => 'R',  
	'TCB' => 'S',  
	'TCF' => 'S',  
	'TCH' => 'S',  
	'TCK' => 'S',  
	'TCM' => 'S',  
	'TCN' => 'S',  
	'TCR' => 'S',  
	'TCS' => 'S',  
	'TCV' => 'S',  
	'TCW' => 'S',  
	'TCY' => 'S',  
	'ACB' => 'T',  
	'ACF' => 'T',  
	'ACH' => 'T',  
	'ACK' => 'T',  
	'ACM' => 'T',  
	'ACN' => 'T',  
	'ACR' => 'T',  
	'ACS' => 'T',  
	'ACV' => 'T',  
	'ACW' => 'T',  
	'ACY' => 'T',  
	'GTB' => 'V',  
	'GTF' => 'V',  
	'GTH' => 'V',  
	'GTK' => 'V',  
	'GTM' => 'V',  
	'GTN' => 'V',  
	'GTR' => 'V',  
	'GTS' => 'V',  
	'GTV' => 'V',  
	'GTW' => 'V',  
	'GTY' => 'V', 
	'TAR' => '*',
	'TRA' => '*',
	);       
     

    if(exists $genetic_code{$codon})       
    {       
        return $genetic_code{$codon};       
    }      
    else      
    {       
        $genetic_code{$codon} = "\?";   
		return $genetic_code{$codon};			
    }
	
} 






	





=pod
sub triplet_codon2ATCG{
my $amb_triplet = @_;
my %ambiguous = (

'B' => 'C'
'B' => 'G'
'B' => 'T'
'D' => 'A'
'D' => 'G'
'D' => 'T'
'H' => 'A'
'H' => 'C'
'H' => 'T'
'K' => 'G'
'K' => 'T'
'M' => 'A'
'M' => 'C'
'N' => 'A'
'N' => 'C'
'N' => 'G'
'N' => 'T'
'R' => 'A'
'R' => 'G'
'S' => 'C'
'S' => 'G'
'V' => 'A'
'V' => 'C'
'V' => 'G'
'W' => 'A'
'W' => 'T'
'Y' => 'C'
'Y' => 'T'
);

my %ambiguous = (

'B' => 'CGT'
'D' => 'AGT'
'H' => 'ACT'
'K' => 'GT'
'M' => 'AC'
'N' => 'ACGT'
'R' => 'AG'
'S' => 'CG'
'V' => 'ACG'
'W' => 'AT'
'Y' => 'CT'
);

sub codon2aa   
{   
       my($codon) = @_;   
    
       if ( $codon =~ /GC./i)        { return 'A' }    # Alanine       
    elsif ( $codon =~ /TG[TCY]/i)     { return 'C' }    # Cysteine   
    elsif ( $codon =~ /GA[TCY]/i)     { return 'D' }    # Aspartic Acid   
    elsif ( $codon =~ /GA[AGR]/i)     { return 'E' }    # Glutamic Acid   
    elsif ( $codon =~ /TT[TCY]/i)     { return 'F' }    # Phenylalanine   
    elsif ( $codon =~ /GG./i)        { return 'G' }    # Glycine   
    elsif ( $codon =~ /CA[TCY]/i)     { return 'H' }    # Histidine   
    elsif ( $codon =~ /AT[TCAHMWY]/i)    { return 'I' }    # Isoleucine   
    elsif ( $codon =~ /AA[AGR]/i)     { return 'K' }    # Lysine   
    elsif ( $codon =~ /TT[AGR]|CT./i) { return 'L' }    # Leucine   
    elsif ( $codon =~ /ATG/i)        { return 'M' }    # Methionine   
    elsif ( $codon =~ /AA[TCY]/i)     { return 'N' }    # Asparagine   
    elsif ( $codon =~ /CC./i)        { return 'P' }    # Proline   
    elsif ( $codon =~ /CA[AGR]/i)     { return 'Q' }    # Glutamine   
    elsif ( $codon =~ /CG.|AG[AGR]/i) { return 'R' }    # Arginine   
    elsif ( $codon =~ /TC.|AG[TCY]/i) { return 'S' }    # Serine   
    elsif ( $codon =~ /AC./i)        { return 'T' }    # Threonine   
    elsif ( $codon =~ /GT./i)        { return 'V' }    # Valine   
    elsif ( $codon =~ /TGG/i)        { return 'W' }    # Tryptophan   
    elsif ( $codon =~ /TA[TCY]/i)     { return 'Y' }    # Tyrosine   
    elsif ( $codon =~ /TA[AGR]|TGA/i) { return '*' }    # Stop   
    else   
    {   
        return '?';   
        exit;   
    }   
}   



sub codon2aa       
{       
    
    #第三种方法      
    #也就是运用哈希      
    #我们将所有的密码子作为hash的key,然后将代表的氨基酸作为hash的value      
    #然后进行匹配      
    # codon2aa       
    # A subroutine to translate a DNA 3-character codon to an amino acid       
    # Version 3, using hash lookup       
    my($codon) = @_;       
       
    $codon = uc $codon;#uc=uppercase;lc=lowercase      
                   #也就是大小写转换,uc表示将所有的小写 转换为大写      
               #lc将所有的大写转换为小写      
        
    my(%genetic_code) = (       
           
    'TCA' => 'S',    # Serine       
    'TCC' => 'S',    # Serine       
    'TCG' => 'S',    # Serine       
    'TCT' => 'S',    # Serine       
    'TTC' => 'F',    # Phenylalanine       
    'TTT' => 'F',    # Phenylalanine       
    'TTA' => 'L',    # Leucine       
    'TTG' => 'L',    # Leucine       
    'TAC' => 'Y',    # Tyrosine        
    'TAT' => 'Y',    # Tyrosine       
    'TAA' => '*',    # Stop       
    'TAG' => '*',    # Stop       
    'TGC' => 'C',    # Cysteine       
    'TGT' => 'C',    # Cysteine       
    'TGA' => '*',    # Stop       
    'TGG' => 'W',    # Tryptophan       
    'CTA' => 'L',    # Leucine       
    'CTC' => 'L',    # Leucine       
    'CTG' => 'L',    # Leucine       
    'CTT' => 'L',    # Leucine       
    'CCA' => 'P',    # Proline       
    'CCC' => 'P',    # Proline       
    'CCG' => 'P',    # Proline       
    'CCT' => 'P',    # Proline       
    'CAC' => 'H',    # Histidine       
    'CAT' => 'H',    # Histidine       
    'CAA' => 'Q',    # Glutamine       
    'CAG' => 'Q',    # Glutamine       
    'CGA' => 'R',    # Arginine       
    'CGC' => 'R',    # Arginine       
    'CGG' => 'R',    # Arginine       
    'CGT' => 'R',    # Arginine       
    'ATA' => 'I',    # Isoleucine       
    'ATC' => 'I',    # Isoleucine       
    'ATT' => 'I',    # Isoleucine       
    'ATG' => 'M',    # Methionine       
    'ACA' => 'T',    # Threonine       
    'ACC' => 'T',    # Threonine       
    'ACG' => 'T',    # Threonine       
    'ACT' => 'T',    # Threonine       
    'AAC' => 'N',    # Asparagine       
    'AAT' => 'N',    # Asparagine       
    'AAA' => 'K',    # Lysine       
    'AAG' => 'K',    # Lysine       
    'AGC' => 'S',    # Serine       
    'AGT' => 'S',    # Serine       
    'AGA' => 'R',    # Arginine       
    'AGG' => 'R',    # Arginine       
    'GTA' => 'V',    # Valine       
    'GTC' => 'V',    # Valine       
    'GTG' => 'V',    # Valine       
    'GTT' => 'V',    # Valine       
    'GCA' => 'A',    # Alanine       
    'GCC' => 'A',    # Alanine       
    'GCG' => 'A',    # Alanine       
    'GCT' => 'A',    # Alanine           
    'GAC' => 'D',    # Aspartic Acid       
    'GAT' => 'D',    # Aspartic Acid       
    'GAA' => 'E',    # Glutamic Acid       
    'GAG' => 'E',    # Glutamic Acid       
    'GGA' => 'G',    # Glycine       
    'GGC' => 'G',    # Glycine       
    'GGG' => 'G',    # Glycine       
    'GGT' => 'G',    # Glycine       
    );       
       
    if(;exists $genetic_code{$codon})       
    {       
        return $genetic_code{$codon};       
    }elsif{
	
	
	}      
    else      
    {       
       
            print STDERR "Bad codon \"$codon\"!!\n";       
            exit;       
    }       
}       



foreach($triplet_codon){
 if($amb =~ m/B/ig)       
    {       
        $triplet_codon =~ s/B/(C|G|T)/);       
    }      
    elsif($triplet_codon =~ m/D/ig)      
    {       
        $triplet_codon =~ s/D/(A|G|T)/);    
    } 
	elsif($triplet_codon =~ m/h/ig)        
    {       
        $triplet_codon =~ s/H/(A|C|T)/);    
    } 
	elsif($triplet_codon =~ m/k/ig)        
    {       
        $triplet_codon =~ s/K/(G|T)/);    
    } 
	elsif($triplet_codon =~ m/m/ig)        
    {       
        $triplet_codon =~ s/M/(A|C)/);    
    } 
	elsif($triplet_codon =~ m/n/ig)        
    {       
        $triplet_codon =~ s/N/[$ambiguous{N}]/);    
    } 
	elsif($triplet_codon =~ m/r/ig)        
    {       
        $triplet_codon =~ s/R/[AG]/);    
    } 
	elsif($triplet_codon =~ m/s/ig)        
    {       
        $triplet_codon =~ s/S/[CG]/);    
    } 
	elsif($triplet_codon =~ m/v/ig)        
    {       
        $triplet_codon =~ s/V/[AGC]/);    
    } 
	elsif($triplet_codon =~ m/w/ig)        
    {       
        $triplet_codon =~ s/W/[AT])/);    
    } 
	elsif($triplet_codon =~ m/y/ig)        
    {       
        $triplet_codon =~ s/Y/[CT]/);    
    } 
}
}


=cut



	