#!/usr/bin/perl
use strict;    
use warnings; 


my $secretnumber=int(1 + rand 100);  
#print "$secretnumber\n";  
{#一个裸块  
print "Please guess a number from 1 to 100!\n";  
my $gussnumber;
chomp($gussnumber=<STDIN>);  
if ($gussnumber=~/quit|exit|\A\s*\z/i)#必须要放在前面，如果放在后面字符串会被当做0  
 {  
	print "Wrong format!!!Be careful!!";
     last;  
 }  
  
elsif ($secretnumber<$gussnumber)  
{  
   print "Too hight\n";  
     redo;  
 }  
 elsif ($secretnumber>$gussnumber)  
 {  
     print "Too low.\n";  
     redo;  
 }  
 elsif ($secretnumber==$gussnumber)  
 {  
     print "You win!\n";  
      
	  
	   last;
 }  
 }
 print "The secretnumber is $secretnumber!\n";
 print "press <Enter> to continue...";
<STDIN> 
 