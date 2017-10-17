#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:02:12 2017

@author: Ja
这个脚本是统计second_pileup文件中各位置的碱基组成
用法
在flutooolkit文件夹下
python readsNTcount.py
"""
import re
import os
#import sys

cdmyoutdir = "cd /home/lei/Desktop/flu_toolkit/flu_out";
os.system(cdmyoutdir)
filelist = os.popen('ls /home/lei/Desktop/flu_toolkit/flu_out/')
filelist = list(filelist)
#print(filelist)
#os.getcwd()


for file in filelist:
    infilepath = '/home/lei/Desktop/flu_toolkit/flu_out/' + file.strip() + r'/second_map/second_pileup'
    outfilepath = '/home/lei/Desktop/flu_toolkit/flu_out/' + file.strip() + r'_readsNTcount'
    f = open(infilepath,'r')
    out = open(outfilepath,'w')
#with out& with f:
    out.write('Segment\tSite\tA count\tC count\tT count\tG count\tN count\n')
    s = f.readlines()
    for line in s:
        out.write(re.sub(r'_<.*$','',line.split('\t')[0]) + '\t' + line.split('\t')[1] + '\t')
        m = line.split('\t')[4].upper()
        pattern = '[+-](\d)'
        num = re.findall(pattern, m)
            
        for number in num:
            numpattern = '[+-]%s[ACTGN]{%s}' %(number,number)
            m = re.sub(numpattern,'', m)
        for nt in 'ACTGN':
            countnumber = m.count(nt)
            out.write('%s is: %d\t' %(nt, countnumber))
        out.write('\n')
    f.close()
    out.close()
