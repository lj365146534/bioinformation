#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:02:12 2017

@author: Ja
这个脚本是用来批量samtools view bam文件的，结果在当前文件夹下的out.sam
需要调用我自己写的PickSeqFromSam包，使用这个包，把out.sam变成fas文件
最后用mafft比对，mega手动选取
用法
python samtoolsview.py HA_\\\<H7N9\\\>:1-2
"""
import os
import sys
from PickSeqFromSam import sam2fas

cdmyoutdir = "cd /home/lei/Desktop/flu_toolkit/flu_out";
os.system(cdmyoutdir)

scaffoldrigion = str(sys.argv[1])

filelist = os.popen('ls /home/lei/Desktop/flu_toolkit/flu_out/')
filelist = list(filelist)
print(filelist)
outfilelist = []
for file in filelist:
    #print(file)
    file = file.strip()
    filepath = file + r'/second_map/'
    outfilelist.append(file.strip() + '-' + scaffoldrigion.replace('\\','') + '.out.sam')
    sort_bam_path = filepath + file +'._second_sort.bam'
    samtoolscommand = 'samtools view /home/lei/Desktop/flu_toolkit/flu_out/' + sort_bam_path + ' ' + scaffoldrigion + ' > ' + '/home/lei/Desktop/flu_toolkit/flu_out/' + file + '-' + scaffoldrigion + '.out.sam'
    os.system(samtoolscommand)
    print(samtoolscommand)
    

print(outfilelist)
sam2fas(outfilelist)
print('Done!!')
