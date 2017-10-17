# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:49:34 2017

@author: Ja
快速运行maff批量比对，在mafft.bat所在文件夹新建input文件夹，fas文件放里面。
输出后缀为out.fas在output文件夹
"""
import os
filelist = os.popen(r'dir /b input')
filelist = list(filelist)
command = "mafft --auto --reorder "

for file in filelist:
    file = file.strip()
    mafftcommand = command + '.\\input\\' + file +' > ' + '.\\output\\' + file +'.out.fas'
    os.system(mafftcommand)
    print(mafftcommand)

print('all %d has be done' % len(filelist))