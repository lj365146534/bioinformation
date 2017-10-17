#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

os.chdir('/home/lei/Desktop/flu_toolkit/flu_out/')
def sam2fas(filelist):
    
    for inline in filelist:
        #inline = line.strip()
        try:
            fastqseq = open(inline,'r')
            fastaseq = open(inline +'.fas','w')
        
            i = 0
            for line in fastqseq:
                i+=1
                seq = ">" + str(i) + "\n" + line.split("\t")[9] + "\n"
                fastaseq.write(seq)
            


            fastqseq.close()
            fastaseq.close()

            print(inline+'\t'+'转换成功！')
        except IOError:
            pass

if __name__ == '__main__':
    sam2fas(sys.argv[1:])
    print('全部完成')
