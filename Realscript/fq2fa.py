#!/usr/bin/env python 
import sys
def fq2fa(filelist):
    
    for fqname in filelist:    
        fastqseq = open(fqname,'r')
        fastaseq = open(fqname +'.fa','w')
        i = 0
        for line in fastqseq:
            i+=1
            if i % 4 == 1:
                fastaseq.write('>' + line)
            if i % 4 == 2:
                fastaseq.write(line)


        fastqseq.close()
        fastaseq.close()

        print(fqname+'\t'+'转换成功！')

if __name__ == '__main__':
    fq2fa(sys.argv[1:])
    print('全部完成')