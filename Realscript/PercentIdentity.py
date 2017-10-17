#!/usr/bin/env python \
#用法python filename
#给入一行的fas文件，生成一个filename.PercentIdentity,计算序列同源性（一致性）
import sys

def OpenFileAsDict(filename):
#读取数据，设置序列名称和序列为键值，创建字典并返回
#缺点：readlines()，必须全部读取序列，大文件慢。

    fasta_file=open(filename,'r')
    
    readlines =fasta_file.readlines()
    t = 0
    seqdict ={}
    for i in range(len(readlines)):
        if readlines[i][0] == '>':
            seqdict[readlines[i].strip()] = readlines[i+1].strip()
    #print(seqdict)
    fasta_file.close()
    return seqdict




def Compare2Seq(seqdict):
#序列提取出来
#建立一个字典1储存seqname和seq，复制这个字典为字典2。for循环遍历字典1中seqa和字典2中seqb互相比较
#for i in dict1:
    #for j in dict1:
    Compareoutfile = open(outfile,'w')
    for seqkeya in seqdict:
        for seqkeyb in seqdict:
            i = 0
            k = 0
            if seqkeya != seqkeyb and len(seqdict[seqkeya]) == len(seqdict[seqkeyb]):
                for i in range(len(seqdict[seqkeyb])):
                    if seqdict[seqkeya][i] == seqdict[seqkeyb][i]:
                        i = i + 1
                    else:
                        i = i + 1
                        k = k + 1
                PercentIdentity = 1-k/len(seqdict[seqkeya])
                Compareoutfile.write('%s--%30s%20.4f\n' % (seqkeya,seqkeyb,PercentIdentity*100))
                #print(seqkeya,'--',seqkeyb,'\t',PercentIdentity*100,'n')
                #print('%s%30s%20.4f' % (seqkeya,seqkeyb,PercentIdentity*100),'\n')
            elif seqkeya == seqkeyb:
                pass
            else:
                print("Please alig these sequences: ",seqkeya,' or ',seqkeyb)


    Compareoutfile.close()


if __name__ == '__main__':
    filename = sys.argv[1]
    outfile = filename+'.PercentIdentity'
    #OpenFileAsDict(filename)
    Compare2Seq(OpenFileAsDict(filename))
    print('全部完成')
