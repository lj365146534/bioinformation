#!/usr/bin/env python 

a = ['']
m = 5
i = 0
seq = a * (m+1)
interleaved_format_phy = open('test8.phy','r')
sequential_format_phy = open('test8.trans.phy','w')
for line in interleaved_format_phy:
    if i == 0:
        i += 1
    elif 0 < i <= m:
        seq[i] = seq[i] + line
        i += 1
    elif i == m:
        seq[i] = seq[i] + line
        i = 0
sequential_format_phy.write(file = seq)
close(interleaved_format_phy)
close(sequential_format_phy)