#!/usr/bin/env python 
fasta_file=open('pickcat.fas','r')
out_file = open('delambigous.derepeat.pickcat.fas','w')
seq=''
UniqueSeq=[]
nucleotide={'A','G','T','C','a','t','c','g'}
smallnucleotide={}
#UniqueHeader=[]
for line in fasta_file:
    line = line.strip()
    if line[0] == '>' and seq == '':
        # process the first line of the input file
        header = line
    elif line [0] != '>' and nucleotide|set(line)==nucleotide:
	    #delete set(line) not in  ATGCagtc 
        # join the lines with sequence
		seq = line
    elif line[0] == '>' and seq != '':
        # in subsequent lines starting with '>',
        # write the previous header and sequence
        # to the output file. Then re-initialize
        # the header and seq variables for the next record
        if seq not in UniqueSeq:   # and header not in UniqueHeader:
		    UniqueSeq.append(seq) # and UniqueHeader.append(header)
		    out_file.write(header+'\n' + seq+'\n')
        seq = ''
        #header = line
fasta_file.close()
out_file.close()

