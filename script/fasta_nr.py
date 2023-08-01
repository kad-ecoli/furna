#!/usr/bin/python2
docstring='''
fasta_nr.py input.fasta output.fasta output.tsv
    make de-duplicated fasta

Input:
    input.fasta - input fasta with duplicated sequence

Output:
    output.fasta - output fasta without duplicated sequence
    output.tsv   - list of name with the same sequence
'''

import sys, os

if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

inputFasta =sys.argv[1]
outputFasta=sys.argv[2]
outputTsv  =sys.argv[3]

fasta2header=dict()
nr_dict=dict()
header_list=[]
txt=''
fp=open(inputFasta)
for block in fp.read().lstrip('>').split('\n>'):
    if len(block.strip())==0:
        continue
    lines=block.splitlines()
    header=lines[0].split()[0]
    sequence=''.join(lines[1:])
    if sequence in fasta2header:
        rep=fasta2header[sequence]
        nr_dict[rep].append(header)
    else:
        fasta2header[sequence]=header
        nr_dict[header]=[header]
        header_list.append(header)
        txt+='>'+header+'\n'+sequence+'\n'
fp.close()

fp=open(outputFasta,'w')
fp.write(txt)
fp.close()

txt=''
for header in header_list:
    txt+='\t'.join(nr_dict[header])+'\n'
fp=open(outputTsv,'w')
fp.write(txt)
fp.close()
