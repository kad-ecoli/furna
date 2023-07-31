#!/usr/bin/python2
docstring=''' 
SortFastaWithResolution.py resolu.idx pdb_atom.fasta pdb_atom.sorted.fasta
    Sort fasta by length. If two sequence has the same length, take the one
    with better (lower) resolution. Output the files to pdb_atom.sorted.fasta
'''

import sys
if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

resolufile=sys.argv[1]
fastafile =sys.argv[2]
outfile   =sys.argv[3]

resolu_dict=dict()
fp=open(resolufile,'r')
lines=fp.read().splitlines()
fp.close()
for line in lines:
    if not "\t;\t" in line:
        continue
    items=line.split("\t;\t")
    if len(items)!=2:
        continue
    idcode=items[0].lower()
    resolu="NA"
    if len(items[1]):
        if float(items[1])>0:
            resolu=float(items[1])
    resolu_dict[idcode]=resolu

fasta_list=[]
fp=open(fastafile,'r')
blocks=fp.read().strip().lstrip('>').split('\n>')
fp.close()
for block in blocks:
    lines=block.splitlines()
    header=lines[0]
    sequence=''.join(lines[1:])
    idcode=header.split(':')[0].lower()
    resolu="NA"
    if idcode in resolu_dict:
        resolu=resolu_dict[idcode]
    L=len(sequence)
    fasta_list.append((resolu,L,header.split('\t')[0].replace(':',''),sequence))

fasta_list=sorted(fasta_list)
fasta_list=[(-L,resolu,header,sequence) for resolu,L,header,sequence in fasta_list]
fasta_list=sorted(fasta_list)
fasta_list=[(-L,resolu,header,sequence) for L,resolu,header,sequence in fasta_list]
sequence_list=set()

txt=''
for L,resolu,header,sequence in fasta_list:
    txt+=">%s\n%s\n"%(header,sequence)
fp=open(outfile,'w')
fp.write(txt)
fp.close()
