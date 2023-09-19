#!/usr/bin/python3
docstring='''
clstr2tsv.py c80.fasta.clstr c80.tsv
'''
import sys

def clstr2dict(filename):
    fp=open(filename)
    blocks=('\n'+fp.read()).split("\n>")
    fp.close()
    clstr_dict=dict()
    clstr_list=[]
    for block in blocks[1:]:
        lines=block.splitlines()
        key=''
        values=[]
        for line in lines:
            if not ', >' in line:
                continue
            header=line.split('>')[1].split('.')[0]
            if line.rstrip()[-1]=='*':
                key=header
            values.append(header)
        clstr_dict[key]=values
        clstr_list.append(key)
    return clstr_list,clstr_dict


if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()
    clstr_list,clstr_dict=clstr2dict(sys.argv[1])
    txt=''
    for key in clstr_list:
        txt+='\t'.join(clstr_dict[key])+'\n'
    fp=open(sys.argv[2],'w')
    fp.write(txt)
    fp.close()

