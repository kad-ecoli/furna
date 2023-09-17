#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess
import textwrap

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT START -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]

def ExitWithError(msg,html_footer):
    print("ERROR!")
    print(msg)
    print("<p></p><a href=.>[Back]</a>")
    if len(html_footer):
        print(html_footer)
    else:
        print("</body> </html>")
    exit()

def splitid(target):
    pdbid=target[:4]
    chainid=target[4:]
    if len(target)>8 and target[:8]==target[:8].lower():
        pdbid=target[:8]
        chainid=target[8:]
    return pdbid,chainid

#### read cgi parameters ####

form = cgi.FieldStorage()
sequence=form.getfirst("sequence",'').strip()
if not sequence:
    sequence=form.getfirst("seq_file",'').strip()
seq_type =form.getfirst("seq_type",'').strip().lower()

print("Content-type: text/html\n")
if len(html_header):
    print(html_header)
else:
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/furna.png ></br>
<p><a href=.>[Back to Home]</a></p>
''')
if not seq_type in ["protein","rna","dna"]:
    ExitWithError("Sequence type must be one of the following: protein, rna, dna, peptide",html_footer)
header=seq_type
txt=''
for line in sequence.splitlines():
    line=line.strip()
    if line.startswith('>'):
        if header[0]=='>':
            print("ERROR! only one sequence allowed per search")
            exit()
        else:
            header=line
    else:
        txt+=line.upper()
if header[0]=='>':
    header=header[1:]
if seq_type in ["rna","dna"]:
    sequence=txt.lower()
else:
    sequence=txt.upper()
if len(set(txt).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
    ExitWithError("Unknown residue type "+' '.join(set(txt
        ).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))),html_footer)
if len(sequence)>1500:
    print("ERROR! Unable to handle sequence with %d &gt; 1500 residues. Only parse the first 1500 residues."%len(sequence))
sequence=sequence[:1500]
print('&gt;'+header+'<br>')
print('<br>'.join(textwrap.wrap(sequence,80))+'<br><br>')
print("The query sequence (length=%d) is searched through a non-redundant set of database sequences <a href=data/%s_nr.fasta>%s_nr.fasta</a> clustered at 100%% identity cutoff to identify representative hits. Homologs that belong to the same sequence cluster of the representative hit are listed in the last column of the table.<p></p>"%(len(sequence),seq_type,seq_type))

blast="blastp"
if seq_type.endswith("na"):
    blast="blastn"
cmd="echo %s|%s/script/%s -db %s/data/%s_nr.fasta -max_target_seqs 1000 -outfmt '6 sacc slen evalue nident length' "%(sequence,rootdir,blast,rootdir,seq_type)
score_name="E-value"
if len(sequence)<30:
    cmd="echo %s | %s/script/NWalign - %s/data/%s_nr.fasta|sort -k3nr|head -1000"%(sequence,rootdir,rootdir,seq_type)
    score_name="Alignment<br>score"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
lines=stdout.decode().splitlines()
print('''
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> Hit </strong></th>
    <th width=10% ALIGN=center><strong> Hit<br>length </strong></th>
    <th width=10% ALIGN=center><strong> Aligned<br>length </strong></th>
    <th width=10% ALIGN=center><strong> Identity<br>(normalized<br>by query)</strong> </th>           
    <th width=10% ALIGN=center><strong> Identity<br>(normalized<br>by hit)</strong> </th>           
    <th width=10% ALIGN=center><strong> Identity (normalized by<br>aligned length)</strong> </th>           
    <th width=10% ALIGN=center><strong> '''+score_name+ '''</strong> </th>           
    <th width=25% ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
''')
hit2clust_dict=dict()
if len(lines):
    fp=open("%s/data/%s_nr.tsv"%(rootdir,seq_type),'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        if len(items)>1:
            hit2clust_dict[items[0]]=items[1:]
    fp.close()

totalNum=0
sacc_list=[]
for line in lines:
    sacc,slen,evalue,nident,Lali=line.split('\t')
    if sacc in sacc_list:
        continue
    totalNum+=1
    bgcolor=''
    if totalNum%2==0:
        bgcolor='BGCOLOR="#DEDEDE"'
    slen=int(slen)
    nident=float(nident)
    Lali=int(Lali)

    
    pdbid,chainid=splitid(sacc)
    hit="<a href=search.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
        pdbid,chainid,pdbid,chainid)
    homolog_list=[]
    if sacc in hit2clust_dict:
        for mem in hit2clust_dict[sacc]:
            pdbid,chainid=splitid(mem)
            homolog_list.append("<a href=search.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
                pdbid,chainid,pdbid,chainid))
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td>%s</td>
    <td>%d</td>
    <td>%d</td>
    <td>%.4f</td>
    <td>%.4f</td>
    <td>%.4f</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    totalNum,
    hit,
    slen,
    Lali,
    nident/len(sequence),
    nident/slen, 
    nident/Lali,
    evalue, 
    ', '.join(homolog_list)))
print("</table><p></p><a href=.>[Back]</a>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
