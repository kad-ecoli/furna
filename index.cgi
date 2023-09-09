#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import subprocess

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT REFRESH START -->')[0]
    html_footer=txt.split('<!-- CONTENT REFRESH END -->')[-1]

print("Content-type: text/html\n")
if len(html_header):
   print(html_header)
else:
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>FURNA</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/furna.png></br>
<p><a href=.>[Back to Home]</a></p>
''')

filename=rootdir+"/data/index.txt"
if os.path.isfile(filename):
    fp=open(filename,'r')
    print(fp.read())
    fp.close()

cmd="ls output/*.cif.gz output/*.ent.gz|shuf|cut -f2 -d/|cut -f1 -d.|head -1"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode().strip()
if len(stdout):
    target=stdout.split('_')[0]
    pdbid=target[:4]
    chain=target[4:]
    if len(target)>=9 and target[:8]==target[:8].lower():
        pdbid=target[:8]
        chain=target[8:]

    if '_' in stdout: # interaction
        items=stdout.split('_')
        lig3=items[1]
        ligCha=items[2]
        ligIdx=items[3]
        print('''
<p>
<h1><span title="PDB $pdbid Chain $chain"><a href=pdb.cgi?pdb=$pdbid&chain=$chain&lig3=$lig3&ligCha=$ligCha&ligIdx=$ligIdx target=_blank>View an example ligand-RNA interaction</a></span></h1>
</p>
'''.replace("$pdbid",pdbid).replace("$chain",chain).replace("$ligCha",ligCha).replace("$ligIdx",ligIdx).replace("$lig3",lig3))
    else: # receptor
        print('''
<p>
<h1><span title="PDB $pdbid Chain $chain"><a href=pdb.cgi?pdb=$pdbid&chain=$chain target=_blank>View an example RNA</a></span></h1>
</p>
'''.replace("$pdbid",pdbid).replace("$chain",chain))


if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
