#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import subprocess
import re

rootdir=os.path.dirname(os.path.abspath(__file__))
bindir=rootdir+"/script"

def ExitWithError(msg):
    print("<br>ERROR!")
    print(msg)
    print("<p></p><a href=.>[Back]</a></body> </html>")
    exit()

print('''Content-type: text/html

<html>
<head>
<title>Structure search</title>
</head>
''')

form = cgi.FieldStorage()
mode = 1
structure=form.getfirst("structure",'').strip()
if not structure:
    structure=form.getfirst("struct_file",'').decode().strip()
if not structure:
    ExitWithError("empty input")
if len(structure) and not "ATOM " in structure:
    if re.match("^\w+:\w+$", structure):
        mode=3
    elif re.match("^[A-Z0-9]+$", structure):
        mode=2
    else:
        ExitWithError("unknown input format")

outdir=rootdir
cmd="date +%N"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
jobID=stdout.decode().strip()
if os.getenv("REMOTE_ADDR"):
    jobID+='.'+os.getenv("REMOTE_ADDR").replace(':','.')
jobID+=".fsearch"
outdir=rootdir+'/output/'+jobID

#### download pdb if necessary ####
if "\nloop_" in structure:
    fp=open(outdir+".cif",'w')
else:
    fp=open(outdir+".pdb",'w')
fp.write(structure+'\n')
fp.close()

filename=jobID+".cif"
cmd="%s/pdb2fasta -ter 2 -split 0 %s.cif"%(bindir,outdir)
if os.path.isfile(outdir+".pdb"):
    cmd="%s/pdb2fasta  -ter 2 -split 0 %s.pdb"%(bindir,outdir)
    filename=jobID+".pdb"
cmd+="|grep -v '>'"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
sequence=stdout.decode().strip()
if len(sequence)<5:
    ExitWithError("Sequence too short: L = %d &lt; 5<br>%s"%(
        len(sequence),sequence))
elif len(sequence)>1500:
    ExitWithError("Sequence too long: L = %d &gt; 1500<br>%s"%(
        len(sequence),sequence))

fp=open(outdir+".html",'w')
fp.write('''<html>
<head>
<title>Structure search</title>
</head>
<body>
&gt;<a href=%s>%s</a> (L=%d)<br>
%s
<p></p>
The structure search will take a few minutes.
This page will be refreshed every 20 seconds.
You may bookmark <a href=%s.html>this page</a> and return later.
<meta http-equiv="refresh" content="20; url='%s.html'"/>
</body>
</html>
'''%(filename,filename,len(sequence),sequence,jobID,jobID))
fp.close()

print('''<meta http-equiv="refresh" content="0; url='output/%s.html'"/>
</body> </html>'''%jobID)
