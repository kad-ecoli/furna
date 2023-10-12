#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import subprocess
import re
import textwrap

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
<title>Infernal search</title>
</head>
''')

form = cgi.FieldStorage()
sequence=form.getfirst("sequence",'').strip()
if not sequence:
    sequence=form.getfirst("seq_file",'').strip()
txt=''
header=''
for line in sequence.splitlines():
    line=line.strip()
    if line.startswith('>'):
        if header:
            print("ERROR! only one sequence allowed per search")
            exit()
        else:
            header=line
    else:
        txt+=line.upper()
if len(header)==0:
    header=">infernal"
if header[0]=='>':
    header=header[1:]
sequence=txt.upper()
if len(set(txt).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
    ExitWithError("Unknown residue type "+' '.join(set(txt
        ).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))),html_footer)
sequence=sequence.replace('T','U')
sequence=sequence[:1500]

outdir=rootdir
cmd="date +%N"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
jobID=stdout.decode().strip()
if os.getenv("REMOTE_ADDR"):
    jobID+='.'+os.getenv("REMOTE_ADDR").replace(':','.')
jobID+=".isearch"
outdir='output/'+jobID

filename=outdir+".fasta"
fp=open(filename,'w')
fp.write(">"+header+'\n'+sequence+'\n')
fp.close()

if len(sequence)<5:
    ExitWithError("Sequence too short: L = %d &lt; 5<br>%s"%(
        len(sequence),sequence))
elif len(sequence)>1500:
    ExitWithError("Sequence too long: L = %d &gt; 1500<br>%s"%(
        len(sequence),sequence))

fp=open(outdir+".html",'w')
fp.write('''<html>
<head>
<title>Infernal search</title>
</head>
<body>
&gt;<a href=%s.fasta>%s.fasta</a> (L=%d)<br>
%s
<p></p>
The infernal search will take a few minutes.
This page will be refreshed every 20 seconds.
You may bookmark <a href=%s.html>this page</a> and return later.
<meta http-equiv="refresh" content="20; url='%s.html'"/>
</body>
</html>
'''%(jobID,jobID,len(sequence),sequence,jobID,jobID))
fp.close()

print('''<meta http-equiv="refresh" content="0; url='output/%s.html'"/>
</body> </html>'''%jobID)
