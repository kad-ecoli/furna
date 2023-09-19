#!/usr/bin/python3

import argparse
import os
from string import Template

#import pandas as pd

#from gscripts import pwm

class MyTemplate(Template):
    delimiter = '&'

parser = argparse.ArgumentParser(description="""takes a cisbp or rbpdb formatted file and converts it to a meme / fimo formatted pwm""")
 
parser.add_argument("--pwm_file", "-p", help="cisbp/rbpdb pwm file", required=True)
parser.add_argument("--out_file", "-o", help="output file", required=True)
parser.add_argument("--freq_file", "-f", help="frequency file", required=False)

args = parser.parse_args()

freq_txt="A 0.250 C 0.250 G 0.250 T 0.250"

template = MyTemplate('''MEME version 4

ALPHABET= ACGT
strands: +

Background letter frequencies
&freq_txt

MOTIF &name
letter-probability matrix: alength= 4 w= &width nsites= &nsites
''')

pwm_name = ".".join(os.path.basename(args.pwm_file).split(".")[:-1])

#x = pd.read_csv(args.pwm_file, sep="\t", index_col=0)
x=[]
fin = open(args.pwm_file, 'r')
for line in fin.read().splitlines()[1:]:
    x.append(line.split('\t')[1:])
fin.close()

if args.freq_file:
    fin = open(args.freq_file, 'r')
    freq_txt=fin.read().strip()
    fin.close()

result = template.substitute(
    name=pwm_name, freq_txt=freq_txt, width=len(x), nsites=len(x))
    
for row in x:
    total=sum([float(e) for e in row])
    if total<1:
        result += "  ".join(["{:.6g}".format(float(e)/total) for e in row]) + "\n"
    else:
        result += "  ".join(["{:.6g}".format(float(e)) for e in row]) + "\n"
    #result += "  ".join(row[1].map(lambda x: "{:.6}".format(x))) + "\n"
    
outfile=open(args.out_file, 'w')
outfile.writelines(result)
outfile.close()
