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
    html_header=txt.split('<!-- CONTENT START -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]
if len(html_header)==0:
    html_header='''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>FURNA</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/furna.png></br>
<p><a href=.>[Back to Home]</a></p>
'''
if len(html_footer)==0:
    html_footer="</body> </html>"

print("Content-type: text/html\n")
print(html_header)
print('''
<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
  text-align: center;
}
</style>

<h4>Download tab-separated tables</h4>
<li>All RNAs: <a href=data/rna.tsv.gz download>rna.tsv.gz</a></li>
<li>All ligand-RNA interactions: <a href=data/interaction.tsv.gz download>interaction.tsv.gz</a></li>
<li>All DNAs that interact with RNAs: <a href=data/dna.tsv.gz download>dna.tsv.gz</a></li>
<li>All proteins that interact with RNAs: <a href=data/protein.tsv.gz download>protein.tsv.gz</a></li>
<li>All small molecule ligands (regular ligands and metal ions) that interact with RNAs: <a href=data/ligand.tsv.gz download>ligand.tsv.gz</a></li>
<li>Match to ATtRACT motifs for protein binding: <a href=data/attract_fimo.tsv.gz download>attract_fimo.tsv.gz</a></li>
<li>Full GO terms for each RNA, including the parent terms: <a href=data/parent.tsv.gz download>parent.tsv.gz</a></li>
<li>Documentation for the above tables: <a href=download/readme.md target=_blank>readme.md</a></li>

<p></p>

<h4>Download FASTA format sequences</h4>
<table width=100%>
<tr>
    <th>Molecule type</th>
    <th>All sequences</th>
    <th>Non-redundant sequences<br>(clustered at 100% sequence identity cutoff)</th>
    <th>Members of each cluster</th>
</tr>
<tr>
    <td>RNA (length &ge;10)</td>
    <td><a href=data/rna.fasta download>rna.fasta</a></td>
    <td><a href=data/rna_nr.fasta download>rna_nr.fasta</a></td>
    <td><a href=data/rna.tsv.gz download>rna.tsv.gz</a></td>
</tr>
<tr>
    <td>RNA fragment (length &le;9)</td>
    <td><a href=data/fragment.fasta download>fragment.fasta</a></td>
    <td> </td>
    <td> </td>
</tr>
<tr>
    <td>DNA</td>
    <td><a href=data/dna.fasta download>dna.fasta</a></td>
    <td><a href=data/dna_nr.fasta download>dna_nr.fasta</a></td>
    <td><a href=data/dna.tsv.gz download>dna.tsv.gz</a></td>
</tr>
<tr>
    <td>Protein</td>
    <td><a href=data/protein.fasta download>protein.fasta</a></td>
    <td><a href=data/protein_nr.fasta download>protein_nr.fasta</a></td>
    <td><a href=data/protein.tsv.gz download>protein.tsv.gz</a></td>
</tr>
</table>

<p></p>

<h4>Download PDB format structures</h4>
<li>All RNA chains: <a href=download/chain.tar download>chain.tar</a></li>
<li>All interaction partners (i.e., ligands) for each RNA chain: <a href=download/ligand.tar download>ligand.tar</a></li>

<p></p>

<h4>Download source code</h4>
<li>The source code for database curation and web interface display: <a href=https://github.com/kad-ecoli/furna target=_blank>https://github.com/kad-ecoli/furna</a></li>
''')
print(html_footer)
