#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import textwrap
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

#### read cgi parameters ####

form = cgi.FieldStorage()
page =form.getfirst("page",'').strip().strip("'")
if not page:
    page='1'
elif page=='0':
    page='1'
order=form.getfirst("order",'').lower().strip().strip("'")
if not order:
    order="pdbid"
pdbid=form.getfirst("pdbid",'').lower().strip().strip("'")
chain=form.getfirst("chain",'').strip().strip("'")
lig3 =form.getfirst("lig3",'').strip().strip("'")
if not lig3:
    lig3=form.getfirst("code",'').strip().strip("'")
rcl  =form.getfirst("rcl",'').upper().strip().strip("'")
rfm  =form.getfirst("rfm",'').upper().strip().strip("'")
ecn  =form.getfirst("ecn",'').strip().strip("'")
got  =form.getfirst("got",'').upper().strip().strip("'")
txn  =form.getfirst("txn",'').strip()
if got:
    got=got.split()[0]
    if got.startswith("GO:"):
        got=got[3:]
ligname=form.getfirst("ligname",'').upper().strip().strip('"').replace("'",'')
pubmed =form.getfirst("pubmed",'').strip("'")
outfmt =form.getfirst("outfmt",'').strip().strip("'")

para_list=[]
if order:
    para_list.append("order=%s"%order)
if pdbid:
    para_list.append("pdbid=%s"%pdbid)
if chain:
    para_list.append("chain=%s"%chain)
if lig3:
    para_list.append("lig3=%s"%lig3)
elif ligname:
    para_list.append('ligname="%s"'%ligname)
if rcl:
    para_list.append("rcl=%s"%rcl)
if rfm:
    para_list.append("rfm=%s"%rfm)
if txn:
    para_list.append("txn=%s"%txn)
if ecn:
    para_list.append("ecn=%s"%ecn)
if got:
    para_list.append("got=%s"%got)
if pubmed:
    para_list.append("pubmed=%s"%pubmed)
para='&'.join(para_list)

#### read database data ####
# pdb:recCha => L resolution rfam rnacentral pubmed ec go_mf go_bp go_cc taxon sequence cssr dssr title
fp=gzip.open(rootdir+"/data/rna.tsv.gz",'rt')
chain_dict=dict()
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    chain_dict[':'.join(items[:2])]=items[2:]
fp.close()

pubmed_dict=dict()
fp=gzip.open(rootdir+"/data/pubmed.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    pubmed_dict[items[0]]=items[1]
fp.close()

ligand_dict=dict()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    ccd  =items[0]
    name =items[-1]
    ligand_dict[ccd]=name
fp.close()
lig_set=set()
if lig3:
    lig_set=set([lig3])
    if lig3=="metal" or lig3=="regular":
        fp=gzip.open(rootdir+"/data/metal.tsv.gz",'rt')
        metal_set=set([line.split()[0] for line in fp.read().splitlines()])
        if lig3=="metal":
            lig_set=metal_set
        else:
            lig_set=set([ccd for ccd in ligand_dict if not ccd in metal_set])
elif ligname:
    lig_set=set([ccd for ccd in ligand_dict \
        if ligname in ligand_dict[ccd].replace("'",'').upper()])

if outfmt!='txt':
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
<img src=images/furna.png ></br>
<p><a href=.>[Back to Home]</a></p>
''')

enzyme_dict=dict()
fp=gzip.open(rootdir+"/data/enzyme.tsv.gz",'rt')
for line in fp.read().splitlines():
    e,name=line.split('\t')
    enzyme_dict[e]=name
fp.close()

go2name_dict=dict()
fp=gzip.open(rootdir+"/data/go2name.tsv.gz",'rt')
for line in fp.read().splitlines():
    g,a,name=line.split('\t')
    go2name_dict[g]='('+a+') '+name
fp.close()

#### parse page ####
pageLimit=200
if lig3 in ["peptide","rna","dna"]:
    pageLimit=100
html_txt=''
sort_line=[]
fp=gzip.open(rootdir+"/data/interaction.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
#for line in stdout.decode().splitlines():
    items=line.split('\t')
    pdb       =items[0]
    recCha    =items[1]
    assembly  =items[2]
    ccd       =items[3]
    ligCha    =items[4]
    ligIdx    =items[5]
    resOrig   =items[6]
    resRenu   =items[7]
    resSeq    =items[8]
    ligSeq    =items[9]
    
    target=pdb+':'+recCha
    if not target in chain_dict:
        continue
    
    if pdbid and pdb!=pdbid:
        continue
    if chain and recCha!=chain and ligCha!=chain:
        continue
    if len(lig_set) and not ccd in lig_set:
        continue

    items     =chain_dict[pdb+':'+recCha]
    L         =items[0]
    reso      =items[1]
    rfam      =items[2]
    rnacentral=items[3]
    pmid      =items[4]
    ec        =items[5]
    go_mf     =items[6].replace("GO:","")
    go_bp     =items[7].replace("GO:","")
    go_cc     =items[8].replace("GO:","")
    taxon     =items[9]
    sequence  =items[10]
    cssr      =items[11]
    dssr      =items[12]
    title     =items[13]

    if got:
        if got=='0' and not go_mf+go_bp+go_cc:
            continue
        elif got!='0' and not got in go_mf+go_bp+go_cc:
            continue
    if ecn:
        if ecn=='0' and not ecn:
            continue
        elif ecn!='0' and not ecn in ec:
            continue
    if rcl and not rcl in rnacentral:
        continue
    if rfm and not rfm in rfam:
        continue
    if txn and not txn in taxon.split(','):
        continue

    items=(pdb,recCha,assembly,ccd,ligCha,ligIdx,
        resOrig,resRenu,resSeq,ligSeq,
        L,reso,rfam,rnacentral,pmid,ec,go_mf,go_bp,go_cc,
        taxon,sequence,cssr,dssr,title)
    if outfmt=='txt':
        html_txt+='\t'.join(items)+'\n'
    else:
        if order=="reso":
            sort_line.append((reso,items))
        elif order=="rnacentral":
            sort_line.append((rnacentral,items))
        elif order=="lig3":
            sort_line.append((ccd,items))
        else:
            sort_line.append((pdb+recCha,items))
fp.close()

if outfmt=="txt":
    print("Content-type: text/plain\n")
    print(html_txt)
    exit()

sort_line.sort()
totalNum=len(sort_line)
totalPage=1+int(totalNum/pageLimit)
if not page:
    page=1
elif page=="last":
    page=totalPage
else:
    page=int(page)
if page<1:
    page=1
elif page>totalPage:
    page=totalPage

for l in range(totalNum):
    if l<pageLimit*(int(page)-1) or l>=pageLimit*(int(page)):
        continue
    items    =sort_line[l][1]

    pdb       =items[0]
    recCha    =items[1]
    assembly  =items[2]
    ccd       =items[3]
    ligCha    =items[4]
    ligIdx    =items[5]
    resOrig   =items[6]
    resRenu   =items[7]
    resSeq    =items[8]
    ligSeq    =items[9]
    L         =items[10]
    reso      =items[11]
    rfam      =items[12]
    rnacentral=items[13]
    pmid      =items[14]
    ec        =items[15]
    go_mf     =items[16]
    go_bp     =items[17]
    go_cc     =items[18]
    taxon     =items[19]
    sequence  =items[20]
    cssr      =items[21]
    dssr      =items[22]
    title     =items[23]
    
    if ec:
        ec_list=ec.split(',')
        ec=''
        for e in ec_list:
            if ec:
                ec+='<br>'
            if not e in enzyme_dict:
                ec+="<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(e,e)
            else:
                ec+='<a href=https://enzyme.expasy.org/EC/%s target=_blank><span title="%s">%s</span></a>'%(e,enzyme_dict[e],e)
    else:
        ec="N/A"
    go=''
    if go_mf+go_bp+go_cc:
        go_list=[]
        if go_mf:
            go_list+=["GO:"+g for g in go_mf.split(',')]
        if go_bp:
            go_list+=["GO:"+g for g in go_bp.split(',')]
        if go_cc:
            go_list+=["GO:"+g for g in go_cc.split(',')]

        go='<span title="'
        for g in go_list:
            if g in go2name_dict:
                go+=g+' '+go2name_dict[g]+'\n'
            else:
                go+=g+'\n'
        go=go[:-1]+'">'+go_list[0]+"...</span>"
    else:
        go="N/A"
    if rnacentral:
        rnacentral_list=[]
        for a in rnacentral.split(','):
            a=a.split('_')[0]
            a="<a href=https://rnacentral.org/rna/%s target=_blank>%s</a>"%(a,a)
            rnacentral_list.append(a)
        rnacentral='<br>'.join(rnacentral_list)
    else:
        rnacentral="N/A"
    if pmid:
        pubmed_list=[]
        for p in pmid.split(','):
            pubmed_list.append("<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(p,p))
            if p in pubmed_dict:
                pubmed_list[-1]='<span title="'+pubmed_dict[p]+ \
                    '">'+pubmed_list[-1]+"</span>"
        pmid='<br>'.join(pubmed_list)
    else:
        pmid="N/A"
    if taxon:
        taxon_list=[]
        for t in taxon.split(','):
            taxon_list.append("<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a>"%(t,t))
        taxon='<br>'.join(taxon_list)
    else:
        taxon="N/A"
    if rfam:
        rfam_list=[]
        for r in rfam.split(','):
            rfam_list.append("<a href=https://rfam.org/family/%s target=_blank>%s</a>"%(r,r))
        rfam='<br>'.join(rfam_list)
    else:
        rfam="N/A"

    name      =""
    ccd_http  =ccd
    reso="("+reso+")"
    if lig3 in ["protein","rna","dna"]:
        ccd_http='<pre>\n'+'\n'.join(textwrap.wrap(ligSeq,50))+'</pre>'
    else:
        if ccd in ligand_dict:
            name=';\n'.join(ligand_dict[ccd].split(';'))
            ccd_http='<span title="%s">%s</span>'%(name,ccd)
        elif ccd in ["protein","rna","dna"]:
            name='\n'.join(textwrap.wrap(ligSeq,50))
            ccd_http='<span title="%s">%s</span>'%(name,ccd)
        ccd_http='<a href="sym.cgi?code=%s" target=_blank>%s</a>'%(
            ccd,ccd_http)
    
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    html_txt+='''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><span title="%s"><a href="https://rcsb.org/structure/%s" target=_blank>%s:%s</a><br>%s</span></td>
    <td><span title="%s"><a href="pdb.cgi?pdb=%s&chain=%s&ligidx=%s&ligcha=%s" target=_blank>%s</span></td>
    <td style="word-wrap: break-word">%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    l+1,
    title,pdb,pdb,recCha,reso,
    resOrig,pdb,recCha,ligIdx,ligCha,ligCha+':'+ligIdx,
    ccd_http,
    go,
    rnacentral,
    rfam,
    taxon,
    pmid,
    )
fp.close()

print('''
Download all results in tab-seperated text for 
<a href="?outfmt=txt&%s" download="BioLiP.txt">%d receptor-ligand interactions</a>, whose format is explained at <a href="download/readme.txt">readme.txt</a>.<br>
<li>Hover over <strong>PDB</strong> to view the title of the structure.
Click <strong>PDB</strong> to view the structure at the RCSB PDB database.
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.</li>
<li>Click <strong>Site #</strong> to view the binding site structure.
Hover over <strong>Site #</strong> to view the binding residues.</li>
'''%(para,totalNum))
if lig3 in ["peptide","rna","dna"]:
    print("<li>The <strong>sequence</strong> is converted from residues with experimentally determined coordinates in the structure; residues not observed in the 3D structure are excluded.</li>")
else:
    print("<li>Hover over <strong>Ligand</strong> to view the full ligand name. Click <strong>Ligand</strong> to view the 2D diagram and other detail information of the ligand.</li>")
print('''
<li>Hover over <strong>GO terms</strong> to view all GO terms.
''')


print(('''<p></p>
<form name="sform" action="qsearch.cgi">
Sort results by
<select name="order" onchange="this.form.submit()">
    <option value="pdbid">PDB ID</option>
    <option value="lig3">Ligand ID</option>
    <option value="rcl">RNAcentral</option>
    <option value="reso">Resolution</option>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=rcl     value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=ligname value='%s'>
<input type=hidden name=pubmed  value='%s'>
</form>'''%(pdbid,lig3,rcl,ecn,got,ligname,pubmed)
).replace('value="%s"'%order,
          'value="%s" selected="selected"'%order))


print('''<center> 
<a class='hover' href='?&page=1&%s'>&lt&lt</a>
<a class='hover' href='?&page=%d&%s'>&lt</a>
'''%(para,page-1,para))
for p in range(page-10,page+11):
    if p<1 or p>totalPage:
        continue
    elif p==page:
        print(' %d '%(p))
    else:
        print('''<a class='hover' href='?&page=%d&%s'>%d</a>'''%(p,para,p))
print('''
<a class='hover' href='?&page=%d&%s'>&gt</a>
<a class='hover' href='?&page=last&%s'>&gt&gt</a>
<form name="pform" action="qsearch.cgi">Go to page <select name="page" onchange="this.form.submit()">
'''%(page+1,para,para))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print('''</select>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=rcl     value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=ligname value='%s'>
<input type=hidden name=pubmed  value='%s'>
</form></center><br>'''%(pdbid,lig3,rcl,ecn,got,ligname,pubmed))


if lig3 in ["protein","rna","dna"]:
    print('''
<style>
div.w {
  word-wrap: break-word;
}
</style>

<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th ALIGN=center><strong> # </strong></th>
    <th ALIGN=center><strong> PDB<br>(resolution) </strong></th>
    <th ALIGN=center><strong> Site<br># </strong></th>
    <th ALIGN=center><strong>''')
    print(lig3.upper() if lig3!="protein" else "Protein")
    print('''<br>sequence</strong> </th>           
    <th ALIGN=center><strong> GO<br>terms </strong> </th>           
    <th ALIGN=center><strong> RNAcentral </strong> </th>           
    <th ALIGN=center><strong> Rfam </strong> </th>           
    <th ALIGN=center><strong> Taxon </strong> </th>           
    <th ALIGN=center><strong> PubMed </strong> </th>           
</tr><tr ALIGN=center>
''')
else:
    print('''  
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th ALIGN=center><strong> # </strong></th>
    <th ALIGN=center><strong> PDB<br>(Resolution &#8491;) </strong></th>
    <th ALIGN=center><strong> Site # </strong></th>
    <th ALIGN=center><strong> Ligand </strong> </th>           
    <th ALIGN=center><strong> GO terms </strong> </th>           
    <th ALIGN=center><strong> RNAcentral </strong> </th>           
    <th ALIGN=center><strong> Rfam </strong> </th>           
    <th ALIGN=center><strong> Taxon </strong> </th>           
    <th ALIGN=center><strong> PubMed </strong> </th>           
</tr><tr ALIGN=center>
''')
print(html_txt)
print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
