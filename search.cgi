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
if lig3:
    para_list.append("lig3=%s"%lig3)
elif ligname:
    para_list.append('ligname=%s'%ligname)
para='&'.join(para_list)

#### read database data ####
pubmed_dict=dict()
fp=gzip.open(rootdir+"/data/pubmed.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    pubmed_dict[items[0]]=items[1]
fp.close()

rfam_dict=dict()
fp=gzip.open(rootdir+"/data/rfam.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    rfam_dict[items[0]]=items[1]
fp.close()

taxon_dict=dict()
fp=gzip.open(rootdir+"/data/taxon.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    taxon_dict[items[0]]=items[1]
fp.close()

ligand_dict=dict()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    ccd  =items[0]
    name =items[5]
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

interaction_dict=dict()
cmd="zcat %s/data/interaction.tsv.gz|grep -v '^#'"%rootdir
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
for line in stdout.decode().splitlines():
    items=line.split('\t')
    key=':'.join(items[:2])
    ccd=items[3]
    if len(lig_set) and not ccd in lig_set:
        continue
    if not key in interaction_dict:
        interaction_dict[key]=[]
    interaction_dict[key].append(items[2:])
fp.close()


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
pageLimit=100
html_txt=''
sort_line=[]
fp=gzip.open(rootdir+"/data/rna.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
#for line in stdout.decode().splitlines():
    items=line.split('\t')
    pdb       =items[0]
    recCha    =items[1]
    L         =items[2]
    reso      =items[3]
    rfam      =items[4]
    rnacentral=items[5]
    pmid      =items[6]
    ec        =items[7].replace("EC:","")
    go_mf     =items[8].replace("GO:","")
    go_bp     =items[9].replace("GO:","")
    go_cc     =items[10].replace("GO:","")
    taxon     =items[11]
    sequence  =items[12]
    cssr      =items[13]
    dssr      =items[14]
    title     =items[15]

    if pdbid and pdb!=pdbid:
        continue
    if chain and recCha!=chain:
        continue

    if got:
        if got=='0' and not go_mf+go_bp+go_cc:
            continue
        elif got!='0' and not got in go_mf+go_bp+go_cc:
            continue
    if ecn:
        if ecn=='0' and not ec:
            continue
        elif ecn!='0' and not ecn in ec:
            continue
    if rcl and not rcl in rnacentral:
        continue
    if rfm and not rfm in rfam:
        continue
    if txn and not txn in taxon.split(','):
        continue
    if len(lig_set) and not pdb+':'+recCha in interaction_dict:
        continue

    items=(pdb,recCha,L,reso,rfam,rnacentral,pmid,ec,go_mf,go_bp,go_cc,
        taxon,sequence,cssr,dssr,title)
    if outfmt=='txt':
        html_txt+='\t'.join(items)+'\n'
    else:
        if order=="reso":
            sort_line.append((reso,items))
        elif order=="rnacentral":
            sort_line.append((rnacentral,items))
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
    L         =items[2]
    reso      =items[3]
    rfam      =items[4]
    rnacentral=items[5]
    pmid      =items[6]
    ec        =items[7]
    go_mf     =items[8].replace("GO:","")
    go_bp     =items[9].replace("GO:","")
    go_cc     =items[10].replace("GO:","")
    taxon     =items[11]
    sequence  =items[12]
    cssr      =items[13]
    dssr      =items[14]
    title     =items[15]
    
    if ec:
        ec_list=ec.split(',')
        ec=''
        for e in ec_list:
            if ec:
                ec+='<br>'
            if not e in enzyme_dict:
                ec+="<a href=https://enzyme.expasy.org/EC/%s target=_blank>EC:%s</a>"%(e,e)
            else:
                ec+='<a href=https://enzyme.expasy.org/EC/%s target=_blank><span title="%s">EC:%s</span></a>'%(e,enzyme_dict[e],e)
    go=''
    if go_mf+go_bp+go_cc:
        go_list=[]
        for g in (','.join((go_mf,go_bp,go_cc))).split(','):
            if not g:
                continue
            g="GO:"+g
            g='<span title="%s"><a href=https://www.ebi.ac.uk/QuickGO/term/%s target=_blank>%s</a></span>'%(go2name_dict[g],g,g)
            go_list+=[g]
        go='<br>'.join(go_list)
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
            taxon="<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a>"%(t,t)
            if t in taxon_dict:
                taxon='<i>'+'<br>'.join(taxon_dict[t].split())+'</i><br>'+taxon
            taxon_list.append(taxon)
        taxon='<br>'.join(taxon_list)
    else:
        taxon="N/A"
    if rfam:
        rfam_list=[]
        for r in rfam.split(','):
            rfam="<a href=https://rfam.org/family/%s target=_blank>%s</a>"%(r,r)
            if r in rfam_dict:
                rfam='<span title="'+rfam_dict[r]+'">'+rfam+"</span>"
            rfam_list.append(rfam)
        rfam='<br>'.join(rfam_list)
    else:
        rfam="N/A"

    name      =""
    ccd_http  =ccd
    if reso=="NA":
        reso=''
    else:
        reso="("+reso+"&#8491;)"
    ligand_list=[]
    key=pdb+':'+recCha
    if key in interaction_dict:
        interaction_list=interaction_dict[key]
        assembly=interaction_list[0][0]
        if assembly=='0':
            assembly=''
        else:
            assembly="Assembly"+assembly
        ligand_list.append(assembly)

        ccd_dict=dict()
        ccd_list=[]
        for items in interaction_list:
            ccd=items[1]
            if not ccd in ccd_dict:
                ccd_list.append(ccd)
                ccd_dict[ccd]=[]
            spantitle="Binding residues: "+items[4]
            if ccd in ["protein","rna","dna"]:
                spantitle=" sequence: "+items[-1]+'\n'+spantitle
                if ccd=="protein":
                    spantitle="Protein"+spantitle
                else:
                    spantitle=ccd.upper()+spantitle
            elif ccd in ligand_dict:
                spantitle=ligand_dict[ccd]+'\n'+spantitle
            ccd_dict[ccd].append('''
                <a href=pdb.cgi?pdbid=%s&chain=%s&lig3=%s&ligCha=%s&ligIdx=%s>
                <span title="%s">%s:%s:%s</span></a>
            '''%(pdb,recCha,ccd,items[2],items[3],
                spantitle,ccd,items[2],items[6].replace(' ','')))

        for ccd in ccd_list:
            ligand_list.append('<br>'.join(ccd_dict[ccd]))
    ccd_http='<br>'.join(ligand_list)

    ecgo=''
    if ec:
        ecgo=ec+'<br>'+go
    elif go:
        ecgo=go
    else:
        ecgo="N/A"
    
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    html_txt+='''
<tr %s>
    <td><a href=pdb.cgi?pdbid=%s&chain=%s>%d</a></td>
    <td style="word-wrap: break-word">
        <span title="%s"><a href="https://rcsb.org/structure/%s" target=_blank>%s</a>:%s<br>%s</span><br>
    </td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td style="word-wrap: break-word"><span title="Length=%s">%s</span></td>
</tr>
'''%(bgcolor,
    pdb,recCha,l+1,
    title,pdb,pdb,recCha,reso,
    ecgo,
    rnacentral,
    rfam,
    taxon,
    pmid,
    ccd_http,
    L,'<br>'.join(textwrap.wrap(sequence,50)),
    )
fp.close()

print('''
Download all results in tab-seperated text for 
<a href="?outfmt=txt&%s" download="database.txt">%d receptor-ligand interactions</a>.<br>
<li>Hover over <strong>PDB</strong> to view the title of the structure.</li>
<li>Hover over <strong>EC number &amp; GO term</strong> to view names of the Enzyme Commission numbers and Gene Ontology terms.</li>
<li>Hover over <strong>Rfam</strong> to view names of Rfam families.</li>
<li>Hover over <strong>PubMed</strong> to view title of the PubMed publications.</li>
<li>Hover over <strong>Ligand</strong> to view the name of the ligand and the ligand-binding nucleotides on the RNA.</li>
'''%(para,totalNum))


print(('''<p></p>
<form name="sform" action="search.cgi">
Sort results by
<select name="order" onchange="this.form.submit()">
    <option value="pdbid">PDB ID</option>
    <option value="lig3">Ligand ID</option>
    <option value="rcl">RNAcentral</option>
    <option value="reso">Resolution</option>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=chain   value='%s'>
<input type=hidden name=rcl     value='%s'>
<input type=hidden name=rfm     value='%s'>
<input type=hidden name=txn     value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=pubmed  value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=ligname value='%s'>
</form>'''%(pdbid,chain,rcl,rfm,txn,ecn,got,pubmed,lig3,ligname)
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
<form name="pform" action="search.cgi">Go to page <select name="page" onchange="this.form.submit()">
'''%(page+1,para,para))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print('''</select>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=chain   value='%s'>
<input type=hidden name=rcl     value='%s'>
<input type=hidden name=rfm     value='%s'>
<input type=hidden name=txn     value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=pubmed  value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=ligname value='%s'>
</form></center><br>'''%(pdbid,chain,rcl,rfm,txn,ecn,got,pubmed,lig3,ligname))


print('''  
<table border="0" width=100%>    
<tr BGCOLOR="#FF9900" align=left>
    <th><strong> # </strong></th>
    <th><strong> PDB<br>(Resolution)</strong></th>
    <th><strong> EC number<br>&amp; GO term </strong> </th>           
    <th><strong> RNAcentral </strong> </th>           
    <th><strong> Rfam </strong> </th>           
    <th><strong> Taxon </strong> </th>           
    <th><strong> PubMed </strong> </th>           
    <th><strong> Ligand </strong> </th>           
    <th><strong> RNA<br>sequence </strong> </th>           
</tr><tr>
''')
print(html_txt)
print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
