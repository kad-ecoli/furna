#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess
import textwrap
import tarfile
import gzip
import re

rootdir=os.path.dirname(os.path.abspath(__file__))

def get_svg_ratio(svgfile):
    fp=open(svgfile)
    for line in fp:
        if line.startswith('<svg width='):
            items=line.split('"')
            width=float(items[1].replace("pt",''))
            height=float(items[3].replace("pt",''))
            return width/height
    fp.close()
    return 0

def display_receptor(rna_info_list,taxon_dict,parent_dict,ec_dict,go_dict,rfam_dict):
    pdbid     =rna_info_list[0]
    asym_id   =rna_info_list[1]
    L     =int(rna_info_list[2])
    reso      =rna_info_list[3]
    rfam_list =rna_info_list[4].split(',')
    ec_list   =rna_info_list[7].strip().split(',')
    mf_list   =rna_info_list[8].split(',')
    bp_list   =rna_info_list[9].split(',')
    cc_list   =rna_info_list[10].split(',')
    taxon_list=rna_info_list[11].split(',')
    sequence  =rna_info_list[12]
    cssr      =rna_info_list[13]
    dssr      =rna_info_list[14]
    title     =rna_info_list[15]

    mf_parent =''
    bp_parent =''
    cc_parent =''
    target=pdbid+':'+asym_id
    if target in parent_dict:
        mf_parent=parent_dict[target][0]
        bp_parent=parent_dict[target][1]
        cc_parent=parent_dict[target][2]

    rfam_table=''
    if len(rfam_list):
        rfam_table+="<tr><td align=center><strong>Rfam<br>families</strong></td><td>"
        for r,rfam in enumerate(rfam_list):
            if r:
                rfam_table+="<br>"
            rfam_table+="<a href=https://rfam.org/family/%s target=_blank>%s</a> "%(rfam,rfam)
            if rfam in rfam_dict:
                rfam_table+=rfam_dict[rfam]
        rfam_table+="</td></tr>"

    go_table=''
    if mf_parent+bp_parent+cc_parent:
        go_table+='''<tr BGCOLOR="#DEDEDE"><td align=center><strong>GO<br>terms</strong></td><td>
        <table><tr><td><table>'''
        for go in mf_list+bp_list+cc_list:
            go_table+="<tr><td><a href=https://www.ebi.ac.uk/QuickGO/term/%s target=_blank>%s</a>"%(go,go)
            if go in go_dict:
                go_table+=" (%s) %s"%(go_dict[go][0],go_dict[go][1])
            go_table+="</td></tr>\n"
        go_table+="</table><tr><td><table><tr valign=bottom>"

        parent2idx=dict()
        fp=open("%s/data/gosvg/list"%rootdir)
        for i,line in enumerate(fp.read().splitlines()):
            parent2idx[line]=i+1
        fp.close()

        width_list=[]
        if mf_parent:
            mf_svg="data/gosvg/%d.svg"%parent2idx[mf_parent]
            width_list.append(get_svg_ratio(mf_svg))
        if bp_parent:
            bp_svg="data/gosvg/%d.svg"%parent2idx[bp_parent]
            width_list.append(get_svg_ratio(bp_svg))
        if cc_parent:
            cc_svg="data/gosvg/%d.svg"%parent2idx[cc_parent]
            width_list.append(get_svg_ratio(cc_svg))
        if len(width_list)>=2:
            sumwidth=sum(width_list)
            width_list=[int(100.*w/sumwidth) for w in width_list]
        else:
            width_list=[]
            
        width_idx=0
        width=''
        if mf_parent:
            if len(width_list):
                width="width=%d%%"%width_list[width_idx]
                width_idx+=1
            go_table+='<td align=center %s><a href="%s" target=_blank><img src="%s" style="display:block;" width="100%%"><br>View graph for<br>Molecular Function (F)</a></td>'%(width,mf_svg,mf_svg)
        if bp_parent:
            if len(width_list):
                width="width=%d%%"%width_list[width_idx]
                width_idx+=1
            go_table+='<td align=center %s><a href="%s" target=_blank><img src="%s" style="display:block;" width="100%%"><br>View graph for<br>Biological Process (P)</a></td>'%(width,bp_svg,bp_svg)
        
        if cc_parent:
            if len(width_list):
                width="width=%d%%"%width_list[width_idx]
                width_idx+=1
            go_table+='<td align=center %s><a href="%s" target=_blank><img src="%s" style="display:block;" width="100%%"><br>View graph for<br>Cellular Component (C)</a></td>'%(width,cc_svg,cc_svg)

        
        go_table+="</tr></table></td></tr></table></td></tr>"

    ec_table=''
    if len(ec_list) and ec_list[0]:
        ec_table+='<tr><td align=center><strong>EC<br>numbers</strong></td><td>'
        for e,ec in enumerate(ec_list):
            if e:
                ec_table+="<br>"
            ec=ec.replace('EC:','')
            ec_table+="<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(ec,ec)
            if ec in ec_dict:
                ec_table+=" "+ec_dict[ec]
            
        ec_table+="</td></tr>"
        

    prefix=pdbid+asym_id

    if reso=="-1.00" or reso=="NA":
        reso="N/A"
    else:
        reso=reso+" &#8491;"
    divided=pdbid[-3:-1]
    chainfile="chain/%s/%s.pdb.gz"%(divided,prefix)
    arenafile="arena/%s/%s.pdb.gz"%(divided,prefix)

    seq_txt=''
    width=100
    for i in range(0,L,width):
        seq_txt+=sequence[i:(i+width)]+' - %d<br>'%(i+width if i+width<L else L)
        seq_txt+='<span title="CSSR secondary structure assignment">'+cssr[i:(i+width)]+'</span><br>'
        seq_txt+='<span title="DSSR secondary structure assignment">'+dssr[i:(i+width)]+'</span><br>'

    taxon_txt=''
    for t,taxon in enumerate(taxon_list):
        if t:
            taxon_txt+="<br>"
        if not taxon in taxon_dict:
            taxon_dict[taxon]=''
        taxon_txt+="<i>%s</i> (<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a>)"%(
            taxon_dict[taxon],taxon,taxon)

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Receptor RNA</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>Sequence &amp;<br>secondary<br>structure</strong></td><td>$seq_txt</td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><strong>PDB</strong></td><td><span title="Search other entries from the same structure"><a href=search.cgi?pdbid=$pdbid target=_blank>$pdbid</a></span> $title</td></tr>
    <tr><td align=center><strong>Chain</strong></td><td><span title="Search other entries for the same chain"><a href=search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$asym_id</a></span></td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><strong>Resolution</strong></td><td>$reso</a></td></tr>
    <tr align=left><td align=center><strong>3D<br>structure</strong></td><td>
    <table><tr align=left>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "load $chainfile; color background black; cartoons; color group; spacefill off; wireframe off;"
                }
                $("#mychain").html(Jmol.getAppletHtml("jmolApplet0",Info))
            });
            </script>
            <span id=mychain></span>
        </td>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "load $arenafile; color background black; cartoons; color group; spacefill off; wireframe off;"
                }
                $("#myarena").html(Jmol.getAppletHtml("jmolApplet1",Info))
            });
            </script>
            <span id=myarena></span>
        </td>
    </tr>
    <tr align=left>
        <td>
            Original structue of $pdbid$asym_id<br>
            [<a href="javascript:Jmol.script(jmolApplet0, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet0, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet0, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet0, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet0, 'color background black')">Black background</a>]<br>
            [<a href=$chainfile>Download</a>]<br>
        </td>
        <td>
            $pdbid$asym_id with missing atoms filled by <a href=https://zhanggroup.org/Arena>Arena</a><br>
            [<a href="javascript:Jmol.script(jmolApplet1, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet1, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet1, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet1, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet1, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet1, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet1, 'color background black')">Black background</a>]<br>
            [<a href=$arenafile>Download</a>]<br>
        </td>
    </tr></table>


    </td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><strong>Species</strong></td><td>$taxon_txt</td></tr>
    $rfam_table
    $go_table
    $ec_table
    </table>
</div>
</td></tr>
'''.replace("$pdbid"    ,pdbid
  ).replace("$title"    ,title
  ).replace("$asym_id"  ,asym_id
  ).replace("$seq_txt"  ,seq_txt
  ).replace("$reso"     ,reso
  ).replace("$chainfile",chainfile
  ).replace("$arenafile",arenafile
  ).replace("$taxon_txt",taxon_txt
  ).replace("$rfam_table",rfam_table
  ).replace("$go_table",go_table
  ).replace("$ec_table",ec_table
  ))

    #if go:
        #display_go(go,uniprot,pdbid,asym_id)
    
    #if ec:
        #display_ec(ec,csaOrig,csaRenu)
    return

def sanitize(inString):
    alphabet=set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890-")
    return ''.join([s for s in inString if s in alphabet])

def get_rna_info(pdbid,asym_id):
    rna_info_list=[]
    fp=gzip.open("%s/data/rna.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        if not line.startswith(pdbid):
            continue
        items=line.split('\t')
        if items[0]==pdbid and items[1]==asym_id:
            rna_info_list=items
    fp.close()
    return rna_info_list

def get_ligand_info(pdbid,asym_id,lig3,ligCha,ligIdx):
    ligand_info_list=[]
    fp=gzip.open("%s/data/interaction.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        if not line.startswith(pdbid):
            continue
        items=line.split('\t')
        if items[0]!=pdbid or items[1]!=asym_id or items[3]!=lig3 and items[4]!=ligCha:
            continue
        if items[5]==ligIdx or lig3 in ["protein","rna","dna"]:
            ligand_info_list=items
    fp.close()
    return ligand_info_list

def read_taxon():
    taxon_dict=dict()
    fp=gzip.open("%s/data/taxon.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        taxon,name=line.split('\t')
        taxon_dict[taxon]=name
    fp.close()
    return taxon_dict

def read_parent():
    parent_dict=dict()
    fp=gzip.open("%s/data/parent.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        parent_dict[':'.join(items[:2])]=items[2:]
    fp.close()
    return parent_dict

def read_go():
    go_dict=dict()
    fp=gzip.open("%s/data/go2name.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        go_dict[items[0]]=items[1:]
    fp.close()
    return go_dict

def read_ec():
    ec_dict=dict()
    fp=gzip.open("%s/data/enzyme.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        ec_dict[items[0]]=items[1]
    fp.close()
    return ec_dict

def read_rfam():
    rfam_dict=dict()
    fp=gzip.open("%s/data/rfam.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        rfam_dict[items[0]]=items[1]
    fp.close()
    return rfam_dict

if __name__=="__main__":
    form   =cgi.FieldStorage()
    pdbid  =form.getfirst("pdbid",'').lower()
    if not pdbid:
        pdbid=form.getfirst("pdb",'').lower()
    asym_id=form.getfirst("chain",'')
    if not asym_id:
        asym_id=form.getfirst("asym_id",'')
    ligIdx =form.getfirst("idx",'')
    if not ligIdx:
        ligIdx =form.getfirst("ligIdx",'')
    ligCha =form.getfirst("ligCha",'')
    lig3   =form.getfirst("lig3",'')
    outfmt =form.getfirst("outfmt",'')

    pdbid  =sanitize(pdbid)
    asym_id=sanitize(asym_id)
    ligIdx =sanitize(ligIdx)
    ligCha =sanitize(ligCha)
    lig3   =sanitize(lig3)
    outfmt =sanitize(outfmt)

    if outfmt=='1':
        download_pdb1(pdbid,asym_id,lig3,ligCha,ligIdx)
        exit(0)
        
    print("Content-type: text/html\n")
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>FURNA</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/furna.png ></br>

<script type="text/javascript" src="jsmol/JSmol.min.js"></script>
<table style="table-layout:fixed;" width="100%" cellpadding="2" cellspacing="0">
<table width=100%>
''')
    fp=open("%s/index.html"%rootdir)
    txt=fp.read()
    if '<!-- MENU START -->' in txt and '<!-- MENU END -->' in txt:
        print(txt.split('<!-- MENU START -->')[1].split('<!-- MENU END -->')[0])
    fp.close()

    rna_info_list=[]
    if pdbid and asym_id:
        rna_info_list=get_rna_info(pdbid,asym_id)
    if len(rna_info_list)==0:
        print('''</table>
Unknown pdb chain %s:%s
<p></p>
[<a href=.>Back to HOME</a>]
</body> </html>'''%(pdbid,asym_id))
        exit()
    
    title="PDB "+pdbid+" Chain "+asym_id
    if (lig3 and ligCha and ligIdx):
        if lig3 in ["protein","dna","rna"]:
            title+=" interacting with "+lig3+" in Chain "+ligCha
        else:
            title+=" interacting with "+lig3+" "+ligIdx
            if asym_id!=ligCha:
                title+=" in Chain "+ligCha
    print("<tr><td><h1 align=center>"+title+"</h1></td></tr>")

    taxon_dict=read_taxon()
    parent_dict=read_parent()
    go_dict=read_go()
    ec_dict=read_ec()
    rfam_dict=read_rfam()

    display_receptor(rna_info_list,taxon_dict,parent_dict,ec_dict,go_dict,rfam_dict)
    if lig3 and ligCha and ligIdx:
        ligand_info_list=get_ligand_info(pdbid,asym_id,lig3,ligCha,ligIdx)
        #display_interaction(rna_info_list)
   
    rnacentral=''
    if len(rna_info_list[5]):
        rnacentral_list=[]
        for r in rna_info_list[5].split(','):
            rnacentral_list.append("<a href=https://rnacentral.org/rna/%s target=_blank>%s</a>"%(r,r))
        rnacentral+='''<tr BGCOLOR="#DEDEDE"><td align=center><strong>RNAcentral</strong></td><td>%s</td></tr>'''%(''.join(rnacentral_list))
    
    pubmed=''
    if len(rna_info_list[6]):
        pubmed_dict=dict()
        fp=gzip.open("%s/data/pubmed.tsv.gz"%rootdir,'rt')
        for line in fp.read().splitlines():
            items=line.split('\t')
            pubmed_dict[items[0]]=items[1]
        fp.close()
        pubmed_list=[]
        for p in rna_info_list[6].split(','):
            pubmed="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(p,p)
            if p in pubmed_dict:
                pubmed=pubmed_dict[p]+' '+pubmed
            pubmed_list.append("<li>"+pubmed+"</li>")
        if rnacentral:
            pubmed='''<tr>'''
        else:
            pubmed='''<tr BGCOLOR="#DEDEDE">'''
        pubmed+='''<td align=center><strong>PubMed</strong></td><td>%s</td></tr>'''%(''.join(pubmed_list))
            
            
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">External links</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>PDB</strong></td>
        <td width=90%>
            <a href=https://rcsb.org/structure/$pdbid target=_blank>RCSB</a>,
            <a href=https://ebi.ac.uk/pdbe/entry/pdb/$pdbid target=_blank>PDBe</a>,
            <a href=https://pdbj.org/mine/summary/$pdbid target=_blank>PDBj</a>,
            <a href=http://ebi.ac.uk/pdbsum/$pdbid target=_blank>PDBsum</a>,
            <a href=https://nakb.org/atlas=$pdbid target=_blank>NAKB</a>,
            <a href=https://dnatco.datmos.org/conformers_cif.php?cifcode=$pdbid target=_blank>DNATCO</a>,
            <a href=http://rna.bgsu.edu/rna3dhub/pdb/$pdbid target=_blank>BGSU RNA</a>
        </td>
    </tr>
    $rnacentral
    $pubmed
    </tr>
'''.replace("$pdbid",pdbid
  ).replace("$pubmed",pubmed
  ).replace("$rnacentral",rnacentral
  ))
    
    print('''</table>
<p></p>
</body> </html>''')
