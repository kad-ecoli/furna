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

def display_receptor(rna_info_list):
    pdbid     =rna_info_list[0]
    asym_id   =rna_info_list[1]
    L     =int(rna_info_list[2])
    reso      =rna_info_list[3]
    rfam_list =rna_info_list[4].split(',')
    ec_list   =rna_info_list[7].split(',')
    mf_list   =rna_info_list[8].split(',')
    bp_list   =rna_info_list[9].split(',')
    cc_list   =rna_info_list[10].split(',')
    taxon_list=rna_info_list[11].split(',')
    sequence  =rna_info_list[12]
    cssr      =rna_info_list[13]
    dssr      =rna_info_list[14]
    title     =rna_info_list[15]

    prefix=pdbid+asym_id
    species="Species: "

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
  ))

    #if ec:
        #display_ec(ec,csaOrig,csaRenu)

    #cmd="zcat %s/data/lig_all.tsv.gz|grep -P '^%s\\t%s\\t'"%(
        #rootdir,pdbid,asym_id)
    #p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    #stdout,stderr=p.communicate()
    #lines=stdout.decode().splitlines()
    #if len(lines):
        #print('''
#<tr><td>
#<div id="headerDiv">
    #<div id="titleText">Interaction with ligand</div>
#</div>
#<div style="clear:both;"></div>
#<div id="contentDiv">
    #<div id="RContent" style="display: block;">
    #<table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    #<tr BGCOLOR="#DEDEDE" align=center>
        #<th width=5%><strong>Site<br>#</strong></th>
        #<th width=5%><strong>Ligand</strong></th>
        #<th width=5%><strong>Ligand<br>chain</strong></th>
        #<th width=29%><strong>Binding residues on receptor<br>(original residue number in PDB)</strong></th>
        #<th width=29%><strong>Binding residues on receptor<br>(residue number reindexed from 1)</strong></th>
        #<th width=27%><strong>Binding affinity</strong></th>
    #</tr>
#''')
        #for l,line in enumerate(lines):
            #items    =line.split('\t')
            #recCha   =items[1]
            #bs       =items[2]
            #ccd      =items[3]
            #ligCha   =items[4]
            #ligIdx   =items[5]
            #resOrig  =items[6]
            #resRenu  =items[7]
            #manual   =items[8]
            #moad     =items[9]
            #pdbbind  =items[10]
            #bindingdb=items[11]

            #baff_list=[]
            #if manual:
                #baff_list.append("Manual survey: "+manual)
            #if moad:
                #baff_list.append("<a href=http://bindingmoad.org/pdbrecords/index/%s target=_blank>MOAD</a>: %s"%(pdbid,moad))
            #if pdbbind:
                #baff_list.append("<a href=http://pdbbind.org.cn/quickpdb.php?quickpdb=%s target=_blank>PDBbind-CN</a>: %s"%(pdbid,pdbbind))
            #if bindingdb:
                #baff_list.append("BindingDB: "+bindingdb)

            #bgcolor=''
            #if l%2==1:
                #bgcolor=' BGCOLOR="#DEDEDE" '
            #print('''
    #<tr $bgcolor align=center>
        #<td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$bs</a></span></td>
        #<td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$ccd</a></span></td>
        #<td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$ligCha</a></span></td>
        #<td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resOrig</a></span></td>
        #<td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resRenu</a></span></td>
        #<td>$baff</td>
    #</tr>
            #'''.replace("$bgcolor",bgcolor
              #).replace("$pdbid",pdbid
              #).replace("$recCha",recCha
              #).replace("$ccd",ccd
              #).replace("$ligCha",recCha
              #).replace("$bs",bs
              #).replace("$resOrig",resOrig
              #).replace("$resRenu",resRenu
              #).replace("$baff",'<br>'.join(baff_list)))

        #print('''   </table>
#</div>
#</td></tr>
#''')

    #if go:
        #display_go(go,uniprot,pdbid,asym_id)
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

    display_receptor(rna_info_list)
    #if lig3 and ligCha and ligIdx:
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
