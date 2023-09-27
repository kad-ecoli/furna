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
outdir =os.path.join(rootdir,"output")

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

def display_interaction_list(pdbid,asym_id,ligand_dict):
    idx=0

    assemblyfile="output/"+pdbid+asym_id+".cif.gz"
    assemblyviewscript="load "+assemblyfile+"; color background black; select protein or rna or dna; spacefill off; wireframe off; cartoons; select not chain like "+asym_id+" and not hetero; color grey; select rna and chain like "+asym_id+" and not hetero; color group; select hetero; spacefill 70%; wireframe on;"

    fp=gzip.open("%s/data/interaction.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        if pdbid!=items[0] or asym_id!=items[1]:
            continue
        if idx==0:
            assembly=''
            if items[2]!="0":
                assembly="in assembly"+items[2]
            print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Interaction partners</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td>

    <table><tr align=left>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$assemblyviewscript;"
                }
                $("#assemblyview").html(Jmol.getAppletHtml("jmolApplet2",Info))
            });
            </script>
            <span id=assemblyview></span>
        </td>
        <td>
            All interaction partners of the RNA<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background black')">Black background</a>]<br>
            [<a href=$assemblyfile>Download all interaction partners</a>]<br>
        </td>
    </tr></table>


    </td></tr>
    <tr><td><table width=100%>
    <tr align=center bgcolor="#DEDEDE">
        <th>Ligand name<br>(Click to view)</th>
        <th>Ligand chain<br>$assembly</th>
        <th>Residue sequence<br>number of the ligand</th>
        <th>Ligand binding nucleotides on the RNA<br>
            (original residue number in PDB)</th>
        <th>Ligand binding nucleotides on the RNA<br>
            (residue number reindexed from 1)</th>
    </tr>
    '''.replace("$assemblyviewscript",assemblyviewscript
       ).replace("$assembly",assembly))
        idx+=1
        ccd=items[3]
        ligCha=items[4]
        ligIdx=items[5]
        residueOriginal=items[6]
        residueRenumber=items[7]
        resSeq=items[8]
        ccd_http=ccd
        if ccd in ligand_dict:
            ccd_http='<span title="%s">%s</span>'%(
                ligand_dict[ccd][4],ccd)
        bgcolor=''
        if idx % 2 == 0:
            bgcolor='bgcolor="#DEDEDE"'
        print('''
    <tr align=center $bgcolor>
        <td><a href=pdb.cgi?pdbid=$pdbid&chain=$chainid&lig3=$lig3&ligCha=$ligCha&ligIdx=$ligIdx>$ccd_http</a></td>
        <td>$ligCha</td>
        <td>$resSeq</td>
        <td>$residueOriginal</td>
        <td>$residueRenumber</td>
    </tr>
        '''.replace("$lig3",ccd
          ).replace("$pdbid",pdbid
          ).replace("$bgcolor",bgcolor
          ).replace("$chainid",asym_id
          ).replace("$ccd_http",ccd_http
          ).replace("$ligCha",ligCha
          ).replace("$ligIdx",ligIdx
          ).replace("$residueOriginal",residueOriginal
          ).replace("$residueRenumber",residueRenumber
          ).replace("$resSeq",resSeq))

    fp.close()
    if idx:
        print("</table></td></tr>")
        print("</table></div></td></tr>")
    return

def display_regular_ligand(ligand_info_list,ligand_dict):
    pdbid          =ligand_info_list[0]
    asym_id        =ligand_info_list[1]
    assembly       =ligand_info_list[2]
    ccd            =ligand_info_list[3]
    ligCha         =ligand_info_list[4]
    ligIdx         =ligand_info_list[5]
    residueOriginal=ligand_info_list[6]
    residueRenumber=ligand_info_list[7]
    resSeq         =ligand_info_list[8]
    prefix=pdbid+asym_id

    formula =''
    InChI   =''
    InChIKey=''
    SMILES  =''
    name    =''
    ChEMBL  =''
    DrugBank=''
    ZINC    =''
    if ccd in ligand_dict:
        formula,InChI,InChIKey,SMILES,name,ChEMBL,DrugBank,ZINC=ligand_dict[ccd]

    ligandID_txt="PDB: <a href=https://www.rcsb.org/ligand/%s target=_blank>%s</a>"%(ccd,ccd)
    if ChEMBL:
        ligandID_txt+="<br>ChEMBL: <a href=https://ebi.ac.uk/chembl/compound_report_card/%s target=_blank>%s</a>"%(ChEMBL,ChEMBL)
    if ZINC:
        ligandID_txt+="<br>ZINC: <a href=http://zinc.docking.org/substances/%s target=_blank>%s</a>"%(ZINC,ZINC)
    if DrugBank:
        ligandID_txt+="<br>DrugBank: <a href=http://go.drugbank.com/drugs/%s target=_blank>%s</a>"%(DrugBank,DrugBank)

    svgfile="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(ccd[0],ccd)
    assemblyfile="output/"+prefix+".cif.gz"
    complexfile="output/"+'_'.join((prefix,lig3,ligCha,ligIdx+".ent.gz"))
    ligandfile="output/"+'_'.join((prefix,lig3,ligCha,ligIdx+".pdb.gz"))

    lbs=[]
    if residueOriginal:
        lbs=[r[1:] for r in residueOriginal.split()]
    ligandselect="chain like "+ligCha+" and resno "
    if not '~' in resSeq:
        ligandselect+="="+resSeq
    else:
        resi1,resi2=resSeq.split('~')
        ligandselect+=">="+resi1+" and resno <="+resi2
        
        
    
    assemblyviewscript="load "+assemblyfile+"; color background black; select protein or rna or dna; spacefill off; wireframe off; cartoons; select not chain like "+asym_id+" and not ("+ligandselect+"); color grey; select rna and chain like "+asym_id+" and not hetero; color group; "
    localviewscript="load "+complexfile+"; color background black; spacefill off; wireframe off; "
    globalviewscript="load "+complexfile+"; color background black; spacefill off; wireframe off; "
    if len(lbs):
        assemblyviewscript+=' select chain like '+asym_id+' and '+ \
            '('+','.join(lbs)+'); spacefill 25%; wireframe 50;'
        assemblyviewscript+=' select chain like '+asym_id+" and *.C1' and "+ \
            '('+','.join(lbs)+'); label %m%R; color label magenta;'
        localviewscript+=' select chain=R and '+ \
            '('+','.join(lbs)+'); spacefill 25%; wireframe 50;'
        localviewscript+=" select chain=R and *.C1' and "+ \
            '('+','.join(lbs)+'); label %m%R; color label magenta;'
        globalviewscript+=' select chain=R and '+ \
            '('+','.join(lbs)+'); spacefill 25%; wireframe 50;'
        globalviewscript+=" select chain=R and *.C1' and "+ \
            '('+','.join(lbs)+'); label %m%R; color label magenta;'
    assemblyviewscript+="select "+ligandselect+"; spacefill 70%; frame all;"
    globalviewscript+="select chain=R; cartoons; color group; select chain=L; spacefill 70%; frame all;"
    localviewscript+="select chain=R; color group; select chain=L; spacefill 70%; frame all; zoom {chain=L}"

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Ligand information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>Ligand ID</strong></td><td>$ligandID_txt</td></tr>
    <tr bgcolor="#DEDEDE"><td align=center><strong>Ligand 2D<br>structure</strong></td><td><a href="$svgfile" target=_blank><img src="$svgfile" style="display:block;" width=400></a></td></tr>
    <tr><td align=center width=10%><strong>Formula</strong></td><td>$formula</td></tr>

    <tr bgcolor="#DEDEDE"><td align=center><strong>SMILES</strong></td><td>$smiles</td></tr>
    <tr><td align=center width=10%><strong>InCHI</strong></td><td>$inchi<br>(InChIKey=$inchikey)</td></tr>


    <tr bgcolor="#DEDEDE"><td align=center><strong>Ligand<br>chain</strong></td><td>Chain $ligCha $assembly</td></tr>
    <tr><td align=center width=10%><strong>Residue<br>sequence<br>number of<br>the ligand</strong></td><td>$resSeq</td></tr>
    </table>
    </div>
</div>
</tr></td>

<tr><td>
<div id="headerDiv">
    <div id="titleText">Interactions</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>Complex<br>structure</strong></td><td>

    <table><tr align=left>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$assemblyviewscript;"
                }
                $("#assemblyview").html(Jmol.getAppletHtml("jmolApplet2",Info))
            });
            </script>
            <span id=assemblyview></span>
        </td>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$globalviewscript"
                }
                $("#globalview").html(Jmol.getAppletHtml("jmolApplet3",Info))
            });
            </script>
            <span id=globalview></span>
        </td>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$localviewscript"
                }
                $("#localview").html(Jmol.getAppletHtml("jmolApplet4",Info))
            });
            </script>
            <span id=localview></span>
        </td>
    </tr>
    <tr align=left>
        <td>
            All interaction partners of the RNA<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background black')">Black background</a>]<br>
            [<a href=$assemblyfile>Download all interaction partners</a>]<br>
        </td>
        <td>
            Global structure of the ligand-RNA pair<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'color background black')">Black background</a>]<br>
            [<a href=$complexfile>Download ligand-RNA complex</a>]<br>
        </td>
        <td>
            Local structure of the ligand binding site<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'color background black')">Black background</a>]<br>
            [<a href=$ligandfile>Download ligand</a>]<br>
        </td>
    </tr></table>
    </td></tr>

    <tr bgcolor="#DEDEDE"><td align=center><strong>Ligand<br>binding<br>nucleotides<br>on the RNA</strong></td><td><li>Original residue number in PDB: $residueOriginal</li><li>Residue number reindexed from 1: $residueRenumber</li></td></tr>
    '''.replace("$ligandID_txt",ligandID_txt
      ).replace("$svgfile",svgfile
      ).replace("$formula",formula
      ).replace("$smiles",SMILES.replace('; ',';<br>')
      ).replace("$inchikey",InChIKey
      ).replace("$inchi",InChI
      ).replace("$ligCha",ligCha
      ).replace("$resSeq",resSeq
      ).replace("$residueOriginal",residueOriginal
      ).replace("$residueRenumber",residueRenumber
      ).replace("$assemblyfile",assemblyfile
      ).replace("$complexfile",complexfile
      ).replace("$ligandfile",ligandfile
      ).replace("$assemblyviewscript",assemblyviewscript
      ).replace("$globalviewscript",globalviewscript
      ).replace("$localviewscript",localviewscript
      ).replace("$assembly",'' if assembly=='0' else " in assembly"+assembly
    ))

    print('''
    </table>
    </div>
</div>
</tr></td>
    ''')

    return

def display_polymer_ligand(ligand_info_list,taxon_dict,parent_dict,ec_dict,
    go_dict,rfam_dict,rna_info_list,rnacentral_dict,uniprot_dict):
    pdbid          =ligand_info_list[0]
    asym_id        =ligand_info_list[1]
    assembly       =ligand_info_list[2]
    ccd            =ligand_info_list[3]
    ligCha         =ligand_info_list[4]
    ligIdx         =ligand_info_list[5]
    residueOriginal=ligand_info_list[6]
    residueRenumber=ligand_info_list[7]
    resSeq         =ligand_info_list[8]
    ligSequence    =ligand_info_list[9]
    receptor_taxon_list=rna_info_list[11].split(',')
    sequence       =rna_info_list[12]
    prefix=pdbid+asym_id

    ligseq_txt=''
    width=100
    L2=len(ligSequence)
    for i in range(0,L2,width):
        ligseq_txt+=ligSequence[i:(i+width)]+' - %d<br>'%(i+width if i+width<L2 else L2)
    ligseq_txt="<tr><td align=center width=10%><strong>$ccd<br>sequence</strong></td><td><font color=magenta>$ligseq_txt</font></td></tr>".replace("$ligseq_txt",ligseq_txt
        ).replace("$ccd",ccd.replace("rna","RNA").replace("dna","DNA").replace("protein","Protein"))


    assemblyfile="output/"+prefix+".cif.gz"
    complexfile="output/"+'_'.join((prefix,lig3,ligCha,ligIdx+".ent.gz"))
    ligandfile="output/"+'_'.join((prefix,lig3,ligCha,ligIdx+".pdb.gz"))

    lbs=[]
    if residueOriginal:
        lbs=[r[1:] for r in residueOriginal.split()]
    ligandselect="chain like '"+ligCha+"' and resno "
    if not '~' in resSeq:
        ligandselect+="="+resSeq
    else:
        resi1,resi2=resSeq.split('~')
        ligandselect+=">="+resi1+" and resno <="+resi2
    
    assemblyviewscript="load "+assemblyfile+"; color background black; select protein or rna or dna; spacefill off; wireframe off; cartoons; select not chain like "+asym_id+" and not ("+ligandselect+"); color grey; select rna and chain like "+asym_id+" and not hetero; color group; "
    localviewscript="load "+complexfile+"; color background black; spacefill off; wireframe off; "
    globalviewscript="load "+complexfile+"; color background black; spacefill off; wireframe off; "
    if len(lbs):
        localviewscript+=' select chain=R and '+ \
            '('+','.join(lbs)+'); spacefill 25%; wireframe 50;'
        localviewscript+=" select chain=R and *.C1' and "+ \
            '('+','.join(lbs)+'); label %m%R; color label grey;'
        globalviewscript+=' select chain=R and '+ \
            '('+','.join(lbs)+'); spacefill 25%; wireframe 50;'
        globalviewscript+=" select chain=R and *.C1' and "+ \
            '('+','.join(lbs)+'); label %m%R; color label grey;'
    assemblyviewscript+="select "+ligandselect+"; cartoons; color magenta;"
    globalviewscript+="select chain=R; cartoons; color group; select chain=L; cartoons; color magenta;"
    localviewscript+="select chain=R; color group; select chain=L; cartoons; color magenta; zoom {chain=L}"


    polymer_table=''
    name=''
    if ccd=="protein":
        fp=gzip.open("data/protein.tsv.gz",'rt')
        protein_lines=fp.read().splitlines()
        fp.close()
        for line in protein_lines:
            items=line.split('\t')
            if items[0]!=pdbid or items[1]!=ligCha.split('-')[0]:
                continue
            uniprot_list=items[3].split(',')
            ec_list     =items[4].split(',')
            mf_list     =items[5].split(',')
            bp_list     =items[6].split(',')
            cc_list     =items[7].split(',')
            taxon_list  =items[8].split(',')
            name        =items[10]
            
            if len(uniprot_list) and uniprot_list[0]:
                polymer_table+='<tr bgcolor="#DEDEDE"><td align=center><strong>UniProt</strong></td><td>'
                for u,uniprot in enumerate(uniprot_list):
                    if u:
                        polymer_table+="<br>"
                    uniprot_html="<a href=https://uniprot.org/uniprotkb/%s target=_blank>%s</a>"%(uniprot,uniprot)
                    if uniprot in uniprot_dict:
                        uniprot_html=uniprot_dict[uniprot].replace('|'+uniprot+'|','|'+uniprot_html+'|')
                    polymer_table+=uniprot_html
                polymer_table+="</td></tr>"


            taxon_txt=''
            for t,taxon in enumerate(taxon_list):
                if t:
                    taxon_txt+="<br>"
                if not taxon in taxon_dict:
                    taxon_dict[taxon]=''
                taxon_txt+="<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a> <i>%s</i>"%(
                    taxon,taxon,taxon_dict[taxon])
            if len(taxon_txt):
                polymer_table+='<tr><td align=center><strong>Species</strong></td><td>'+taxon_txt+"</td></tr>"

            mf_parent =''
            bp_parent =''
            cc_parent =''
            target=pdbid+':'+ligCha.split('-')[0]
            if target in parent_dict:
                mf_parent=parent_dict[target][0]
                bp_parent=parent_dict[target][1]
                cc_parent=parent_dict[target][2]

            go_table=''
            if mf_parent+bp_parent+cc_parent:
                go_table+='''<tr bgcolor="#DEDEDE"><td align=center><strong>GO<br>terms</strong></td><td>
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
                bgcolor='bgcolor="#DEDEDE"'
                if mf_parent+bp_parent+cc_parent:
                    bgcolor=''
                ec_table+='<tr '+bgcolor+'><td align=center><strong>EC<br>numbers</strong></td><td>'
                for e,ec in enumerate(ec_list):
                    if e:
                        ec_table+="<br>"
                    ec=ec.replace('EC:','')
                    ec_table+="<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(ec,ec)
                    if ec in ec_dict:
                        ec_table+=" "+ec_dict[ec]
                ec_table+="</td></tr>"

            polymer_table+=go_table+ec_table

    elif ccd=="rna":
        fp=gzip.open("data/rna.tsv.gz",'rt')
        rna_lines=fp.read().splitlines()
        fp.close()
        for line in rna_lines:
            items=line.split('\t')
            if items[0]!=pdbid or items[1]!=ligCha.split('-')[0]:
                continue
            rfam_list   =items[4].split(',')
            rc_list     =items[5].split(',')
            ec_list     =items[7].split(',')
            mf_list     =items[8].split(',')
            bp_list     =items[9].split(',')
            cc_list     =items[10].split(',')
            taxon_list  =items[11].split(',')
            cssr        =items[13]
            dssr        =items[14].replace('&','')
            if len(items)>15:
                name    =items[15]
           
            ligseq_txt=''
            width=100
            L2=len(ligSequence)
            for i in range(0,L2,width):
                ligseq_txt+=ligSequence[i:(i+width)]+' - %d<br>'%(i+width if i+width<L2 else L2)
                ligseq_txt+='<span title="CSSR secondary structure assignment">'+cssr[i:(i+width)]+'</span><br>'
                ligseq_txt+='<span title="DSSR secondary structure assignment">'+dssr[i:(i+width)]+'</span><br>'

            ligseq_txt="<tr><td align=center width=10%><strong>RNA<br>sequence &amp;<br>secondary<br>structure</strong></td><td><font color=magenta>$ligseq_txt</font></td></tr>".replace("$ligseq_txt",ligseq_txt)

            taxon_txt=''
            for t,taxon in enumerate(taxon_list):
                if t:
                    taxon_txt+="<br>"
                if not taxon in taxon_dict:
                    taxon_dict[taxon]=''
                taxon_txt+="<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a> <i>%s</i>"%(
                    taxon,taxon,taxon_dict[taxon])
            if len(taxon_txt):
                polymer_table+='<tr bgcolor="#DEDEDE"><td align=center><strong>Species</strong></td><td>'+taxon_txt+"</td></tr>"

            if len(rc_list) and rc_list[0]:
                polymer_table+='<tr><td align=center><strong>RNAcentral</strong></td><td>'
                for r,rnacentral in enumerate(rc_list):
                    if r:
                        polymer_table+="<br>"
                    polymer_table+="<a href=https://rnacentral.org/rna/%s target=_blank>%s</a>"%(rnacentral,rnacentral)
                    if rnacentral in rnacentral_dict:
                        polymer_table+=rnacentral_dict[rnacentral]
                polymer_table+="</td></tr>"

            if len(rfam_list) and rfam_list[0]:
                polymer_table+='<tr bgcolor="#DEDEDE"><td align=center><strong>Rfam<br>families</strong></td><td>'
                for r,rfam in enumerate(rfam_list):
                    if r:
                        polymer_table+="<br>"
                    polymer_table+="<a href=https://rfam.org/family/%s target=_blank>%s</a> "%(rfam,rfam)
                    if rfam in rfam_dict:
                        polymer_table+=rfam_dict[rfam]
                polymer_table+="</td></tr>"

            mf_parent =''
            bp_parent =''
            cc_parent =''
            target=pdbid+':'+ligCha.split('-')[0]
            if target in parent_dict:
                mf_parent=parent_dict[target][0]
                bp_parent=parent_dict[target][1]
                cc_parent=parent_dict[target][2]

            go_table=''
            if mf_parent+bp_parent+cc_parent:
                go_table+='''<tr><td align=center><strong>GO<br>terms</strong></td><td>
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

            polymer_table+=go_table+ec_table
    elif ccd=="dna":
        fp=gzip.open("data/dna.tsv.gz",'rt')
        dna_lines=fp.read().splitlines()
        fp.close()
        for line in dna_lines:
            items=line.split('\t')
            if items[0]!=pdbid or items[1]!=ligCha.split('-')[0]:
                continue
            taxon_list  =items[3].split(',')
            ligSequence =items[4]
            name        =items[5]

            if (len(taxon_list)==0 or len(taxon_list[0])==0) and (
                len(receptor_taxon_list) and len(receptor_taxon_list[0])) and (
                sequence.lower().replace('u','t')==ligSequence.upper()[::-1
                    ].replace('A','t').replace('U','a').replace('T','a'
                    ).replace('C','g').replace('G','c')):
                taxon_list=receptor_taxon_list
           
            ligseq_txt=''
            width=100
            L2=len(ligSequence)
            for i in range(0,L2,width):
                ligseq_txt+=ligSequence[i:(i+width)]+' - %d<br>'%(i+width if i+width<L2 else L2)

            ligseq_txt="<tr><td align=center width=10%><strong>DNA<br>sequence</td><td><font color=magenta>$ligseq_txt</font></td></tr>".replace("$ligseq_txt",ligseq_txt)

            taxon_txt=''
            for t,taxon in enumerate(taxon_list):
                if t:
                    taxon_txt+="<br>"
                if not taxon in taxon_dict:
                    taxon_dict[taxon]=''
                taxon_txt+="<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a> <i>%s</i>"%(
                    taxon,taxon,taxon_dict[taxon])
            if len(taxon_txt):
                polymer_table+='<tr bgcolor="#DEDEDE"><td align=center><strong>Species</strong></td><td>'+taxon_txt+"</td></tr>"

    name=name.strip()
    if name.startswith('(') and name.endswith(')'):
        name=name[1:-1]
    if len(name):
        name='('+name+')'

    cssr_table=''
    if ccd in ["dna","rna"]:
        cmd="script/CSSR output/%s%s_%s_%s_0.ent.gz - -o 1"%(pdbid,asym_id,ccd,ligCha)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        stdout=stdout.decode()
        width=100
        L1=len(sequence)
        L2=len(ligSequence)
        cssr1=stdout[:L1]
        cssr2=stdout[L1:]
        colorbar=[]
        if len(stdout):
            if ccd in "dna":
                cssr_table="<tr><td align=center><strong>RNA/DNA<br>sequence &amp;<br>secondary<br>structure</strong></td><td>"
            else:
                cssr_table="<tr><td align=center><strong>RNA<br>sequence &amp;<br>secondary<br>structure</strong></td><td>"
            seq_list=list(sequence)
            cssr_list=list(cssr1)
            #colorbar=["red","cyan","green","yellow","orange","red"]
            #for i in range(L1):
                #color=colorbar[int(i*len(colorbar)/L1)]
                #seq_list[i]="<font color="+color+">"+seq_list[i]+"</font>"
                #cssr_list[i]="<font color="+color+">"+cssr_list[i]+"</font>"

            for i in range(0,L1,width):
                seq=sequence[i:(i+width)]
                ssr=cssr1[i:(i+width)]
                cssr_table+=seq+' - %d<br>'%(i+width if i+width<L1 else L1
                    )+'<span title="CSSR secondary structure assignment">'+ssr+'</span><br>'
            cssr_table+="<font color=magenta>"
            for i in range(0,L2,width):
                cssr_table+=ligSequence[i:(i+width)]+' - %d<br>'%(i+width if i+width<L2 else L2)+'<span title="CSSR secondary structure assignment">'+cssr2[i:(i+width)]+'</span><br>'
            cssr_table+="</font></td></tr>"

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Interaction partner information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    $ligseq_txt
    <tr bgcolor="#DEDEDE"><td align=center><strong>Ligand<br>chain</strong></td><td>Chain $ligCha $assembly<br>$name</td></tr>
    <tr><td align=center width=10%><strong>Residue<br>sequence<br>number of<br>the ligand</strong></td><td>$resSeq</td></tr>
    $polymer_table
    </table>
    </div>
</div>
</tr></td>

<tr><td>
<div id="headerDiv">
    <div id="titleText">Interactions</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>Complex<br>structure</strong></td><td>

    <table><tr align=left>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$assemblyviewscript;"
                }
                $("#assemblyview").html(Jmol.getAppletHtml("jmolApplet2",Info))
            });
            </script>
            <span id=assemblyview></span>
        </td>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$globalviewscript"
                }
                $("#globalview").html(Jmol.getAppletHtml("jmolApplet3",Info))
            });
            </script>
            <span id=globalview></span>
        </td>
        <td>
            <script type="text/javascript"> 
            $(document).ready(function()
            {
                Info = {
                    width: 400,
                    height: 400,
                    j2sPath: "jsmol/j2s",
                    script: "$localviewscript"
                }
                $("#localview").html(Jmol.getAppletHtml("jmolApplet4",Info))
            });
            </script>
            <span id=localview></span>
        </td>
    </tr>
    <tr align=left>
        <td>
            All interaction partners of the RNA<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet2, 'color background black')">Black background</a>]<br>
            [<a href=$assemblyfile>Download all interaction partners</a>]<br>
        </td>
        <td>
            Global structure of the ligand-RNA pair<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet3, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet3, 'color background black')">Black background</a>]<br>
            [<a href=$complexfile>Download ligand-RNA complex</a>]<br>
        </td>
        <td>
            Local structure of the ligand binding site<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'spin on')">Spin on</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'spin off')">Spin off</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'Reset')">Reset orientation</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'set antialiasDisplay true')">High quality</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'set antialiasDisplay false')">Low quality</a>]<br>
            [<a href="javascript:Jmol.script(jmolApplet4, 'color background white')">White background</a>]
            [<a href="javascript:Jmol.script(jmolApplet4, 'color background black')">Black background</a>]<br>
            [<a href=$ligandfile>Download ligand</a>]<br>
        </td>
    </tr></table>
    </td></tr>

    <tr bgcolor="#DEDEDE"><td align=center><strong>Ligand<br>binding<br>nucleotides<br>on the RNA</strong></td><td><li>Original residue number in PDB: $residueOriginal</li><li>Residue number reindexed from 1: $residueRenumber</li></td></tr>
    $cssr_table
    '''.replace("$ligseq_txt",ligseq_txt
      ).replace("$ligCha",ligCha
      ).replace("$resSeq",resSeq
      ).replace("$polymer_table",polymer_table
      ).replace("$residueOriginal",residueOriginal
      ).replace("$residueRenumber",residueRenumber
      ).replace("$assemblyfile",assemblyfile
      ).replace("$complexfile",complexfile
      ).replace("$ligandfile",ligandfile
      ).replace("$assemblyviewscript",assemblyviewscript
      ).replace("$globalviewscript",globalviewscript
      ).replace("$localviewscript",localviewscript
      ).replace("$assembly",'' if assembly=='0' else " in assembly"+assembly
      ).replace("$name",name
      ).replace("$cssr_table",cssr_table
    ))


    print('''
    </table>
    </div>
</div>
</tr></td>
    ''')

    return

def display_receptor(rna_info_list,taxon_dict,parent_dict,ec_dict,
    go_dict,rfam_dict,rnacentral_dict,fimo_dict,pubmed_dict):
    pdbid     =rna_info_list[0]
    asym_id   =rna_info_list[1]
    L     =int(rna_info_list[2])
    reso      =rna_info_list[3]
    rfam_list =rna_info_list[4].split(',')
    rc_list   =rna_info_list[5].split(',')
    ec_list   =rna_info_list[7].strip().split(',')
    mf_list   =rna_info_list[8].split(',')
    bp_list   =rna_info_list[9].split(',')
    cc_list   =rna_info_list[10].split(',')
    taxon_list=rna_info_list[11].split(',')
    sequence  =rna_info_list[12]
    cssr      =rna_info_list[13]
    dssr      =rna_info_list[14].replace('&','')
    title     =rna_info_list[15]
    if len(rna_info_list)>16:
        name  =rna_info_list[16]

    mf_parent =''
    bp_parent =''
    cc_parent =''
    target=pdbid+':'+asym_id
    if target in parent_dict:
        mf_parent=parent_dict[target][0]
        bp_parent=parent_dict[target][1]
        cc_parent=parent_dict[target][2]
    
    taxon_txt=''
    for t,taxon in enumerate(taxon_list):
        if t:
            taxon_txt+="<br>"
        if not taxon in taxon_dict:
            taxon_dict[taxon]=''
        taxon_txt+="<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a> <i>%s</i>"%(
            taxon,taxon,taxon_dict[taxon])

    rnacentral_table=''
    bgcolor='bgcolor="#DEDEDE"'
    if len(rc_list) and rc_list[0]:
        if bgcolor:
            bgcolor=''
        else:
            bgcolor='bgcolor="#DEDEDE"'
        rnacentral_list=[]
        for r in rc_list:
            rnacentral_list.append("<a href=https://rnacentral.org/rna/%s target=_blank>%s</a>"%(r,r))
            if r in rnacentral_dict:
                rnacentral_list[-1]+=rnacentral_dict[r]
        rnacentral_table+='<tr '+bgcolor+'''><td align=center><strong>RNAcentral</strong></td><td>%s</td></tr>'''%(''.join(rnacentral_list))

    rfam_table=''
    if len(rfam_list) and rfam_list[0]:
        if bgcolor:
            bgcolor=''
        else:
            bgcolor='bgcolor="#DEDEDE"'
        rfam_table+='<tr '+bgcolor+'><td align=center><strong>Rfam<br>families</strong></td><td>'
        for r,rfam in enumerate(rfam_list):
            if r:
                rfam_table+="<br>"
            rfam_table+="<a href=https://rfam.org/family/%s target=_blank>%s</a> "%(rfam,rfam)
            if rfam in rfam_dict:
                rfam_table+=rfam_dict[rfam]
        rfam_table+="</td></tr>"

    go_table=''
    if mf_parent+bp_parent+cc_parent:
        if bgcolor:
            bgcolor=''
        else:
            bgcolor='bgcolor="#DEDEDE"'
        go_table+="<tr "+bgcolor+'''><td align=center><strong>GO<br>terms</strong></td><td>
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
        if bgcolor:
            bgcolor=''
        else:
            bgcolor='bgcolor="#DEDEDE"'
        ec_table+='<tr '+bgcolor+'><td align=center><strong>EC<br>numbers</strong></td><td>'
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


    if len(name):
        if name.startswith('(') and name.endswith(')'):
            name=name[1:-1]
        name='('+name+')'

    fimo_table=''
    if target in fimo_dict:
        if bgcolor:
            bgcolor=''
        else:
            bgcolor='bgcolor="#DEDEDE"'
        
        cmd="zcat %s|grep -P '^ATOM  '|cut -c23-27|uniq"%chainfile
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        resi_list=stdout.decode().splitlines()
        #print(resi_list)

        fimo_table='<tr><td '+bgcolor+''' align=center><strong>Motif for<br>protein<br>binding</strong>
    </td><td><table width=100%>
    <tr align=center bgcolor="#DEDEDE">
        <th>Motif<br>ID</th>
        <th>Logo</th>
        <th>Matched sequence<br>(residue number) qvalue</th>
        <th>Binding<br>protein</th>
        <th>PubMed<br>for motif</th>
    </tr>'''
        motif_dict=dict()
        motif_list=[]
        for line in fimo_dict[target]:
            db       =line[0]
            resi1    =line[1]
            resi2    =line[2]
            qvalue   =line[3]
            Motif_seq=line[4].replace('T','U').lower()
            Motif_ID =line[5]
            Gene_name=line[6]
            Gene_id  =line[7]
            #citation =line[8]
            citation_list=[]
            if not "not available" in line[8]:
                for citation in line[8].split(';'):
                    if citation=="NULL":
                        continue
                    citation_list.append('''<a href=https://pubmed.ncbi.nlm.nih.gov/$citation target=_blank>$citation</a>'''.replace("$citation",citation))
                    if citation in pubmed_dict:
                        citation_list[-1]='<span title="'+pubmed_dict[citation] \
                            +'">'+citation_list[-1]+'</span>'
            pngfile="data/"+db+"/logo"+Motif_ID.replace('.','_').replace("_pfm",".pfm")+".png"
            if not db+'\t'+Motif_ID in motif_dict:
                motif_list.append((db,Motif_ID,Gene_name,Gene_id,
                    pngfile,'<br>'.join(citation_list)))
                motif_dict[db+'\t'+Motif_ID]=[]

            if len(resi_list) and int(resi1)<=len(resi_list
                            ) and int(resi2)<=len(resi_list):
                resi1=resi_list[int(resi1)-1].strip()
                resi2=resi_list[int(resi2)-1].strip()
            motif_dict[db+'\t'+Motif_ID].append("%s (%s ~ %s) qvalue=%s"%(
                Motif_seq,resi1,resi2,qvalue))

        for l,items in enumerate(motif_list):
            bgcolor=''
            if l % 2==1:
                bgcolor='bgcolor="#DEDEDE"'
            db,Motif_ID,Gene_name,Gene_id,pngfile,citation=items
            Motif_match=''.join(['<li>'+m+'</li>' for m in motif_dict[db+'\t'+Motif_ID]])

            if db=="attract":
                Motif_ID="ATtRACT:"+Motif_ID
            elif db=="cisbp":
                Motif_ID="CISBP-RNA:"+Motif_ID
                
            fimo_table+='''
    <tr align=center $bgcolor>
        <td>$Motif_ID</td>
        <td><a href=$pngfile><img src=$pngfile height=50></a></td>
        <td>$Motif_match</td>
        <td>$Gene_name<br>($Gene_id)</td>
        <td>$citation</td>
    </tr>'''.replace("$Motif_ID",Motif_ID
           ).replace("$resi1",resi1
           ).replace("$resi2",resi2
           ).replace("$bgcolor",bgcolor
           ).replace("$Motif_match",Motif_match
           ).replace("$Gene_name",Gene_name
           ).replace("$Gene_id",Gene_id
           ).replace("$citation",citation
           ).replace("$pngfile",pngfile
           )
            
        fimo_table+="</table></td></tr>"

    

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">RNA information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>Sequence &amp;<br>secondary<br>structure</strong></td><td>$seq_txt</td></tr>
    <tr bgcolor="#DEDEDE"><td align=center><strong>PDB</strong></td><td><span title="Search other entries from the same structure"><a href=search.cgi?pdbid=$pdbid target=_blank>$pdbid</a></span> $title</td></tr>
    <tr><td align=center><strong>Chain</strong></td><td><span title="Search other entries for the same chain"><a href=search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$asym_id</a></span> $name</td></tr>
    <tr bgcolor="#DEDEDE"><td align=center><strong>Resolution</strong></td><td>$reso</a></td></tr>
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
        <!--
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
    -->
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
        <!--
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
        -->
    </tr></table>


    </td></tr>
    <tr bgcolor="#DEDEDE"><td align=center><strong>Species</strong></td><td>$taxon_txt</td></tr>
    $rnacentral_table
    $rfam_table
    $go_table
    $ec_table
    $fimo_table
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
  ).replace("$go_table" ,go_table
  ).replace("$ec_table" ,ec_table
  ).replace("$fimo_table",fimo_table
  ).replace("$rnacentral_table",rnacentral_table
  ).replace("$name"     ,name
  ))
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
            break
    fp.close()
    return rna_info_list

def get_ligand_info(pdbid,asym_id,lig3,ligCha,ligIdx):
    ligand_info_list=[]
    fp=gzip.open("%s/data/interaction.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        if not line.startswith(pdbid):
            continue
        items=line.split('\t')
        if items[0]!=pdbid or items[1]!=asym_id or items[3]!=lig3 or items[4]!=ligCha:
            continue
        if items[5]==ligIdx or lig3 in ["protein","rna","dna"]:
            ligand_info_list=items
            break
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

def read_ligand():
    ligand_dict=dict()
    fp=gzip.open("%s/data/ligand.tsv.gz"%rootdir,'rt')
    for line in fp.read().splitlines():
        if not line.startswith('#'):
            items=line.split('\t')
            ligand_dict[items[0]]=items[1:]
    fp.close()
    return ligand_dict

def read_rnacentral():
    rnacentral_dict=dict()
    fp=gzip.open("data/RNAcentral.fasta.gz",'rt')
    for line in fp.read().splitlines():
        if not line.startswith('>'):
            continue
        line=line[1:]
        r=line.split()[0]
        name=line[len(r):]
        rnacentral_dict[r]=name
    fp.close()
    return rnacentral_dict

def read_uniprot():
    uniprot_dict=dict()
    fp=gzip.open(rootdir+"/data/uniprot.tsv.gz",'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        uniprot_dict[items[0]]=items[1]
    fp.close()
    return uniprot_dict

def read_fimo():
    fimo_dict=dict()
    if os.path.isfile(rootdir+"/data/attract_fimo.tsv.gz"):
        fp=gzip.open(rootdir+"/data/attract_fimo.tsv.gz",'rt')
        for line in fp.read().splitlines():
            items=line.split('\t')
            key=':'.join(items[:2])
            if not key in fimo_dict:
                fimo_dict[key]=[]
            fimo_dict[key].append(["attract"]+items[2:])
        fp.close()
    if os.path.isfile(rootdir+"/data/cisbp_fimo.tsv.gz"):
        fp=gzip.open(rootdir+"/data/cisbp_fimo.tsv.gz",'rt')
        for line in fp.read().splitlines():
            items=line.split('\t')
            key=':'.join(items[:2])
            if not key in fimo_dict:
                fimo_dict[key]=[]
            fimo_dict[key].append(["cisbp"]+items[2:])
        fp.close()
    return fimo_dict

def extract_ligand(pdbid,asym_id,lig3,ligCha,ligIdx):
    divided=pdbid[-3:-1]
    prefix=pdbid+asym_id
    interimfile=os.path.join(rootdir,"interim",divided,prefix+".tar.gz")
    outdir=os.path.join(rootdir,"output")
    
    ligandbase='_'.join((prefix,lig3,ligCha,ligIdx+".pdb"))
    ligandfile=os.path.join(outdir,ligandbase+".gz")
    if not os.path.isfile(ligandfile):
        txt=''
        tar = tarfile.open(interimfile)
        fin=tar.extractfile(ligandbase)
        for line in fin.read().decode().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" L"+line[22:]+'\n'
            else:
                txt+=line+'\n'
        fout=gzip.open(ligandfile,'wt')
        fout.write(txt)
        fout.close()
        fin.close()
        tar.close()

    complexbase='_'.join((prefix,lig3,ligCha,ligIdx+".ent"))
    complexfile=os.path.join(outdir,complexbase+".gz")
    if not os.path.isfile(complexfile):
        txt=''
        chainfile=os.path.join(rootdir,"chain",divided,prefix+".pdb.gz")
        fin=gzip.open(chainfile,'rt')
        for line in fin.read().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" R"+line[22:]+'\n'
        fin.close()
        txt+='TER\n'
        fin=gzip.open(ligandfile,'rt')
        for line in fin.read().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" L"+line[22:]+'\n'
        fin.close()
        txt+='END\n'
        fout=gzip.open(complexfile,'wt')
        fout.write(txt)
        fout.close()
    return

def extract_assembly(pdbid,asym_id):
    divided=pdbid[-3:-1]
    prefix=pdbid+asym_id
    interimfile=os.path.join(rootdir,"interim",divided,prefix+".tar.gz")
    if not os.path.isfile(interimfile):
        return
    outdir=os.path.join(rootdir,"output")
    
    assemblyfile=os.path.join(outdir,prefix+".cif.gz")
    if not os.path.isfile(assemblyfile):
        txt='''data_%s
#
_entry.id %s
#
loop_
_atom_site.group_PDB 
_atom_site.auth_atom_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
'''%(prefix,prefix)
        chainfile=os.path.join(rootdir,"chain",divided,prefix+".pdb.gz")
        fin=gzip.open(chainfile,'rt')
        for line in fin.read().splitlines():
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            txt+=' '.join((
                line[:6].strip(),          # group_PDB
                line[12:16].strip(),       # auth_atom_id
                line[17:20].strip(),       # auth_comp_id
                asym_id,                   # auth_asym_id
                line[22:26].strip(),       # auth_seq_id
                line[26].replace(' ','?'), # pdbx_PDB_ins_code
                line[30:38].strip(),       # Cartn_x
                line[38:46].strip(),       # Cartn_y
                line[46:54].strip(),       # Cartn_z
                ))+'\n'
        fin.close()

        indexfile=os.path.join(rootdir,"interim",divided,prefix+".txt")
        member_list=[]
        fin = open(indexfile)
        for line in fin.read().split('#CCD\t')[-1].splitlines()[1:]:
            items=line.split('\t')
            member_list.append(prefix+'_'+items[0]+'_'+items[1]+'_'+items[2]+".pdb")
        fin.close()

        tar = tarfile.open(interimfile)
        for member in member_list:
            chainID=member.split('_')[2]
            fin=tar.extractfile(member)
            for line in fin.read().decode().splitlines():
                if not line.startswith("ATOM") and not line.startswith("HETATM"):
                    continue
                txt+=' '.join((
                    line[:6].strip(),          # group_PDB
                    line[12:16].strip(),       # auth_atom_id
                    line[17:20].strip(),       # auth_comp_id
                    chainID,                   # auth_asym_id
                    line[22:26].strip(),       # auth_seq_id
                    line[26].replace(' ','?'), # pdbx_PDB_ins_code
                    line[30:38].strip(),       # Cartn_x
                    line[38:46].strip(),       # Cartn_y
                    line[46:54].strip(),       # Cartn_z
                ))+'\n'
            fin.close()

        txt+='#\n'
        fout=gzip.open(assemblyfile,'wt')
        fout.write(txt)
        fout.close()
    return

def display_asymetric_unit(pdbid):
    print("<tr><td><h1 align=center> PDB "+pdbid+"</h1></td></tr>")
    html_end='''
</table>
<p></p>[<a href=.>Back to home</a>]
</body> </html>'''


    rna_info_mat=[]
    fp=gzip.open("data/rna.tsv.gz",'rt')
    for line in fp.read().splitlines():
        if not line.startswith(pdbid):
            continue
        items=line.split('\t')
        if pdbid!=items[0]:
            continue
        rna_info_mat.append(items)
    fp.close()

    if len(rna_info_mat)==0:
        print(html_end)
        exit()
    
    asym_id_list=[]
    reso=''
    pubmed_list=[]
    title=''
    for rna_info_list in rna_info_mat:
        asym_id_list.append(rna_info_list[1])
        if len(reso)==0:
            reso=rna_info_list[3]
            pubmed_list=rna_info_list[6].split(',')
            title=rna_info_list[15]
    if reso=="-1.00" or reso=="NA":
        reso="N/A"
    else:
        reso=reso+" &#8491;"

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Structure information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>PDB</strong></td><td><span title="Search other entries from the same structure"><a href=search.cgi?pdbid=$pdbid target=_blank>$pdbid</a></span> $title</td></tr>
    <tr bgcolor="#DEDEDE"><td align=center><strong>Resolution</strong></td><td>$reso</a></td></tr>
    </table>
</div>
</td></tr>

<tr><td>
<div id="headerDiv">
    <div id="titleText">Chain information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr bgcolor="#DEDEDE" align=left>
        <th>Chain</th>
        <th>Rfam</th>
        <th>RNAcentral</th>
        <th>EC number &amp; GO term</th>
        <th>Species</th>
        <th>Ligand</th>
        <th>Sequence &amp; length</th>
    </tr>
'''.replace("$pdbid"    ,pdbid
  ).replace("$title"    ,title
  ).replace("$reso"     ,reso
  ))


    taxon_dict=read_taxon()
    parent_dict=read_parent()
    go_dict=read_go()
    ec_dict=read_ec()
    rfam_dict=read_rfam()
    ligand_dict=read_ligand()

    ligand_info_dict=dict()
    assembly=''
    fp=gzip.open("data/interaction.tsv.gz",'rt')
    for line in fp.read().splitlines():
        if not line.startswith(pdbid):
            continue
        items=line.split('\t')
        if pdbid!=items[0]:
            continue
        asym_id=items[1]
        ccd   =items[3]
        ligCha=items[4]
        ligIdx=items[5]
        residueOriginal=items[6]
        resSeq=items[8].replace(' ','')
        ligSequence=items[9]
        if not asym_id in ligand_info_dict:
            ligand_info_dict[asym_id]="Assembly"+items[2]


        if ccd in ["protein","rna","dna"]:
            if ccd=="protein":
                spantitle="Protein"
            else:
                spantitle=ccd.upper()
            spantitle+=" sequence: "
            if len(ligSequence)<100:
                spantitle+=ligSequence
            else:
                spantitle+=ligSequence[:100]+"..."
        elif ccd in ligand_dict:
            spantitle=ligand_dict[ccd][4]
        else:
            spantitle=ccd
        spantitle+="\nChain ID of the ligand: %s\nResidue sequence number of the ligand: %s\nBinding nucleotides on the RNA receptor: %s\n"%(ligCha,resSeq,residueOriginal)


        bs_list=[bs[1:] for bs in residueOriginal.split(' ')]
        bs_txt=''
        write_prev=0
        maxWidth=30
        curWidth=len(ccd)+len(ligCha)+len(resSeq)+3
        for b,bs in enumerate(bs_list):
            if b==0:
                bs_txt+=bs
            else:
                if not bs[-1] in "1234567890" or not bs_list[b-1][-1] in "1234567890":
                    if curWidth>maxWidth:
                        bs_txt+=",<br>"+bs
                        curWidth=len(bs)
                    else:
                        bs_txt+=","+bs
                        curWidth+=1+len(bs)
                else:
                    if int(bs_list[b-1])+1==int(bs):
                        if b+1==len(bs_list) or not bs_list[b+1][-1] in "1234567890" or int(bs)+1!=int(bs_list[b+1]):
                            if curWidth>maxWidth:
                                bs_txt+="~<br>"+bs
                                curWidth=len(bs)
                            else:
                                bs_txt+="~"+bs
                                curWidth+=1+len(bs)
                    else:
                        if curWidth>maxWidth:
                            bs_txt+=",<br>"+bs
                            curWidth=len(bs)
                        else:
                            bs_txt+=","+bs
                            curWidth+=1+len(bs)

        ligand_info_dict[asym_id]+='''<br><span title="$spantitle"><a href=pdb.cgi?pdbid=$pdbid&chain=$asym_id&lig3=$ccd&ligCha=$ligCha&ligIdx=$ligIdx>$ccd:$ligCha:$resSeq</a> $bs_txt</span>
        '''.replace("$spantitle",spantitle
          ).replace("$pdbid",pdbid
          ).replace("$asym_id",asym_id
          ).replace("$ccd",ccd
          ).replace("$ligCha",ligCha
          ).replace("$resSeq",resSeq
          ).replace("$bs_txt",bs_txt
          ).replace("$ligIdx",ligIdx)
    fp.close()


    for r,rna_info_list in enumerate(rna_info_mat):
        bgcolor='bgcolor="#DEDEDE"'
        if r%2==0:
            bgcolor=''
        asym_id   =rna_info_list[1]
        L     =int(rna_info_list[2])
        rfam_list =rna_info_list[4].split(',')
        rc_list   =rna_info_list[5].split(',')
        ec_list   =rna_info_list[7].strip().split(',')
        mf_list   =rna_info_list[8].split(',')
        bp_list   =rna_info_list[9].split(',')
        cc_list   =rna_info_list[10].split(',')
        taxon_list=rna_info_list[11].split(',')
        sequence  =rna_info_list[12]
        cssr      =rna_info_list[13]
        dssr      =rna_info_list[14].replace('&','')
        title     =rna_info_list[15]
        if len(rna_info_list)>16:
            name  =rna_info_list[16]

        for r,rfam in enumerate(rfam_list):
            rfam_list[r]="<a href=https://rfam.org/family/%s target=_blank>%s</a>"%(rfam,rfam)
            if rfam in rfam_dict:
                rfam_list[r]+=' '+rfam_dict[rfam]
        rfam_table='<br>'.join(rfam_list)

        for r,rnacentral in enumerate(rc_list):
            rc_list[r]="<a href=https://rnacentral.org/rna/%s target=_blank>%s</a> "%(rnacentral,rnacentral.split('_')[0])
        rnacentral_table='<br>'.join(rc_list)

        go_list=mf_list+bp_list+cc_list
        for g,go in enumerate(go_list):
            go_list[g]="<a href=https://www.ebi.ac.uk/QuickGO/term/%s target=_blank>%s</a>"%(go,go)
            if go in go_dict:
                go_list[g]+=" (%s) %s"%(go_dict[go][0],go_dict[go][1])
        for e,ec in enumerate(ec_list):
            ec=ec.replace('EC:','')
            ec_list[e]="<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(ec,ec)
            if ec in ec_dict:
                ec_list[e]+=" "+ec_dict[ec]
        ecgo='<br>'.join(go_list+ec_list)

        for t,taxon in enumerate(taxon_list):
            if not taxon in taxon_dict:
                taxon_dict[taxon]=''
            taxon_list[t]="<i>%s</i> (<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s target=_blank>%s</a>)"%(
                    taxon_dict[taxon],taxon,taxon)
        taxon_table='<br>'.join(taxon_list)

        seq_list=[]
        width=50
        for i in range(0,L,width):
            seq_list.append(sequence[i:(i+width)])
        seq_txt="&gt;Length=$len<br>".replace("$pdbid",pdbid).replace("$asym_id",asym_id).replace("$len",str(L))
        seq_txt+='<br>'.join(seq_list)

        if not asym_id in ligand_info_dict:
            ligand_info_dict[asym_id]=''

        print('''
    <tr $bgcolor>
        <td><a href=pdb.cgi?pdbid=$pdbid&chain=$asym_id>$asym_id</a> $name</td>
        <td>$rfam</td>
        <td>$RNAcentral</td>
        <td>$ecgo</td>
        <td>$taxon</td>
        <td>$ligand_info</td>
        <td>$seq_txt</td>
    </tr>'''.replace("$bgcolor",bgcolor
           ).replace("$pdbid",pdbid
           ).replace("$asym_id",asym_id
           ).replace("$name",name
           ).replace("$rfam",rfam_table
           ).replace("$ecgo",ecgo
           ).replace("$RNAcentral",rnacentral_table
           ).replace("$taxon",taxon_table
           ).replace("$ligand_info",ligand_info_dict[asym_id]
           ).replace("$L",str(L)
           ).replace("$seq_txt",seq_txt
           ))


    pubmed_html=''
    if len(pubmed_list):
        pubmed_dict=dict()
        fp=gzip.open("%s/data/pubmed.tsv.gz"%rootdir,'rt')
        for line in fp.read().splitlines():
            items=line.split('\t')
            pubmed_dict[items[0]]=items[1]
        fp.close()
        pubmed_html='''<tr bgcolor="#DEDEDE"><td align=center><strong>PubMed</strong></td><td>'''
        for p in pubmed_list:
            pubmed="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(p,p)
            if p in pubmed_dict:
                pubmed=pubmed_dict[p]+' '+pubmed
            pubmed_html+="<li>"+pubmed+"</li>"
        pubmed_html+="</td></tr>"

    print('''
</table></td></tr>
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
    $pubmed
    </tr>
'''.replace("$pdbid",pdbid
  ).replace("$pubmed",pubmed_html
  ))
    
    print(html_end)
    return

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
   
    if pdbid and not asym_id:
        display_asymetric_unit(pdbid)
        exit()

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
    ligand_dict=read_ligand()
    rnacentral_dict=read_rnacentral()
    uniprot_dict=read_uniprot()
    fimo_dict=read_fimo()

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
        pubmed='''<tr bgcolor="#DEDEDE">'''
        pubmed+='''<td align=center><strong>PubMed</strong></td><td>%s</td></tr>'''%(''.join(pubmed_list))

    display_receptor(rna_info_list,taxon_dict,parent_dict,ec_dict,
        go_dict,rfam_dict,rnacentral_dict,fimo_dict,pubmed_dict)
    extract_assembly(pdbid,asym_id)
    if lig3 and ligCha and ligIdx:
        ligand_info_list=get_ligand_info(pdbid,asym_id,lig3,ligCha,ligIdx)
        if ligand_info_list:
            extract_ligand(pdbid,asym_id,lig3,ligCha,ligIdx)
            if not lig3 in ["protein","dna","rna"]:
                display_regular_ligand(ligand_info_list,ligand_dict)
            else:
                display_polymer_ligand(ligand_info_list,taxon_dict,
                    parent_dict,ec_dict,go_dict,rfam_dict,rna_info_list,
                    rnacentral_dict,uniprot_dict)
    else:
        display_interaction_list(pdbid,asym_id,ligand_dict)
   
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
    $pubmed
    </tr>
'''.replace("$pdbid",pdbid
  ).replace("$pubmed",pubmed
  ))
    
    print('''</table>
<p></p>[<a href=.>Back to home</a>]
</body> </html>''')
