#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

my $date=`date +%Y-%m-%d`;
chomp($date);
my $RNAs=`zcat $rootdir/data/rna.tsv.gz |cut -f1,2|grep -v '^#'|wc -l`+0;
my $ecn=`zcat $rootdir/data/rna.tsv.gz |cut -f8 |grep -F EC:|wc -l`+0;
my $mf =`zcat $rootdir/data/rna.tsv.gz |cut -f9 |grep -F GO:|wc -l`+0;
my $bp =`zcat $rootdir/data/rna.tsv.gz |cut -f10|grep -F GO:|wc -l`+0;
my $cc =`zcat $rootdir/data/rna.tsv.gz |cut -f11|grep -F GO:|wc -l`+0;
my $got=`zcat $rootdir/data/rna.tsv.gz |cut -f9-11|grep -F GO:|wc -l`+0;
my $attract=`zcat $rootdir/data/attract_fimo.tsv.gz |cut -f1,2|sort|uniq|wc -l`+0;
my $rfam=`zcat $rootdir/data/rna.tsv.gz |cut -f5|grep -F RF|wc -l`+0;

my %metal_dict;
foreach my $line(`zcat $rootdir/data/metal.tsv.gz|cut -f1`)
{
    chomp($line);
    $metal_dict{$line}=1;
}
my $protein=0;
my $dna=0;
my $rna=0;
my $metal=0;
my $regular=0;
my $interaction=0;
foreach my $line(`zcat $rootdir/data/interaction.tsv.gz |grep -v '^#'|cut -f4`)
{
    chomp($line);
    if ($line eq "protein")
    {
        $protein++;
    }
    elsif ($line eq "dna")
    {
        $dna++;
    }
    elsif ($line eq "rna")
    {
        $rna++;
    }
    elsif (exists($metal_dict{$line}))
    {
        $metal++;
    }
    else
    {
        $regular++;
    }
    $interaction++;
}

my $ligand=`zcat $rootdir/data/ligand.tsv.gz |cut -f1|grep -v '^#'|wc -l`+0;

my $index_txt=<<EOF
<p>
<h1><u>Database statistics</u></h1>
</p>

The database is updated weekly and the current version ($date) contains:
<li>Number of ligand-RNA interactions: <a href=search.cgi>$interaction</a></li>
<li>Number of regular ligand-RNA interactions: <a href=search.cgi?lig3=regular>$regular</a></li>
<li>Number of metal-ion-RNA interactions: <a href=search.cgi?lig3=metal>$metal</a></li>
<li>Number of protein-RNA interactions: <a href=search.cgi?lig3=protein>$protein</a></li>
<li>Number of DNA-RNA interactions: <a href=search.cgi?lig3=dna>$dna</a></li>
<li>Number of RNA-RNA interactions: <a href=search.cgi?lig3=rna>$rna</a></li>
<li>Number of types of ligands: <a href=ligand.cgi>$ligand</a></li>

<li>Number of RNAs: <a href=search.cgi>$RNAs</a></li>
<li>Number of RNAs with Enzyme Commission (EC) numbers: <a href=search.cgi?ecn=0>$ecn</a></li>
<li>Number of RNAs with Gene Ontology (GO) terms: <a href=search.cgi?&got=0>$got</a>
(<a href=search.cgi?&got=0003674>$mf</a> with Molecular Function, 
 <a href=search.cgi?&got=0008150>$bp</a> with Biological Process, and 
 <a href=search.cgi?&got=0005575>$cc</a> with Cellular Component)</li>
<li>Number of RNAs with ATtRACT motifs for protein binding: <a href=search.cgi?mot=attract>$attract</a></li>
<li>Number of RNAs with Rfam families: <a href=search.cgi?rfm=RF>$rfam</a></li>
EOF
;

open(FP,">$rootdir/data/index.txt");
print FP $index_txt;
close(FP);

#my $cmd="cd $rootdir; tar -cJf $rootdir/download/chain.tar.xz chain/";
my $cmd="cd $rootdir; tar -cf $rootdir/download/chain.tar chain/";
print  "$cmd\n";
system("$cmd");
#$cmd="cd $rootdir; tar -cJf $rootdir/download/ligand.tar.xz interim/*/*.tar.gz";
$cmd="cd $rootdir; tar -cf $rootdir/download/ligand.tar interim/*/*.tar.gz";
print  "$cmd\n";
system("$cmd");

exit();
