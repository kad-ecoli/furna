#!/usr/bin/perl
my $docstring=<<EOF
fsearch.pl jobID
    perform USalign on output/jobID.cif or output/jobID.pdb

fsearch.pl
    perform USalign on all jobs at output/
EOF
;

use strict;
use File::Basename;
use Cwd 'abs_path';
use Fcntl qw(:flock);
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

if (scalar @ARGV==0)
{
    flock DATA, LOCK_EX|LOCK_NB or die "Unable to lock file $!";
    foreach my $line(`ls $rootdir/output/|grep .fsearch.html`)
    {
        chomp($line);
        if ($line=~/(\S+).html$/)
        {
            my $input="$1";
            next if (-s "$rootdir/output/$input/index.html");
            my $cmd=abs_path(__FILE__)." $input";
            #print "timeout 1h $cmd\n";
            system("timeout 1h $cmd");
        }
    }
    exit(0);
}

my $input =$ARGV[0];
my $output="$rootdir/output/$input";

if (!-d "$output/")
{
    system("mkdir -p $output/");
}

my $infile="input.pdb";
my $sequence;
if (!-s "$output/input.pdb" && !-s "$output/input.cif")
{
    if (-s "$output.pdb")
    {
        system("mv $output.pdb $output/input.pdb");
        $sequence=`$bindir/pdb2fasta $output/input.pdb |grep -v '^>'`;
    }
    else
    {
        if (!-s "$output.cif")
        {
            print "no such file $output.cif\n";
            exit(1);
        }
        system("mv $output.cif $output/input.cif");
        $sequence=`$bindir/pdb2fasta $output/input.cif |grep -v '^>'`;
    }
}
$infile="input.cif" if (-s "$output/input.cif");
chomp($sequence);
my $L1=length $sequence;

if (!-s "$output/aln.m8")
{
    my $cmd="cd $output; cat $infile | $bindir/USalign - -dir2 $rootdir/data/USalign/ $rootdir/data/USalign/cluster.list -outfmt 2 -fast > $output/aln.m8";
    system("$cmd");
    if (!-f "$output/aln.m8")
    {
        print "no usalign hit\n";
        exit(2);
    }
}

#### prepare output ####
my $txt=`cat $rootdir/index.html`;
my @items=split(/<!-- CONTENT START -->/,$txt);
my $html_header=$items[0];
@items   =split(/<!-- CONTENT END -->/,$txt);
my $html_end=$items[1];
@items   =split(/<!-- ==BBB== ======ending of head======== -->/,$html_header);
$html_header=~s/href="/href="..\/..\//g;
$html_header=~s/href="..\/..\/http/href="http/g;
$html_header=~s/href="..\/..\/\//href="\//g;
open(FP,">$output/index.html");
print FP <<EOF
$html_header
Search <a href=$infile>$infile</a> (L=$L1) through the database by US-align (-fast mode). Results are ranked in descending order of TM-score normalized by query.<br>

<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th ALIGN=center><strong> # </strong></th>
    <th ALIGN=center><strong> Hit</strong></th>
    <th ALIGN=center><strong> TM-score<br> normalized <br>by query</strong> </th>           
    <th ALIGN=center><strong> TM-score<br> normalized <br>by hit</strong> </th>           
    <th ALIGN=center><strong> RMSD </strong> </th>           
    <th ALIGN=center><strong> Identity<br> normalized <br>by query</strong> </th>           
    <th ALIGN=center><strong> Identity<br> normalized <br>by hit</strong> </th>           
    <th ALIGN=center><strong> Identity<br> normalized <br> by aligned <br>length</strong> </th>           
    <th ALIGN=center><strong> Hit<br>length</strong> </th>           
    <th ALIGN=center><strong> Aligned<br>length</strong></th>
    <th ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
EOF
;


my %target2id;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2`)
{
    if ($line=~/(\w+)\t(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        $target2id{"$pdbid$chainid"}="$pdbid:$chainid";
    }
}


my %nr_dict;
foreach my $line(`cat $rootdir/data/rna_nr.tsv`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $key="$items[0]";
    if (scalar @items>1)
    {
        $nr_dict{$key}=join(",",@items);
    }
}



my %hit2clust_dict;
foreach my $line(`cat $rootdir/data/USalign/cluster.txt`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $key="$items[0]";
    if (scalar @items>1)
    {
        my $last=scalar @items;
        my $txt=$nr_dict{$items[0]};
        for (my $i=1;$i<scalar @items;$i++)
        {
            $txt.=",$nr_dict{$items[$i]}";
        }
        $hit2clust_dict{$key}=$txt;
    }
}

my $totalNum=0;
foreach my $line(`cat $output/aln.m8|sort -k3nr|grep -v '^#'`)
{
    my @items=split(/\t/,$line);
    my $sacc =$items[1];
    my $TM1  =$items[2];
    my $TM2  =$items[3];
    my $RMSD =$items[4];
    my $ID1  =$items[5];
    my $ID2  =$items[6];
    my $IDali=$items[7];
    my $L1   =$items[8];
    my $L2   =$items[9];
    my $Lali =$items[10];
    $totalNum++;
    last if ($totalNum>100 && $TM1<0.5);
    my $bgcolor='';
       $bgcolor='BGCOLOR="#DEDEDE"' if ($totalNum % 2==0);

    if ($sacc=~/(\w+):\w*/)
    {
        $sacc="$1";
    }
    my $pdbid="$sacc";
    my $asym_id="";
    if (exists($target2id{$sacc}) && $target2id{$sacc}=~/(\w+):(\w+)/)
    {
        $pdbid="$1";
        $asym_id="$2";
    }
    my $hit="<a href=../../search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a>";
    my $homolog_line;
    if (exists $hit2clust_dict{$sacc})
    {
        foreach my $mem (split(/,/,$hit2clust_dict{$sacc}))
        {
            if ($mem ne $sacc && $target2id{$mem}=~/(\w+):(\w+)/)
            {
                my $pdbid="$1";
                my $asym_id="$2";
                $homolog_line.=", <a href=../../search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a>";
            }
        }
    }
    $homolog_line=substr($homolog_line,2) if (length $homolog_line);
    print FP <<EOF
<tr $bgcolor ALIGN=center>
    <td>$totalNum</td>
    <td>$hit</td>
    <td>$TM1</td>
    <td>$TM2</td>
    <td>$RMSD</td>
    <td>$ID1</td>
    <td>$ID2</td>
    <td>$IDali</td>
    <td>$L2</td>
    <td>$Lali</td>
    <td>$homolog_line</td>
</tr>
EOF
;
}

print FP <<EOF
</table><p></p><a href=../..>[Back]</a>
$html_end
EOF
;
close(FP);

system("mv $output.html $output/old.html");
open(FP,">$output.html");
print FP <<EOF
<html>
<head>
<title>Structure search</title>
</head>
<body>
<meta http-equiv="refresh" content="0; url='$input'"/>
</body>
</html>
EOF
;
close(FP);

exit(0);
__DATA__
