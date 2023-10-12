#!/usr/bin/perl
my $docstring=<<EOF
isearch.pl jobID
    perform Infernal on output/jobID.fasta

isearch.pl
    perform Infernal on all jobs at output/
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
    foreach my $line(`ls $rootdir/output/|grep .isearch.html`)
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

if (!-s "$output/input.fasta")
{
    if (-s "$output.fasta")
    {
        system("mv $output.fasta $output/input.fasta");
    }
}
my $sequence="";
foreach my $line(`grep -v '>' $output/input.fasta`)
{
    chomp($line);
    $sequence.="$line";
}
my $L1=length $sequence;

if (!-s "$output/infernal.tblout")
{
    my $cmd="$bindir/cmscan -Z 1 --cpu 1 --tblout $output/infernal.tblout $rootdir/data/rna.cm $output/input.fasta > $output/infernal.out";
    system("$cmd");
    if (!-f "$output/infernal.tblout")
    {
        print "no infernal hit\n";
        exit(2);
    }
}
if (!-s "$output/hmmer.tblout")
{
    my $cmd="$bindir/hmmscan -Z 1 --cpu 1 --tblout $output/hmmer.tblout $rootdir/data/rna.hmm $output/input.fasta > $output/hmmer.out";
    system("$cmd");
    if (!-f "$output/hmmer.tblout")
    {
        print "no hmmer hit\n";
        exit(2);
    }
}

my %rna_dict;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz|cut -f1,2,16,17`)
{
    chomp($line);
    my @items  =split(/\t/,$line);
    my $pdbid  =$items[0];
    my $asym_id=$items[1];
    my $title  =$items[2];
    my $name   =$items[3];
    $rna_dict{"$pdbid:$asym_id"}="$title;\n$name";
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
Search <a href=input.fasta>input.fasta</a> (L=$L1) through the database by Infernal (cmscan). Results are ranked in ascending order of E-value.<br>

<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th ALIGN=center><strong> # </strong></th>
    <th ALIGN=center><strong> Hit</strong></th>
    <th ALIGN=center><strong> Hit<br>length</strong></th>
    <th ALIGN=center><strong> Strand</strong></th>
    <th ALIGN=center><strong> E-value</strong></th>
    <th ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
EOF
;

my %len_dict;
foreach my $line(`$bindir/fasta2len.py $rootdir/data/rna_nr.fasta`)
{
    if ($line=~/(\S+)\t(\d+)/)
    {
        my $target="$1";
        my $L2="$2";
        $len_dict{$target}=$L2;
    }
}

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



my $totalNum=0;
$txt="";
foreach my $line(`grep -v '^#' $output/infernal.tblout`)
{
    my @items =split(/\s+/,$line);
    my $sacc  =$items[0];
    my $strand=$items[9];
    my $evalue=$items[15];
    $txt.="$evalue\t$sacc\t$strand\n";
}
foreach my $line(`grep -v '^#' $output/hmmer.tblout`)
{
    my @items =split(/\s+/,$line);
    my $sacc  =$items[0];
    my $evalue=$items[4];
    $txt.="$evalue\t$sacc\t+\n";
}

open(FOUT,">$output/unsorted.tsv");
print FOUT $txt;
close(FOUT);
system("sort -gk1 $output/unsorted.tsv > $output/sorted.tsv");

foreach my $line(`cat $output/sorted.tsv`)
{
    chomp($line);
    my @items =split(/\s+/,$line);
    my $evalue=$items[0];
    my $sacc  =$items[1];
    my $strand=$items[2];

    $totalNum++;
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
    my $name="";
    if (exists($rna_dict{"$pdbid:$asym_id"}))
    {
        $name=$rna_dict{"$pdbid:$asym_id"};
    }
    my $hit="<span title=\"$name\"><a href=../../search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a></span>";
    my $homolog_line;
    if (exists $nr_dict{$sacc})
    {
        foreach my $mem (split(/,/,$nr_dict{$sacc}))
        {
            if ($mem ne $sacc && $target2id{$mem}=~/(\w+):(\w+)/)
            {
                my $pdbid="$1";
                my $asym_id="$2";
                my $name="";
                if (exists($rna_dict{"$pdbid:$asym_id"}))
                {
                    $name=$rna_dict{"$pdbid:$asym_id"};
                }
                $homolog_line.=", <span title=\"$name\"><a href=../../search.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$pdbid:$asym_id</a></span>";
            }
        }
    }
    $homolog_line=substr($homolog_line,2) if (length $homolog_line);
    my $L2=0;
    $L2=$len_dict{$sacc} if (exists $len_dict{$sacc});
    print FP <<EOF
<tr $bgcolor ALIGN=center>
    <td>$totalNum</td>
    <td>$hit</td>
    <td>$L2</td>
    <td>$strand</td>
    <td>$evalue</td>
    <td>$homolog_line</td>
</tr>
EOF
;
}

my $infernal_out=`cat $output/infernal.out`;
my $hmmer_out=`cat $output/hmmer.out`;

print FP <<EOF
</table><p></p>
Infernal result:
<pre>
$infernal_out
</pre>
HMMER result:
<pre>
$hmmer_out
</pre>
<a href=../..>[Back]</a>
$html_end
EOF
;
close(FP);

system("mv $output.html $output/old.html");
open(FP,">$output.html");
print FP <<EOF
<html>
<head>
<title>Infernal search</title>
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
