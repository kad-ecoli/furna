#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("mkdir -p $rootdir/data/USalign");
my $txt="";
foreach my $line(`cut -f1 $rootdir/data/rna_nr.tsv`)
{
    chomp($line);
    print "$line\n";
    my $pdbid=substr($line,0,4);
    my $chainid=substr($line,4);
    if (length($line)>8 && substr($line,0,8)==lc(substr($line,0,8)))
    {
        $pdbid=substr($line,0,8);
        $chainid=substr($line,8);
    }
    my $divided=substr($pdbid,length($pdbid)-3,2);
    if (!-s "$rootdir/data/USalign/$line")
    {
        system("zcat $rootdir/arena/$divided/$line.pdb.gz|cut -c1-54|grep ATOM|grep \" C3'\" > $rootdir/data/USalign/$line");
    }
    if (!-s "$rootdir/data/USalign/$line")
    {
        system("zcat $rootdir/chain/$divided/$line.pdb.gz|cut -c1-54|grep ATOM|grep \" C3'\" > $rootdir/data/USalign/$line");
    }
    $txt.="$line\n";
}
open(FP,">$rootdir/data/USalign/list");
print FP $txt;
close(FP);

my %fasta_dict;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2,13`)
{
    if ($line=~/^(\w+)\t(\w+)\t(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence="$3";
        my $target="$pdbid$chainid";
        $fasta_dict{$target}=$sequence;
    }
}

$txt="";
foreach my $line(`cat $rootdir/data/USalign/list`)
{
    chomp($line);
    $txt.=">$line\n$fasta_dict{$line}\n";
}
open(FP,">$rootdir/data/USalign/list.fasta");
print FP $txt;
close(FP);

system("$bindir/cd-hit-est -i $rootdir/data/USalign/list.fasta -o $rootdir/data/USalign/c80.fasta -c 0.8 -s 0.5");
system("$bindir/clstr2tsv.py $rootdir/data/USalign/c80.fasta.clstr $rootdir/data/USalign/c80.tsv");
system("$bindir/qTMclust -fast -dir $rootdir/data/USalign/  $rootdir/data/USalign/list -o  $rootdir/data/USalign/cluster.txt -init $rootdir/data/USalign/c80.tsv");
system("cut -f1 $rootdir/data/USalign/cluster.txt > $rootdir/data/USalign/cluster.list");

exit();

sub gzipFile
{
    my ($filename)=@_;
    my $oldNum=`zcat $filename.gz 2>/dev/null|wc -l`+0;
    my $newNum=` cat $filename               |wc -l`+0;
    if (0.8*$oldNum>$newNum)
    {
        print "WARNING! do not update $filename from $oldNum to $newNum entries\n";
        return;
    }
    print "update $filename from $oldNum to $newNum entries\n";
    system("gzip -f $filename");
    return;
}
