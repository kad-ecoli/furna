#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "combine dna level annotation\n";
my @target_list;
my %target_dict;
print "$rootdir/data/interaction.tsv.gz\n";
foreach my $line(`zcat $rootdir/data/interaction.tsv.gz |cut -f1,4,5,10|grep -P "\\tdna\\t"|uniq|sort|uniq`)
{
    if ($line=~/^(\w+)\tdna\t([-\w]+)\t(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence="$3";
        if ($chainid=~/(\w+)\-\w+/)
        {
            $chainid="$1";
        }
        my $target="$pdbid:$chainid";
        if (!exists($target_dict{$target}))
        {
            push(@target_list,($target));
            $target_dict{$target}=$sequence;
        }
    }
}
my $target_num=scalar @target_list;
print "$target_num dna\n";

my %taxon_dict;
my $taxon_num=0;
foreach my $line(`zcat $rootdir/data/dna.tsv.gz|grep -v '^#'`)
{
    if ($line=~/^(\w+)\t(\w+)\t\d+\t(\d+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $taxon="$3";
        my $target="$pdbid:$chainid";
        $taxon_dict{$taxon}="$taxon";
        $taxon_num++;
    }
}
print "$taxon_num dna with taxon\n";
foreach my $target(@target_list)
{
    next if (exists($taxon_dict{$target}));
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        my $cmd="$bindir/cif2chain $rootdir/pdb/data/structures/divided/mmCIF/$divided/$pdbid.cif.gz - $chainid|grep ORGANISM_TAXID";
        print "$cmd\n";
        if (`$cmd`=~/ORGANISM_TAXID:\s*(\d+)/)
        {
            my $taxon="$1";
            $taxon_dict{$target}="$taxon";
        }
    }
}

my %name_dict;
foreach my $line(`grep '>' $rootdir/pdb/derived_data/pdb_seqres.txt |grep -F mol:na`)
{
    chomp($line);
    if ($line=~/>(\w+)_(\w+)\s+mol:na\s+length:\d+\s+([\s\S]+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $name="$3";
        my $target="$pdbid:$chainid";
        $name_dict{$target}=$name;
    }
}

my $txt="#pdb\tchain\tL\ttaxon\tsequence\tname\n";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence=$target_dict{$target};
        my $L=length $sequence;
        my $taxon="";
        my $name="";
        if (exists($taxon_dict{$target}))
        {
            $taxon=$taxon_dict{$target};
        }
        if (exists($name_dict{$target}))
        {
            $name=$name_dict{$target};
        }
        $txt.="$pdbid\t$chainid\t$L\t$taxon\t$sequence\t$name\n";
    }
}
open(FP,">$rootdir/data/dna.tsv");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/dna.tsv");

$txt="";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence=$target_dict{$target};
        $txt.=">$pdbid$chainid\n$sequence\n";
    }
}
open(FP,">$rootdir/data/dna.fasta");
print FP "$txt";
close(FP);
system("$bindir/fasta_nr.py $rootdir/data/dna.fasta $rootdir/data/dna_nr.fasta $rootdir/data/dna_nr.tsv");

exit();

sub gzipFile
{
    my ($filename)=@_;
    my $oldNum=`zcat $filename.gz 2>/dev/null|wc -l`+0;
    my $newNum=` cat $filename   |wc -l`+0;
    if (0.8*$oldNum>$newNum)
    {
        print "WARNING! do not update $filename from $oldNum to $newNum entries\n";
        return;
    }
    print "update $filename from $oldNum to $newNum entries\n";
    system("gzip -f $filename");
    return;
}
