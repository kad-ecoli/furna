#!/usr/bin/perl
my $docstring=<<EOF
fasta2taxon.pl
    extract taxon from fasta
EOF
;
use strict;
use File::Basename;
use Cwd 'abs_path';

my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);
system("mkdir -p $rootdir/pdb/fasta/entry");
my %target2taxon;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,12|grep -P '\\t\$'|uniq`)
{
    if ($line=~/^(\S+)/)
    {
        my $pdbid="$1";
        my $filename="$rootdir/pdb/fasta/entry/$pdbid";
        if (!-s "$filename")
        {
            print "$filename\n";
            system("curl -s https://www.rcsb.org/fasta/entry/$pdbid > $filename");
        }
        foreach my $header(`grep '^>' $filename`)
        {
            chomp($header);
            #print "$header\n";
            my $taxon="";
            if ($header=~/ \((\d+)\)$/)
            {
                $taxon="$1";
            }
            else
            {
                next;
            }
            my @items=split(/\|/,$header);
            foreach my $chain(split(/,/,$items[1]))
            {
                my $chainid;
                if ($chain=~/\w+\[auth (\w+)\]/)
                {
                    $chainid="$1";
                }
                elsif ($chain=~/(\w+)$/)
                {
                    $chainid="$1";
                }
                my $target="$pdbid:$chainid";
                print "$target => $taxon\n";
                $target2taxon{$target}=$taxon;
            }
        }
    }
}

my $txt;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz`)
{
    my @items=split(/\t/,$line);
    if ($items[11] eq "")
    {
        my $pdbid=$items[0];
        my $chainid=$items[1];
        my $target="$pdbid:$chainid";
        if (exists($target2taxon{$target}))
        {
            $items[11]=$target2taxon{$target};
            $line=join("\t",@items);
        }
    }
    $txt.="$line";
}
open(FP,">$rootdir/data/rna.tsv");
print FP $txt;
close(FP);
#&gzipFile("$rootdir/data/rna.tsv");

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
