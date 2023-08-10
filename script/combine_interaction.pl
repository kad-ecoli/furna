#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "combine interaction annotation\n";

my %pubmed_dict;
print "$rootdir/sifts/pdb_pubmed.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_pubmed.tsv.gz`)
{
    if ($line=~/^(\w+)\t\d+\t(\d+)/)
    {
        my $pdbid=lc("$1");
        my $pubmed="$2";
        if (exists ($pubmed_dict{$pdbid}))
        {
            $pubmed_dict{$pdbid}.=",$pubmed";
        }
        else
        {
            $pubmed_dict{$pdbid}="$pubmed";
        }
    }
}

my $txt="#pdb\tchain\tassembly\tCCD\tligCha\tligIdx\t";
$txt.="residueOriginal\tresidueRenumber\tresSeq\tligandSequence\n";
print "$rootdir/data/rna.tsv\n";
foreach my $line(`grep -v '^#' $rootdir/data/rna.tsv|cut -f1,2,3`)
{
    if ($line=~/(\w+)\t(\w+)\t(\d+)/)
    {
        my $pdbid="$1";
        my $chainID="$2";
        my $L="$2";
        my $target="$pdbid$chainID";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        my %fasta_dict;
        next if (!-s "$rootdir/interim/$divided/$target.txt" ||
                 !-s "$rootdir/interim/$divided/$target.tar.gz");

        my $startLig=0;
        my $header="";
        print "$rootdir/interim/$divided/$target.txt\n";
        
        my $assembly="1";
        my $indir ="$rootdir/pdb/data/assemblies/mmCIF/divided/$divided";
        if (-s "$indir/${pdbid}-assembly2.cif.gz")
        {
            $assembly="0";
            foreach my $a(`ls $indir|grep ${pdbid}-assembly|grep -ohP "\\d+.cif.gz\$"|cut -f1 -d.|sort -n`)
            {
                chomp($a);
                if (`$bindir/cif2fasta $indir/${pdbid}-assembly$a.cif.gz|grep '^>'`=~/>XXXX:$chainID\t/)
                {
                    $assembly="$a";
                    print "$indir/${pdbid}-assembly$assembly.cif.gz => $target\n";
                    last;
                }
            }
        }
        foreach my $line(`tail -n+3 $rootdir/interim/$divided/$target.txt`)
        {
            chomp($line);
            if (length $line==0)
            {
                next;
            }
            elsif ($line=~/^#CCD/)
            {
                $startLig=1;
                next;
            }
            if ($startLig==0)
            {
                if ($line=~/^>([-\w]+)/)
                {
                    $header="$1";
                }
                else
                {
                    $fasta_dict{$header}="$line";
                }
            }
            else
            {
                my @items=split(/\t/,$line);
                my $ligCha=$items[1];
                my $sequence="";
                if (exists($fasta_dict{$ligCha}))
                {
                    $sequence=$fasta_dict{$ligCha};
                }
                $txt.="$pdbid\t$chainID\t$assembly";
                for (my $i=0;$i<scalar @items;$i++)
                {
                    $txt.="\t$items[$i]";
                }
                $txt.="\t$sequence\n";
            }
        }

    }
}
print "$rootdir/data/interaction.tsv\n";
open(FP,">$rootdir/data/interaction.tsv");
print FP "$txt";
close(FP);

&gzipFile("$rootdir/data/interaction.tsv");

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
