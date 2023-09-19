#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "curate cisbp\n";
system("mkdir -p $rootdir/fimo/");
my @taxon_list;
my %taxon2fasta;
foreach my $line(`cat $rootdir/data/TF_Information_all_motifs.txt|cut -f8|grep -v TF_Species|uniq|sort|uniq`)
{
    chomp($line);
    push(@taxon_list,($line));
    $taxon2fasta{$line}="";
}

foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2,12,13|grep -v '^#'`)
{
    if ($line=~/^(\w+)\t(\w+)\t(\d+)\t(\w+)/)    
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $taxon="$3";
        my $sequence="$4";
        $sequence=~s/u/t/g;
        $sequence=~s/i/a/g;
        $sequence=uc($sequence);
        next if (!exists($taxon2fasta{$taxon}));
        $taxon2fasta{$taxon}.=">${pdbid}_$chainid\n$sequence\n";
    }
}
foreach my $taxon(@taxon_list)
{
    next if (-s "$rootdir/fimo/$taxon.tsv");

    open(FP,">$rootdir/fimo/$taxon.fasta");
    print FP "$taxon2fasta{$taxon}";
    close(FP);

    my $txt=<<EOF
MEME version 4

ALPHABET= ACGT
strands: +

Background letter frequencies
EOF
;
    $txt.=`cat $rootdir/genomes/$taxon.freq`;

    my @Motif_list=`cat $rootdir/data/TF_Information_all_motifs.txt|cut -f4,8|grep -P "\\t$taxon\$"|cut -f1|sort|uniq`;
    foreach my $Motif_ID(@Motif_list)
    {
        chomp($Motif_ID);
        print "$Motif_ID\t$taxon\n";
        my $filename="$rootdir/cisbp/meme/$Motif_ID.txt";
        #$txt.=`cat $rootdir/cisbp/meme/$Motif_ID.txt` if (-s "$filename");
        if (-s "$filename")
        {
            $txt.="\n";
            $txt.=`grep -PA1000 "^MOTIF " $filename`;
            #$txt.="\n";
        }
    }
    open(FP,">$rootdir/fimo/$taxon.meme.txt");
    print FP $txt;
    close(FP);

    my $cmd="$bindir/fimo -o $rootdir/fimo/$taxon $rootdir/fimo/$taxon.meme.txt $rootdir/fimo/$taxon.fasta";
    print  "$cmd\n";
    system("$cmd");
    $cmd="cat $rootdir/fimo/$taxon/fimo.tsv |grep -vP '\\t\\d+\\t\\d+\\t-\\t' > $rootdir/fimo/$taxon.tsv";
    print  "$cmd\n";
    system("$cmd");
    if (-s "$rootdir/fimo/$taxon.tsv")
    {
        system("rm -rf $rootdir/fimo/$taxon");
    }
}

my $txt="";
foreach my $taxon(@taxon_list)
{
    foreach my $line(`cat $rootdir/fimo/$taxon.tsv|cut -f1,3,4,5,10|sort -k2|grep -vP '^#'`)
    {
        if ($line=~/^(\S+)\t(\w+)_(\w+)\t(\d+)\t(\d+)\t(\w+)/)
        {
            my $Motif_ID="$1";
            my $pdbid   ="$2";
            my $chainid ="$3";
            my $start   ="$4";
            my $end     ="$5";
            my $sequence="$6";
            $txt.="$pdbid\t$chainid\t$Motif_ID\t$start\t$end\t$sequence\n";
        }
    }
}
open(FP,">$rootdir/data/fimo.tsv");
print FP $txt;
close(FP);

exit();
