#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("$bindir/cif2fasta $rootdir/pdb/data/structures/divided/mmCIF/ $rootdir/pdb/derived_data/nuc.list > $rootdir/pdb/derived_data/nuc.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d{2,}' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/rna.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d\$' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/fragment.fasta");

system("mkdir -p $rootdir/data") if (!-d "$rootdir/data");
system("$bindir/SortFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $rootdir/pdb/derived_data/rna.fasta $rootdir/data/rna.fasta");
system("$bindir/SortFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $rootdir/pdb/derived_data/fragment.fasta $rootdir/data/fragment.fasta");
system("$bindir/fasta_nr.py $rootdir/data/rna.fasta $rootdir/data/rna_nr.fasta $rootdir/data/rna_nr.tsv");

foreach my $line(`cat $rootdir/pdb/derived_data/rna.fasta |grep '^>'`)
{
    if ($line=~/>(\w+):(\w+)\tRNA\t(\d+)/)
    {
        my $pdbid="$1";
        my $chainID="$2";
        my $L="$3";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        #print "$pdbid $chainID $L\n";
        
        system("mkdir -p $rootdir/chain/$divided") if (!-d "$rootdir/chain/$divided");
        my $filename="$rootdir/chain/$divided/$pdbid$chainID.pdb";
        if (!-s "$filename.gz")
        {
            print "$filename.gz\n";
            system("$bindir/cif2chain $rootdir/pdb/data/structures/divided/mmCIF/$divided/$pdbid.cif.gz $filename $chainID");
            system("gzip -f $filename");
            if (-s "$filename.gz")
            {
                system("zcat $filename.gz | cut -c1-54 > $filename");
                system("gzip -f $filename");
            }
        }

        system("mkdir -p $rootdir/cssr/$divided") if (!-d "$rootdir/cssr/$divided");
        my $cssrfile="$rootdir/cssr/$divided/$pdbid$chainID.cssr";
        if (!-s "$cssrfile")
        {
            print "$cssrfile\n";
            system("$bindir/CSSR $filename.gz $cssrfile -o 3");
        }

        system("mkdir -p $rootdir/dssr/$divided") if (!-d "$rootdir/dssr/$divided");
        system("mkdir -p $rootdir/arena/$divided") if (!-d "$rootdir/arena/$divided");
        my $arenafile="$rootdir/arena/$divided/$pdbid$chainID.pdb";
        my $dssrfile="$rootdir/dssr/$divided/$pdbid$chainID.dssr";

        if (!-s "$dssrfile" || !-s "$arenafile.gz")
        {
            print "$dssrfile\n";
            my $tmpdir="$rootdir/tmp/$pdbid$chainID";
            system("mkdir -p $tmpdir");
            system("$bindir/Arena $filename.gz $tmpdir/arena.pdb 5");
            system("cd $tmpdir; $bindir/x3dna-dssr -i=arena.pdb >/dev/null");
            my $sequence=`cat $tmpdir/dssr-2ndstrs.dbn |head -2|tail -1|sed 's/&//g'`;
            chomp($sequence);
            if (length $sequence ne $L)
            {
                system("$bindir/Arena $tmpdir/arena.pdb $tmpdir/arena.pdb 6");
                system("cd $tmpdir; $bindir/x3dna-dssr -i=arena.pdb >/dev/null");
            }
            if (`cat $tmpdir/dssr-2ndstrs.dbn|wc -l`+0==3)
            {
                my $txt=`tail -1 $tmpdir/dssr-2ndstrs.dbn|sed 's/&//g'`;
                $txt.=`$bindir/x3dna-dssr -i=$tmpdir/arena.pdb --pair-only`;
                open(FP,">$dssrfile");
                print FP "$txt";
                close(FP);
            }
            system("cat $tmpdir/arena.pdb |cut -c1-54|gzip - > $arenafile.gz");
            if (!-s "$tmpdir/arena.pdb")
            {
                system("$bindir/Arena $filename.gz $tmpdir/arena.pdb 4");
                system("cd $tmpdir; $bindir/x3dna-dssr -i=arena.pdb >/dev/null");
                my $txt=`tail -1 $tmpdir/dssr-2ndstrs.dbn|sed 's/&//g'`;
                $txt.=`$bindir/x3dna-dssr -i=$tmpdir/arena.pdb --pair-only`;
                open(FP,">$dssrfile");
                print FP "$txt";
                close(FP);
                system("cat $tmpdir/arena.pdb |cut -c1-54|gzip - > $arenafile.gz");
            }
            system("rm -rf $tmpdir");
        }
    }
}

exit();

