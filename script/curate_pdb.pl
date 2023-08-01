#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

#system("$bindir/cif2fasta $rootdir/pdb/data/structures/divided/mmCIF/ $rootdir/pdb/derived_data/nuc.list > $rootdir/pdb/derived_data/nuc.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d{2,}' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/rna.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d\$' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/fragment.fasta");

system("mkdir -p $rootdir/data") if (!-d "$rootdir/data");
system("$bindir/SortFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $rootdir/pdb/derived_data/rna.fasta $rootdir/data/rna.fasta");
system("$bindir/SortFastaWithResolution.py $rootdir/pdb/derived_data/index/resolu.idx $rootdir/pdb/derived_data/fragment.fasta $rootdir/data/fragment.fasta");
system("$bindir/fasta_nr.py $rootdir/data/rna.fasta $rootdir/data/rna_nr.fasta $rootdir/data/rna_nr.tsv");

foreach my $chain (`cat $rootdir/pdb/derived_data/rna.fasta |grep '>'|cut -f1|sed 's/^>//g'`)
{
    chomp($chain);
    if ($chain=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainID="$2";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        system("mkdir -p $rootdir/chain/$divided") if (!-d "$rootdir/chain/$divided");
        my $filename="$rootdir/chain/$divided/$pdbid$chainID.pdb";
        next if (-s "$filename.gz");
        print "$filename.gz\n";
        system("$bindir/cif2chain $rootdir/pdb/data/structures/divided/mmCIF/$divided/$pdbid.cif.gz $filename $chainID");
        system("gzip -f $filename");
    }
}

exit();

