#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

#system("$bindir/cif2fasta $rootdir/pdb/data/structures/divided/mmCIF/ $rootdir/pdb/derived_data/nuc.list > $rootdir/pdb/derived_data/nuc.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d{2,}' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/rna.fasta");
system("grep -PA1 --no-group-separator '\\tRNA\\t\\d\$' $rootdir/pdb/derived_data/nuc.fasta > $rootdir/pdb/derived_data/fragment.fasta");

exit();

