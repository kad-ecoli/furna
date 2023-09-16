#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "formatdb\n";
system("cd $rootdir/data; $bindir/makeblastdb -in rna_nr.fasta -dbtype nucl");
system("cd $rootdir/data; $bindir/makeblastdb -in dna_nr.fasta -dbtype nucl");
system("cd $rootdir/data; $bindir/makeblastdb -in protein_nr.fasta -dbtype prot");

exit();
