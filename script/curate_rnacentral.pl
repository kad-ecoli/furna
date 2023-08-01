#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "curate RNAcentral\n";
system("mkdir -p $rootdir/RNAcentral");
system("cat $rootdir/RNAcentral/pdb.tsv |cut -f1,4|sort|uniq|sed 's/\\t/_/g' > $rootdir/RNAcentral/pdb.list");
system("zcat $rootdir/RNAcentral/rnacentral_species_specific_ids.fasta.gz |cut -f1 -d' ' |  $bindir/fasta2miss $rootdir/RNAcentral/pdb.list - - $rootdir/RNAcentral/pdb.fasta");
exit();
