#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "curate GOA\n";
system("mkdir -p $rootdir/goa/UNIPROT");
system("zgrep '^RNAcentral' $rootdir/goa/UNIPROT/goa_uniprot_all.gaf.gz|gzip - > $rootdir/goa/goa_RNAcentral_all.gaf.gz");
system("$bindir/obo2csv $rootdir/goa/go-basic.obo $rootdir/goa/is_a.csv $rootdir/goa/name.csv $rootdir/goa/alt_id.csv 0");
system("zcat $rootdir/goa/goa_RNAcentral_all.gaf.gz | $bindir/subset_GOA - $rootdir/RNAcentral/pdb.list $rootdir/goa/goa_RNAcentral_subset.gaf"); 
system("cat $rootdir/goa/name.csv |gzip - > $rootdir/data/go2name.tsv.gz");
exit();
