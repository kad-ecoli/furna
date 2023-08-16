#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download sifts\n";
system("mkdir -p $rootdir/sifts");
&download_from_sifts("pdb_pubmed.tsv.gz",         "sifts/pdb_pubmed.tsv.gz");
&download_from_sifts("pdb_chain_uniprot.tsv.gz",  "sifts/pdb_chain_uniprot.tsv.gz");
&download_from_sifts("pdb_chain_taxonomy.tsv.gz", "sifts/pdb_chain_taxonomy.tsv.gz");
&download_from_sifts("pdb_chain_go.tsv.gz",       "sifts/pdb_chain_go.tsv.gz");
&download_from_sifts("pdb_chain_enzyme.tsv.gz",   "sifts/pdb_chain_enzyme.tsv.gz");

exit();

sub download_from_sifts
{
    my ($url_query,$outfile)=@_;
    system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/$url_query -O $rootdir/$outfile");
    system("wget -q http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
