#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download RNAcentral\n";
system("mkdir -p $rootdir/RNAcentral");
&download_from_rnacentral("id_mapping/database_mappings/pdb.tsv",
                                           "RNAcentral/pdb.tsv");

&download_from_rnacentral("sequences/rnacentral_species_specific_ids.fasta.gz",
                         "RNAcentral/rnacentral_species_specific_ids.fasta.gz");
exit();

sub download_from_rnacentral
{
    my ($url_query,$outfile)=@_;
    system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/$url_query -O $rootdir/$outfile");
    system("wget -q http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
