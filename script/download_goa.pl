#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download GOA\n";
system("mkdir -p $rootdir/goa/UNIPROT");
&download_from_goa("UNIPROT/goa_uniprot_all.gaf.gz", "goa/UNIPROT/goa_uniprot_all.gaf.gz");

exit();

sub download_from_goa
{
    my ($url_query,$outfile)=@_;
    system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/$url_query -O $rootdir/$outfile");
    system("wget -q http://ftp.ebi.ac.uk/pub/databases/GO/goa/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
