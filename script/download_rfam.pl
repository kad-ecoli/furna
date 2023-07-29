#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download Rfam\n";
system("mkdir -p $rootdir/Rfam");
&download_from_rfam("Rfam.pdb.gz", "Rfam/Rfam.pdb.gz");

&download_from_rfam("rfam2go/rfam2go", "Rfam/rfam2go");

exit();

sub download_from_rfam
{
    my ($url_query,$outfile)=@_;
    system("wget -q  ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/$url_query -O $rootdir/$outfile");
    system("wget -q http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
