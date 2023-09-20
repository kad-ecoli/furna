#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download rbpdb\n";
system("mkdir -p $rootdir/rbpdb/");
&download_from_rbpdb("PFMDir.zip", "rbpdb/PFMDir.zip");
&download_from_rbpdb("invivo.zip", "rbpdb/invivo.zip");
&download_from_rbpdb("RBPDB_v1.3.1_2012-11-21_TDT.zip", "rbpdb/RBPDB_TDT.zip");
&download_from_rbpdb("RBPDB_v1.3.1_experiments_2012-11-21.tdt", "rbpdb/RBPDB_experiments.tdt");
&download_from_rbpdb("RBPDB_v1.3.1_proteins_2012-11-21.tdt", "rbpdb/RBPDB_proteins.tdt");
&download_from_rbpdb("RBPDB_v1.3.1_protExp_2012-11-21.tdt", "rbpdb/RBPDB_protExp.tdt");


exit();

sub download_from_rbpdb
{
    my ($url_query,$outfile)=@_;
    my $cmd="wget -q http://rbpdb.ccbr.utoronto.ca/downloads/$url_query -O $rootdir/$outfile";
    print  "$cmd\n";
    system("$cmd") if (!-s "$rootdir/$outfile");
}
