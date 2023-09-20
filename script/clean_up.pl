#!/usr/bin/perl
# clean up intermediate files
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

foreach my $filename(`find $rootdir/output/*.gz -mmin +60 2>/dev/null`)
{
    chomp($filename);
    my $cmd="rm -f $filename";
    system("$cmd");
}
foreach my $filename(`find $rootdir/output/*.fsearch/aln.m8 -mmin +60 2>/dev/null |sed 's/\\/aln.m8//g'`)
{
    chomp($filename);
    my $cmd="rm -rf $filename*";
    system("$cmd");
}


exit();
