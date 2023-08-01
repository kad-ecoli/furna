#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "curate Rfam\n";
system("mkdir -p $rootdir/Rfam");
system("rm $rootdir/Rfam/rna_nr.split.*");
system("split -l 1000 $rootdir/data/rna_nr.fasta $rootdir/Rfam/rna_nr.split.");
system("ls $rootdir/Rfam/|grep rna_nr.split. |cut -f3 -d.|sort|uniq > $rootdir/Rfam/rna_nr.list");
foreach my $batch (`tac $rootdir/Rfam/rna_nr.list`)
{
    chomp($batch);
    my $cmd="$bindir/cmsearch --cpu 4 -Z 1 --toponly --tblout $rootdir/Rfam/rna_nr.tblout.$batch $rootdir/Rfam/Rfam.cm.gz $rootdir/Rfam/rna_nr.split.$batch > /dev/null";
    print  "$cmd\n";
    system("$cmd");
}
system("cat $rootdir/Rfam/rna_nr.tblout.* > $rootdir/Rfam/rna_nr.tblout");

if (-s "$rootdir/Rfam/rna_nr.tblout")
{
    foreach my $batch (`tac $rootdir/Rfam/rna_nr.list`)
    {
        system("rm $rootdir/Rfam/rna_nr.tblout.$batch $rootdir/Rfam/rna_nr.split.$batch");
    }
}

exit();
