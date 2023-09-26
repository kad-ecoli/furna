#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download PDB\n";
system("mkdir -p $rootdir/pdb/derived_data/index/");
system("mkdir -p $rootdir/pdb/data/monomers/");
&download_from_pdb("derived_data/pdb_seqres.txt",
               "pdb/derived_data/pdb_seqres.txt");

&download_from_pdb("data/monomers/components.cif.gz",
               "pdb/data/monomers/components.cif.gz");

&download_from_pdb("derived_data/pdb_entry_type.txt",
               "pdb/derived_data/pdb_entry_type.txt");

&download_from_pdb("derived_data/index/resolu.idx",
               "pdb/derived_data/index/resolu.idx");

&download_from_pdb("derived_data/index/cmpd_res.idx",
               "pdb/derived_data/index/cmpd_res.idx");

system("cat $rootdir/pdb/derived_data/pdb_entry_type.txt |grep -P 'nuc\\t'|cut -f1 > $rootdir/pdb/derived_data/nuc.list");
foreach my $pdbid(`cat $rootdir/pdb/derived_data/nuc.list`)
{
    chomp($pdbid);
    my $divided=substr($pdbid,length($pdbid)-3,2);
    next if (-s "$rootdir/pdb/data/structures/divided/mmCIF/$divided/$pdbid.cif.gz");
    print "$pdbid\n";
    system("mkdir -p $rootdir/pdb/data/structures/divided/mmCIF/$divided"
          ) if (!-d "$rootdir/pdb/data/structures/divided/mmCIF/$divided");
    &download_from_pdb("data/structures/divided/mmCIF/$divided/$pdbid.cif.gz",
                   "pdb/data/structures/divided/mmCIF/$divided/$pdbid.cif.gz");
}

exit();

sub download_from_pdb
{
    my ($url_query,$outfile)=@_;
    system("wget -q http://files.wwpdb.org/pub/pdb/$url_query -O $rootdir/$outfile");
    system("wget -q https://files.rcsb.org/pub/pdb/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
