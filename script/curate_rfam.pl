#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "curate Rfam\n";
system("mkdir -p $rootdir/Rfam");

my @rfam_list;
my %rfam_dict;
my $rfam;
foreach my $line(`zcat $rootdir/Rfam/Rfam.cm.gz |grep -P '(^ACC )|(^DESC )'`)
{
    chomp($line);
    if ($line=~/^ACC\s+(RF\d+)/)
    {
        $rfam="$1";
        push(@rfam_list, ($rfam)) if (!exists($rfam_dict{$rfam}));
    }
    elsif ($line=~/^DESC\s+([\s\S]+)/)
    {
        $rfam_dict{$rfam}="$1";
    }
}
my $txt;
foreach my $rfam(@rfam_list)
{
    $txt.="$rfam\t$rfam_dict{$rfam}\n";
}
open(FP,">$rootdir/data/rfam.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/rfam.tsv");

system("rm $rootdir/Rfam/rna_nr.split.*");
system("zcat $rootdir/Rfam/Rfam.pdb.gz | cut -f2,3|sed 's/\t//g'|sort|uniq | $bindir/fasta2miss - $rootdir/data/rna.fasta - $rootdir/Rfam/Rfam.pdb.fasta");
system("cat $rootdir/data/rna.fasta|grep '>'|sed 's/>//g' | $bindir/fasta2miss - $rootdir/Rfam/Rfam.pdb.fasta $rootdir/Rfam/Rfam.pdb.miss.list - > /dev/null");
system("$bindir/fasta2len.py $rootdir/data/rna.fasta|grep -P '\\t\\d{1,3}\$' |cut -f1 >> $rootdir/Rfam/Rfam.pdb.miss.list");
system("cat $rootdir/Rfam/Rfam.pdb.miss.list|sort|uniq | $bindir/fasta2miss - $rootdir/data/rna.fasta - $rootdir/Rfam/rna.fasta");
system("$bindir/fasta_nr.py $rootdir/Rfam/rna.fasta $rootdir/Rfam/rna_nr.fasta $rootdir/Rfam/rna_nr.tsv");
system("rm $rootdir/Rfam/rna.fasta  $rootdir/Rfam/Rfam.pdb.miss.list  $rootdir/Rfam/Rfam.pdb.fasta");
system("split -l 1024 $rootdir/Rfam/rna_nr.fasta $rootdir/Rfam/rna_nr.split.");

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
        system("rm $rootdir/Rfam/rna_nr.tblout.$batch");
        system("rm $rootdir/Rfam/rna_nr.split.$batch");
    }
}

exit();

sub gzipFile
{
    my ($filename)=@_;
    my $oldNum=`zcat $filename.gz 2>/dev/null|wc -l`+0;
    my $newNum=` cat $filename   |wc -l`+0;
    if (0.8*$oldNum>$newNum)
    {
        print "WARNING! do not update $filename from $oldNum to $newNum entries\n";
        return;
    }
    print "update $filename from $oldNum to $newNum entries\n";
    system("gzip -f $filename");
    return;
}
