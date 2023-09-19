#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download uniprot\n";
my %uniprot_dict;
foreach my $line(`zcat $rootdir/data/uniprot.tsv.gz`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $uniprot=$items[0];
    my $name   =$items[1];
    $uniprot_dict{$uniprot}=$name;
}


my $txt;
foreach my $uniprot(`zcat $rootdir/data/protein.tsv.gz |grep -v '^#'|cut -f4|sed 's/,/\\n/g'|sort|uniq`)
{
    chomp($uniprot);
    if (exists($uniprot_dict{$uniprot}))
    {
        $txt.="$uniprot\t$uniprot_dict{$uniprot}\n";
    }
    else
    {
        print "$uniprot\n";
        my $rst=`curl -s https://rest.uniprot.org/uniprotkb/$uniprot.fasta|grep '^>'|sed 's/^>//g'`;
        chomp($rst);
        $txt.="$uniprot\t$rst\n" if (length $rst);
    }
}
open(FP,">$rootdir/data/uniprot.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/uniprot.tsv");

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
