#!/usr/bin/perl
my $docstring=<<EOF
download_taxon.pl
    download NCBI taxononmy
EOF
;
use strict;
use File::Basename;
use Cwd 'abs_path';

my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);
system("mkdir -p $rootdir/taxononmy");
system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O $rootdir/taxononmy/taxononmy.tar.gz");
if (!-s "$rootdir/taxononmy/taxononmy.tar.gz")
{
    system("wget -q https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O $rootdir/taxononmy/taxononmy.tar.gz");
}

my %has_taxon;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |grep -v '^#'|cut -f12|sort|uniq`)
{
    chomp($line);
    if (length $line)
    {
        $has_taxon{$line}=1;
    }
}
foreach my $line(`zcat $rootdir/data/protein.tsv.gz |grep -v '^#'|cut -f9|sort|uniq`)
{
    chomp($line);
    if (length $line)
    {
        $has_taxon{$line}=1;
    }
}

system("cd $rootdir/taxononmy; tar -xvf  taxononmy.tar.gz names.dmp");
my @taxon_list;
my %taxon_dict;
foreach my $line(`cat $rootdir/taxononmy/names.dmp`)
{
    chomp($line);
    next if (length $line==0);
    my @items=split(/\t\|\t/,$line);
    my $taxon=$items[0];
    next if (!exists($has_taxon{$taxon}));
    my $name =$items[1];
    if (!exists($taxon_dict{$taxon}))
    {
        $taxon_dict{$taxon}=$name;
        push(@taxon_list,("$taxon"));
    }
    elsif ($line=~/\tscientific name\t\|$/)
    {
        $taxon_dict{$taxon}=$name;
    }
}
my $txt;
foreach my $taxon(@taxon_list)
{
    $txt.="$taxon\t$taxon_dict{$taxon}\n";
}
open(FP,">$rootdir/data/taxon.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/taxon.tsv");

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
