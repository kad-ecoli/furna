#!/usr/bin/perl
my $docstring=<<EOF
download_pubmed.pl pmid
    download pubmed abstract to pmid.txt
EOF
;
use strict;
use File::Basename;
use Cwd 'abs_path';

if (scalar @ARGV)
{
    my $pmid=$ARGV[0];
    my $pubmeddir=".";
    $pubmeddir=$ARGV[1] if (scalar @ARGV>1);
    &download_pubmed($ARGV[$a],$pubmeddir);
    exit(0);
}


my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);
my $pubmeddir = "$rootdir/pubmed";
system("mkdir -p $pubmeddir");

my $cmd="cat $rootdir/data/rna.tsv | grep -v '#'  |cut -f7|sort|uniq|sed 's/,/\\n/g'|sort|uniq";
if (-s "$rootdir/data/rna.tsv.gz")
{
    $cmd="z$cmd";
}
my $txt;
foreach my $pmid(`$cmd`)
{
    chomp($pmid);
    if (!-s "$pubmeddir/$pmid.txt")
    {
        &download_pubmed($pmid,$pubmeddir);
    }
    if (-s "$pubmeddir/$pmid.txt")
    {
        my $title=`head -1 $pubmeddir/$pmid.txt`;
        $txt.="$pmid\t$title";
    }
}
open(FP,">$rootdir/data/pubmed.tsv");
print FP "$txt";
close(FP);

&gzipFile("$rootdir/data/pubmed.tsv");

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

sub download_pubmed
{
    my ($pmid,$pubmeddir)=@_;

    return if (-s "$pubmeddir/$pmid.txt");
    my $abstract_txt="";
    my $title_txt="";
    my $ab=0;
    my $ti=0;
    my $cmd="curl -s 'https://pubmed.ncbi.nlm.nih.gov/$pmid/?format=pubmed'";
    print "$cmd\n";
    foreach my $line(`$cmd`)
    {
        $line=~s/\r//g;
        chomp($line);
        if ($line=~/^AB  - /)
        {
            $ab=1;
            $line=substr($line, 6);
            $abstract_txt=$line;
        }
        elsif ($line=~/^TI  - /)
        {
            $ti=1;
            $line=substr($line, 6);
            $title_txt=$line;
        }
        elsif ($line=~/^\w+\s+- /)
        {
            $ab=2;
            $ti=2;
        }
        elsif ($ab==1)
        {
            $line=substr($line, 6);
            $abstract_txt.="$line";
        }
        elsif ($ti==1)
        {
            $line=substr($line, 6);
            $title_txt.="$line";
        }
    }
    if (length $abstract_txt || length $title_txt)
    {
        open(FP,">$pubmeddir/$pmid.txt");
        print FP "$title_txt\n$abstract_txt\n";
        close(FP);
    }
}
