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

foreach my $pmid(`grep -v '#' $rootdir/data/rna.tsv |cut -f7|sort|uniq|sed 's/,/\\n/g'|sort|uniq`)
{
    chomp($pmid);
    next if (-s "$pubmeddir/$pmid.txt");
    &download_pubmed($pmid,$pubmeddir);
}

exit();

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
