#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

system("mkdir -p $rootdir/infernal");

my @rna_list;
my %rna_nr;
my $rna;
my %rna_cm_dict;
foreach my $line(`grep '^NAME ' $rootdir/data/rna.cm|grep -ohP "\\S+\$"|uniq`)
{
    chomp($line);
    $rna_cm_dict{$line}=1;
}

foreach my $line(`sed 's/\\t/,/g' $rootdir/data/rna_nr.tsv`)
{
    chomp($line);
    my @items=split(/,/,$line);
    my $rna=$items[0];
    next if (exists($rna_cm_dict{$rna}));
    push(@rna_list,($rna));
    $rna_nr{$rna}=$line;
}

my %fasta_dict;
foreach my $line(`cat $rootdir/data/rna.fasta`)
{
    chomp($line);
    if ($line=~/^>(\S+)/)
    {
        $rna="$1";
    }
    else
    {
        if (!exists($rna_cm_dict{$rna}))
        {
            $fasta_dict{$rna}="$line";
        }
    }
}



foreach my $rna(@rna_list)
{
    next if ( -s "$rootdir/infernal/$rna.hmm" ||
             !-s "$rootdir/infernal/$rna.sto");
    system("$bindir/hmmbuild --rna $rootdir/infernal/$rna.hmm $rootdir/infernal/$rna.sto");
}

my $txt="";
foreach my $rna(@rna_list)
{
    if ( -s "$rootdir/infernal/$rna.hmm")
    {
        $txt.=`cat $rootdir/infernal/$rna.hmm`;
    }
}
open(FP,">$rootdir/data/rna.hmm");
print FP "$txt";
close(FP);
system("$bindir/hmmpress -f $rootdir/data/rna.hmm");
exit();
