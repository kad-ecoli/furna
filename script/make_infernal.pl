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
foreach my $line(`sed 's/\\t/,/g' $rootdir/data/rna_nr.tsv`)
{
    chomp($line);
    my @items=split(/,/,$line);
    my $rna=$items[0];
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
        $fasta_dict{$rna}="$line";
    }
}

my %cssr_dict;
foreach my $rna(@rna_list)
{
    my $pdbid=substr($rna,0,4);
    if (length $rna>8 && substr($rna,0,8)==lc(substr($rna,0,8)))
    {
        $pdbid=substr($rna,0,8);
    }
    my $divided=substr($pdbid,(length $pdbid)-3,2);
    my $filename="$rootdir/cssr/$divided/$rna.cssr";
    my $cssr=`cat $filename|head -1`;
    chomp($cssr);
    if (length $cssr==0)
    {
        system("$bindir/CSSR $rootdir/chain/$divided/$rna.pdb.gz $filename -o 3");
        $cssr=`cat $filename|head -1`;
        chomp($cssr);
    }
    $cssr=~s/[\[\]\(\)]/./g;
    $cssr_dict{$rna}=$cssr;
}

foreach my $rna(@rna_list)
{
    next if ( -s "$rootdir/infernal/$rna.cm");
    my $sequence=$fasta_dict{$rna};
    my $cssr    =$cssr_dict{$rna};
    my $L1=length $sequence;
    my $L2=length $cssr;
    if ($L1 != $L2)
    {
        print "$rna\t$L1!=$L2\n";
        next;
    }
    open(FP,">$rootdir/infernal/$rna.sto");
    print FP "# STOCKHOLM 1.0\n\n";
    my $line="                  ";
    $line=$rna.(substr($line,length $rna));
    print FP "$line$sequence\n";
    print FP "#=GC SS_cons      $cssr\n//\n";
    close(FP);
    open(FP,">$rootdir/infernal/$rna.sh\n");
    print FP "#!/bin/bash\n";
    print FP "#SBATCH -o $rootdir/infernal/$rna.out\n";
    print FP "#SBATCH -t 48:00:00\n";
    print FP "#SBATCH -p sigbio\n";
    if ($L2<100) { print FP "#SBATCH --mem=1gb\n"; }
    else         { print FP "#SBATCH --mem=".($L2*10)."mb\n"; }
    print FP "cd $rootdir/infernal/\n";
    print FP "$bindir/cmbuild $rootdir/infernal/$rna.cmbuild $rootdir/infernal/$rna.sto\n";
    print FP "$bindir/cmcalibrate --cpu 1 $rootdir/infernal/$rna.cmbuild\n";
    print FP "mv $rootdir/infernal/$rna.cmbuild $rootdir/infernal/$rna.cm\n";
    close(FP);
    system("chmod +x $rootdir/infernal/$rna.sh");
    system("$rootdir/infernal/$rna.sh");
}

my $txt="";
foreach my $rna(@rna_list)
{
    if ( -s "$rootdir/infernal/$rna.cm")
    {
        $txt.=`cat $rootdir/infernal/$rna.cm`;
    }
}
open(FP,">$rootdir/data/rna.cm");
print FP "$txt";
close(FP);
system("$bindir/cmpress -F $rootdir/data/rna.cm");
system("$bindir/make_hmmer.pl");
exit();
