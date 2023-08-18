#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "combine protein level annotation\n";
my @target_list;
my %target_dict;
print "$rootdir/data/interaction.tsv.gz\n";
foreach my $line(`zcat $rootdir/data/interaction.tsv.gz |cut -f1,4,5,10|grep -P "\\tprotein\\t"|uniq|sort|uniq`)
{
    if ($line=~/^(\w+)\tprotein\t([-\w]+)\t(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence="$3";
        if ($chainid=~/(\w+)\-\w+/)
        {
            $chainid="$1";
        }
        my $target="$pdbid:$chainid";
        if (!exists($target_dict{$target}))
        {
            push(@target_list,($target));
            $target_dict{$target}=$sequence;
        }
    }
}

my %ec_dict;
print "$rootdir/sifts/pdb_chain_enzyme.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_chain_enzyme.tsv.gz|cut -f1,2,4|uniq`)
{
    if ($line=~/^(\w+)\t(\w+)\t([.\w]+)$/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $ec="EC:$3";
        my $target="$pdbid:$chainid";
        next if (!exists($target_dict{$target}));
        if (exists($ec_dict{$target}))
        {
            $ec_dict{$target}.=",$ec";
        }
        else
        {
            $ec_dict{$target}="$ec";
        }
    }
}

my %go_dict;
print "$rootdir/sifts/pdb_chain_go.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_chain_go.tsv.gz|cut -f1,2,6|sort|uniq`)
{
    if ($line=~/^(\w+)\t(\w+)\t(GO:\d+)$/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $go="$3";
        my $target="$pdbid:$chainid";
        next if (!exists($target_dict{$target}));
        if (exists($go_dict{$target}))
        {
            $go_dict{$target}.=",$go";
        }
        else
        {
            $go_dict{$target}="$go";
        }
    }
}

my %taxon_dict;
print "$rootdir/sifts/pdb_chain_taxonomy.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_chain_taxonomy.tsv.gz |cut -f1,2,3|uniq`)
{
    if ($line=~/^(\w+)\t(\w+)\t(\w+)$/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $taxon="$3";
        my $target="$pdbid:$chainid";
        next if (!exists($target_dict{$target}));
        if (exists($taxon_dict{$target}))
        {
            $taxon_dict{$target}.=",$taxon";
        }
        else
        {
            $taxon_dict{$target}="$taxon";
        }
    }
}

my %uniprot_dict;
print "$rootdir/sifts/pdb_chain_uniprot.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_chain_uniprot.tsv.gz |cut -f1,2,3|uniq`)
{
    if ($line=~/^(\w+)\t(\w+)\t(\w+)$/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $uniprot="$3";
        my $target="$pdbid:$chainid";
        next if (!exists($target_dict{$target}));
        if (exists($uniprot_dict{$target}))
        {
            $uniprot_dict{$target}.=",$uniprot";
        }
        else
        {
            $uniprot_dict{$target}="$uniprot";
        }
    }
}

my %aspect_dict;
print "$rootdir/goa/name.csv\n";
foreach my $line(`cat $rootdir/goa/name.csv`)
{
    if ($line=~/^(GO:\d+)\t(\w+)/)
    {
        my $GOterm="$1";
        my $Aspect="$2";
        $aspect_dict{$GOterm}=$Aspect;
    }
}

my $txt="#pdb\tchain\tL\tuniprot\tEC\tGO_MF\tGO_BP\tGO_CC\ttaxon\tsequence\n";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence=$target_dict{$target};
        my $L=length $sequence;
        my $uniprot="";
        my $ec="";
        my $go_mf="";
        my $go_bp="";
        my $go_cc="";
        my $taxon="";
        if (exists($uniprot_dict{$target}))
        {
            $uniprot=$uniprot_dict{$target};
        }
        if (exists($ec_dict{$target}))
        {
            $ec=$ec_dict{$target};
        }
        if (exists($go_dict{$target}))
        {
            foreach my $GOterm(split(/,/,$go_dict{$target}))
            {
                next if (!exists($aspect_dict{$GOterm}));
                my $Aspect=$aspect_dict{$GOterm};

                if ($Aspect=~/F/ && $go_mf!~/$GOterm/)
                {
                    $go_mf.=",$GOterm";
                }
                elsif ($Aspect=~/P/ && $go_bp!~/$GOterm/)
                {
                    $go_bp.=",$GOterm";
                }
                elsif ($Aspect=~/C/ && $go_cc!~/$GOterm/)
                {
                    $go_cc.=",$GOterm";
                }
            }
            if (length $go_mf)
            {
                $go_mf=substr($go_mf,1);
            }
            if (length $go_bp)
            {
                $go_bp=substr($go_bp,1);
            }
            if (length $go_cc)
            {
                $go_cc=substr($go_cc,1);
            }
        }
        if (exists($taxon_dict{$target}))
        {
            $taxon=$taxon_dict{$target};
        }
        $txt.="$pdbid\t$chainid\t$L\t$uniprot\t$ec\t$go_mf\t$go_bp\t$go_cc\t$taxon\t$sequence\n";
    }
}
open(FP,">$rootdir/data/protein.tsv");
print FP "$txt";
close(FP);
&gzipFile("$rootdir/data/protein.tsv");

$txt="";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $sequence=$target_dict{$target};
        $txt.=">$pdbid$chainid\n$sequence\n";
    }
}
open(FP,">$rootdir/data/protein.fasta");
print FP "$txt";
close(FP);
system("$bindir/fasta_nr.py $rootdir/data/protein.fasta $rootdir/data/protein_nr.fasta $rootdir/data/protein_nr.tsv");

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
