#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);
my $tmpdir="$rootdir/tmp";
system("mkdir -p $tmpdir");

my @ribocentre_list;
my %ribocentre2lbs;
my %rfam2ribocentre;
my %ribocentre2rfam;
foreach my $line(`grep -v '^#' $rootdir/ribocentre/docs.tsv`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $ribocentre=$items[0];
    my $rfam_list =$items[2];
    my $lbs_list  =$items[4];
    push(@ribocentre_list,($ribocentre));
    foreach my $rfam(split(/,/,$rfam_list))
    {
        $rfam2ribocentre{$rfam}=$ribocentre;
    }
    $ribocentre2lbs{$ribocentre}=$lbs_list;
    $ribocentre2rfam{$ribocentre}=$rfam_list;
}

my $match_txt="#pdb\tchain\tRfam\tTM1\tTM2\tribocentre\tcsa_origin\tcsa_renum\n";
my %match_dict;
foreach my $line(`zcat $rootdir/data/ribocentre.tsv.gz|grep -v '^#'|sort|uniq`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $pdbid=$items[0];
    my $chainid=$items[1];
    my $divided=substr($pdbid,length($pdbid)-3,2);
    my $key="$divided/$pdbid$chainid";
    $match_dict{$key}=1;
    $match_txt.="$line\n";
}

foreach my $ribocentre(@ribocentre_list)
{
    print "$ribocentre\n";
    #next if ("$ribocentre" eq "Ribosome");
    #next if ("$ribocentre" ne "VS-ribozyme");
    my $rfam_list=$ribocentre2rfam{$ribocentre};
    #next if ($rfam_list!~/NA/);
    my $lbs_list =$ribocentre2lbs{$ribocentre};
    my $list_txt;
    my %target_dict;
    my %target2rfam;
    foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2,5|grep -v '^#'|sed 's/\\t\$/\tNA/g'`)
    {
        chomp($line);
        my @items=split(/\t/,$line);
        my $pdbid=$items[0];
        my $chainid=$items[1];
        my $divided=substr($pdbid,length($pdbid)-3,2);
        my $key="$divided/$pdbid$chainid";
        next if (exists($match_dict{$key}));
        foreach my $rfam(split(/,/,$items[2]))
        {
            if ($rfam_list=~/$rfam/)
            {
                $list_txt.="$key\n";
                $target_dict{$key}="$pdbid:$chainid";
                $target2rfam{$key}=$rfam;
                last;
            }
        }
    }
    open(FP,">$tmpdir/${ribocentre}.chain.list");
    print FP "$list_txt";
    close(FP);
    if (!-s "$tmpdir/${ribocentre}.chain.list")
    {
        print "No hit for $ribocentre\n";
        next;
    }
    my $cmd="$bindir/USalign $rootdir/ribocentre/${ribocentre}.ent -dir2 $rootdir/arena/ $tmpdir/${ribocentre}.chain.list -suffix .pdb.gz -outfmt 2 -fast |grep -v '^#'|cut -f2,3,4 |grep -P \"(\\t1\\.0000)|(\\t0\\.[4-9]\\d+)\" > $tmpdir/${ribocentre}.TM.tsv";
    if ($rfam_list eq "NA")
    {
        $cmd="$bindir/USalign $rootdir/ribocentre/${ribocentre}.ent -dir2 $rootdir/arena/ $tmpdir/${ribocentre}.chain.list -suffix .pdb.gz -outfmt 2 -fast |grep -v '^#'|cut -f2,3,4 |grep -P \"(\\t1\\.0000)|(\\t0\\.[4-9]\\d+\\t0\\.[4-9]\\d+)\" > $tmpdir/${ribocentre}.TM.tsv";
    }
    print  "$cmd\n";
    system("$cmd");
    system("cat $tmpdir/${ribocentre}.TM.tsv |cut -f1 -d. > $tmpdir/${ribocentre}.TM.list");
    if (!-s "$tmpdir/${ribocentre}.TM.list")
    {
        print "No hit for $ribocentre\n";
        next;
    }
    my $cmd2="cd $rootdir/ribocentre; $bindir/USalign ${ribocentre}.ent -dir2 $rootdir/arena/ $tmpdir/${ribocentre}.TM.list -suffix .pdb.gz -outfmt 1 > $tmpdir/${ribocentre}.TM.txt";
    if ($ribocentre eq "Ribosome")
    {
        $cmd2="cd $rootdir/ribocentre; $bindir/USalign ${ribocentre}.ent -dir2 $rootdir/arena/ $tmpdir/${ribocentre}.TM.list -suffix .pdb.gz -outfmt 1 -fast > $tmpdir/${ribocentre}.TM.txt";
    }
    print  "$cmd2\n";
    system("$cmd2");

    my $rst=`cat $tmpdir/${ribocentre}.TM.txt`;
    foreach my $block(split(/\$\$\$\$/,$rst))
    {
        my $aln1="";
        my $aln2="";
        my $tm1="";
        my $tm2="";
        my $key="";
        foreach my $line(split(/\n/,$block))
        {
            if ($line=~/^>$ribocentre.ent:\S*\t\S+\t\S+\t\S+\tTM-score=([.\d]+)/)
            {
                $tm1="$1";
            }
            elsif ($line=~/^>(\w+\/\w+).pdb.gz:\S*\t\S+\t\S+\t\S+\tTM-score=([.\d]+)/)
            {
                $key="$1";
                $tm2="$2";
            }
            elsif ($line=~/^([-A-Za-z]+)$/)
            {
                if (length($aln1)==0)
                {
                    $aln1="$1";
                }
                else
                {
                    $aln2="$1";
                }
            }
        }
        next if (length($aln1)*length($tm1)*length($aln2)*length($tm2)*length($key)==0 || !exists($target_dict{$key}));
        next if ($tm1<0.45 && $tm2<0.45);
        next if ($target2rfam{$key} eq "NA" && ($tm1<0.45 || $tm2<0.45));
        my $lbs_renum="";
        my $r1=0;
        my $r2=0;
        for (my $r=0;$r<length($aln1);$r++)
        {
            my $nt1=substr($aln1,$r,1);
            my $nt2=substr($aln2,$r,1);
            $r1+=($nt1 ne '-');
            $r2+=($nt2 ne '-');
            if ($nt1 ne '-' && $nt2 ne '-' && $lbs_list=~/\b$r1\b/)
            {
                if (length($lbs_renum)==0)
                {
                    $lbs_renum="$r2";
                }
                else
                {
                    $lbs_renum.=",$r2";
                }
            }
        }
        next if (length($lbs_renum)==0);
        my $pdbid="";
        my $chainid="";
        if ($target_dict{$key}=~/(\w+):(\w+)/)
        {
            $pdbid="$1";
            $chainid="$2";
        }
        my $lbs_origin="";
        my $filename="$rootdir/arena/$key.pdb.gz";
        foreach my $line(`zcat $rootdir/arena/$key.pdb.gz|grep -F " C3'"|cut -c23-27|cat -n`)
        {
            if ($line=~/(\d+)\s+(\w+)/)
            {
                my $r="$1";
                my $resi="$2";
                if ($lbs_renum=~/\b$r\b/)
                {
                    if (length($lbs_origin)==0)
                    {
                        $lbs_origin="$resi";
                    }
                    else
                    {
                        $lbs_origin.=",$resi";
                    }
                }
            }
        }
        
        $match_txt.="$pdbid\t$chainid\t$target2rfam{$key}\t$tm1\t$tm2\t$ribocentre\t$lbs_origin\t$lbs_renum\n";
        #if ($key eq "ff/1ffk0")
        #{
            #print ">$ribocentre\ttm1=$tm1\tlbs_list=$lbs_list\n$aln1\n";
            #print ">$key\ttm2=$tm2\tlbs_renum=$lbs_renum\n$aln2\n";
        #}
    }
}
open(FP,">$rootdir/data/ribocentre.tsv");
print FP "$match_txt";
close(FP);
&gzipFile("$rootdir/data/ribocentre.tsv");

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
