#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "combine chain level annotation\n";
my @target_list;
my %len_dict;
print "$rootdir/pdb/derived_data/rna.fasta\n";
foreach my $line(`grep '^>' $rootdir/pdb/derived_data/rna.fasta`)
{
    if ($line=~/^>(\w+:\w+)\tRNA\t(\d+)/)
    {
        my $target="$1";
        my $L="$2";
        push(@target_list, ($target));
        $target=~s/://g;
        $len_dict{$target}=$L;
    }
}

my %resolu_dict;
print "$rootdir/pdb/derived_data/index/resolu.idx\n";
foreach my $line(`cat $rootdir/pdb/derived_data/index/resolu.idx`)
{
    if ($line=~/^(\w+)\t;\t([-.\d]+)/)
    {
        my $pdbid=lc("$1");
        my $resolu="$2";
        $resolu="NA" if ($resolu<0);
        $resolu_dict{$pdbid}=$resolu;
    }
}

my %title_dict;
print "$rootdir/pdb/derived_data/index/cmpd_res.idx\n";
foreach my $line(`cat $rootdir/pdb/derived_data/index/cmpd_res.idx`)
{
    chomp($line);
    if ($line=~/^(\w+)\t;\t\S*\t;\t([\s\S]+)/)
    {
        my $pdbid=lc("$1");
        my $title="$2";
        $title_dict{$pdbid}=$title;
    }
}

my %pubmed_dict;
print "$rootdir/sifts/pdb_pubmed.tsv.gz\n";
foreach my $line(`zcat $rootdir/sifts/pdb_pubmed.tsv.gz`)
{
    if ($line=~/^(\w+)\t\d+\t(\d+)/)
    {
        my $pdbid=lc("$1");
        my $pubmed="$2";
        if (exists ($pubmed_dict{$pdbid}))
        {
            $pubmed_dict{$pdbid}.=",$pubmed";
        }
        else
        {
            $pubmed_dict{$pdbid}="$pubmed";
        }
    }
}

my %fasta_dict;
my $target;
print "$rootdir/data/rna.fasta\n";
foreach my $line(`cat $rootdir/data/rna.fasta`)
{
    chomp($line);
    next if (length $line==0);
    if ($line=~/^>(\w+)/)
    {
        $target="$1";
    }
    else
    {
        $fasta_dict{$target}="$line";
    }
}

my %rna_nr_dict;
print "$rootdir/data/rna_nr.tsv\n";
foreach my $line(`cat $rootdir/data/rna_nr.tsv`)
{
    chomp($line);
    if ($line=~/^(\w+)\t/)
    {
        my $target="$1";
        $rna_nr_dict{$target}=$line;
    }
}
my %rna_nr_rfam_dict;
print "$rootdir/data/rna_nr_rfam.tsv\n";
foreach my $line(`cat $rootdir/Rfam/rna_nr.tsv`)
{
    chomp($line);
    if ($line=~/^(\w+)\t/)
    {
        my $target="$1";
        $rna_nr_rfam_dict{$target}=$line;
    }
}

my %alt_id_dict;
print "$rootdir/goa/alt_id.csv\n";
foreach my $line(`cat $rootdir/goa/alt_id.csv`)
{
    if ($line=~/(GO:\d+)\t(GO:\d+)/)
    {
        my $oldGOterm="$1";
        my $newGOterm="$2";
        if ($oldGOterm ne $newGOterm)
        {
            $alt_id_dict{$oldGOterm}=$newGOterm;
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

my %rfam2go_dict;
print "$rootdir/Rfam/rfam2go\n";
foreach my $line(`cat $rootdir/Rfam/rfam2go`)
{
    chomp($line);
    if ($line=~/^Rfam:(RF\d+)[\s\S]+(GO:\d+)$/)
    {
        my $rfam="$1";
        my $GOterm="$2";
        if (exists($alt_id_dict{$GOterm}))
        {
            $GOterm=$alt_id_dict{$GOterm};
        }
        if (!exists($aspect_dict{$GOterm}))
        {
            print "WARNING! rfam2go: no such term $GOterm\n";
            next;
        }
        if (exists($rfam2go_dict{$rfam}))
        {
            $rfam2go_dict{$rfam}.=",$GOterm";
        }
        else
        {
            $rfam2go_dict{$rfam}="$GOterm";
        }
    }
}

my %rfam_dict;
my %rfam_name_dict;
print "$rootdir/Rfam/rna_nr.tblout\n";
foreach my $line(`grep -v '^#' $rootdir/Rfam/rna_nr.tblout`)
{
    if ($line=~/^(\w+)\s+\S+\s+(\S+)\s+(RF\d+)/)
    {
        my $target="$1";
        my $name  ="$2";
        my $rfam  ="$3";
        $rfam_name_dict{$rfam}="$name";
        if (exists($rfam_dict{$target}))
        {
            $rfam_dict{$target}.=",$rfam";
        }
        else
        {
            $rfam_dict{$target}="$rfam";
        }
        if (exists($rna_nr_rfam_dict{$target}))
        {
            foreach my $member(split(/\t/,$rna_nr_rfam_dict{$target}))
            {
                if ($member ne $target)
                {
                    $rfam_dict{$member}=$rfam_dict{$target};
                }
            }
        }
    }
}
my %rfam_origin_dict;
my @rfam_origin_list;
print "$rootdir/Rfam/Rfam.pdb.gz";
foreach my $line(`zcat $rootdir/Rfam/Rfam.pdb.gz|cut -f1,2,3`)
{
    if ($line=~/(RF\d+)\t(\w+)\t(\w+)/)
    {
        my $rfam="$1";
        my $pdbid="$2";
        my $chainid="$3";
        my $target="$pdbid$chainid";
        if (exists($rfam_origin_dict{$target}))
        {
            $rfam_origin_dict{$target}.=",$rfam";
        }
        else
        {
            $rfam_origin_dict{$target}="$rfam";
            push(@rfam_origin_list,($target));
        }
    }
}
foreach my $target(@rfam_origin_list)
{
    if (!exists($rfam_dict{$target}))
    {
        $rfam_dict{$target}=$rfam_origin_dict{$target};
    }
}


my %rnacentral_dict;
my %taxon_dict;
print "$rootdir/RNAcentral/pdb.tsv\n";
foreach my $line(`cat $rootdir/RNAcentral/pdb.tsv`)
{
    if ($line=~/^(\w+)\tPDB\t(\w+)_(\w+)\t(\d+)/)
    {
        my $rnacentral="$1";
        my $pdbid  =lc("$2");
        my $chainid   ="$3";
        my $taxon     ="$4";
        my $chain="$pdbid$chainid";
        $rnacentral.="_$taxon";
        if (exists($rnacentral_dict{$chain}))
        {
            $rnacentral_dict{$chain}.=",$rnacentral";
        }
        else
        {
            $rnacentral_dict{$chain} ="$rnacentral";
        }
        if (exists($taxon_dict{$chain}))
        {
            $taxon_dict{$chain}.=",$taxon";
        }
        else
        {
            $taxon_dict{$chain}="$taxon";
        }
    }
}

my %goa_dict;
print "$rootdir/goa/goa_RNAcentral_subset.gaf\n";
foreach my $line(`cat $rootdir/goa/goa_RNAcentral_subset.gaf`)
{
    if ($line=~/^RNAcentral\t(\S+)\t\S+\t(\S*)\t(GO:\d+)/)
    {
        my $rnacentral="$1";
        my $qualifier="$2";
        my $GOterm="$3";
        next if ($qualifier=~/NOT/);
        if (exists($alt_id_dict{$GOterm}))
        {
            $GOterm=$alt_id_dict{$GOterm};
        }
        if (!exists($aspect_dict{$GOterm}))
        {
            print "WARNING! goa: no such term $GOterm\n";
            next;
        }
        if (exists($goa_dict{$rnacentral}))
        {
            $goa_dict{$rnacentral}.=",$GOterm";
        }
        else
        {
            $goa_dict{$rnacentral}="$GOterm";
        }
    }
}

print "$rootdir/chain/*/*.pdb.gz\n";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid  ="$1";
        my $chainid="$2";
        my $chain  ="$pdbid$chainid";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        next if (exists($taxon_dict{$chain}));
        
        my $filename="$rootdir/chain/$divided/$chain.pdb.gz";
        foreach my $line(`zcat $filename|grep -F 'ORGANISM_TAXID:'`)
        {
            if ($line=~/ORGANISM_TAXID: (\d+);/)
            {
                my $taxon="$1";
                if (exists($taxon_dict{$chain}))
                {
                    $taxon_dict{$chain}.=",$taxon";
                }
                else
                {
                    $taxon_dict{$chain}="$taxon";
                }
                print "$filename $taxon\n";
            }
        }
    }
}

# map taxon, rfam and RNAcentral to identical sequence
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid  ="$1";
        my $chainid="$2";
        my $target ="$pdbid$chainid";
        next if (!exists($rna_nr_dict{$target}));
        my @member_list=split(/\t/,$rna_nr_dict{$target});

        my $has_rfam=0;
        my $has_taxon=0;
        my $has_rnacentral=0;
        my $taxon_list="";
        my $rfam_list="";
        my $rnacentral_list="";
        foreach my $member(@member_list)
        {
            if (exists($rfam_dict{$member}))
            {
                $has_rfam++;
                $rfam_list=$rfam_dict{$member} if (length $rfam_list==0);
            }
            if (exists($taxon_dict{$member}))
            {
                $has_taxon++;
                foreach my $taxon(split(/,/,$taxon_dict{$member}))
                {
                    next if ($taxon eq $taxon_list);
                    if ($taxon_list!~/\b$taxon\b/)
                    {
                        $taxon_list.=",$taxon";
                    }
                }
            }
            if (exists($rnacentral_dict{$member}))
            {
                $has_rnacentral++;
                foreach my $rnacentral(split(/,/,$rnacentral_dict{$member}))
                {
                    next if ($rnacentral eq $rnacentral_list);
                    if ($rnacentral_list!~/\b$rnacentral\b/)
                    {
                        $rnacentral_list.=",$rnacentral";
                    }
                }
            }
        }
        my $member_num=scalar @member_list;
        if (0<$has_rfam && $has_rfam<scalar @member_list)
        {
            foreach my $member(@member_list)
            {
                if (!exists($rfam_dict{$member}))
                {
                    $rfam_dict{$member}=$rfam_list;
                }
            }
        }
        if (0<$has_taxon && $has_taxon<scalar @member_list)
        {
            my $taxon=substr($taxon_list,1);
            if ($taxon!~/,/)
            {
                foreach my $member(@member_list)
                {
                    if (!exists($taxon_dict{$member}))
                    {
                        $taxon_dict{$member}=$taxon;
                        print "$member => $taxon\n";
                    }
                }
            }
        }
        if (0<$has_rnacentral && $has_rnacentral<scalar @member_list)
        {
            my $rnacentral=substr($rnacentral_list,1);
            if ($rnacentral!~/,/)
            {
                foreach my $member(@member_list)
                {
                    if (!exists($rnacentral_dict{$member}))
                    {
                        $rnacentral_dict{$member}=$rnacentral;
                        print "$member => $rnacentral\n";
                    }
                }
            }
        }
    }
}

print "$rootdir/cssr/*/*.cssr\n";
print "$rootdir/dssr/*/*.dssr\n";
my %cssr_dict;
my %dssr_dict;
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid  ="$1";
        my $chainid="$2";
        my $chain  ="$pdbid$chainid";
        my $divided=substr($pdbid,length($pdbid)-3,2);
        my $cssrfile="$rootdir/cssr/$divided/$chain.cssr";
        if (-s "$cssrfile")
        {
            my $cssr=`head -1 $cssrfile`;
            chomp($cssr);
            $cssr_dict{$chain}=$cssr;
        }
        my $dssrfile="$rootdir/dssr/$divided/$chain.dssr";
        if (-s "$dssrfile")
        {
            my $dssr=`head -1 $dssrfile`;
            chomp($dssr);
            $dssr_dict{$chain}=$dssr;
        }
    }
}

print "$rootdir/goa/ec2go\n";
my %go2ec_dict;
foreach my $line(`cat $rootdir/goa/ec2go`)
{
    chomp($line);
    if ($line=~/^(EC:[-.\d]+)/)
    {
        my $ec="$1";
        if ($line=~/(GO:\d+)$/)
        {
            my $GOterm="$1";
            $go2ec_dict{$GOterm}=$ec;
        }
    }
}

print "$rootdir/goa/is_a.csv\n";
my %isa_dict;
foreach my $line(`cat $rootdir/goa/is_a.csv`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $GOterm =$items[0];
    my $Aspect =$items[1];
    my $parent =$items[2];
    my $indirect=$items[3];
    if (length $parent && length $indirect)
    {
        $parent.=",$indirect";
    }
    $isa_dict{$GOterm}=$parent;
}


print "$rootdir/pdb/derived_data/pdb_seqres.txt\n";
my %name_dict;
foreach my $line(`grep '>' $rootdir/pdb/derived_data/pdb_seqres.txt |grep -F mol:na`)
{
    chomp($line);
    if ($line=~/>(\w+)_(\w+)\s+mol:na\s+length:\d+\s+([\s\S]+)/)
    {
        my $pdbid="$1";
        my $chainid="$2";
        my $name="$3";
        my $target="$pdbid:$chainid";
        $name_dict{$target}=$name;
    }
}

print "$rootdir/data/rna.tsv\n";
my $txt="#pdb\tchain\tL\tresolution\tRfam\tRNAcentral\tpubmed\t";
$txt.="EC\tGO_MF\tGO_BP\tGO_CC\ttaxon\tsequence\tcssr\tdssr\ttitle\tname\n";
foreach my $target(@target_list)
{
    if ($target=~/(\w+):(\w+)/)
    {
        my $pdbid  ="$1";
        my $chainid="$2";
        my $chain  ="$pdbid$chainid";
        
        my $resolu ="";
        my $rfam   ="";
        my $rnacentral="";
        my $pubmed ="";
        my $EC_list="";
        my $GO_MF  ="";
        my $GO_BP  ="";
        my $GO_CC  ="";
        my $taxon  ="";
        my $cssr   ="";
        my $dssr   ="";
        my $title  ="";
        my $name="";

        if (exists($resolu_dict{$pdbid}))
        {
            $resolu="$resolu_dict{$pdbid}";
        }
        if (exists($rfam_dict{$chain}))
        {
            my $rfam_list=$rfam_dict{$chain};
            foreach my $rfam(split(/,/,$rfam_list))
            {
                if (exists($rfam2go_dict{$rfam}))
                {
                    my $GOterm_list=$rfam2go_dict{$rfam};
                    foreach my $GOterm(split(/,/,$GOterm_list))
                    {
                        my $Aspect=$aspect_dict{$GOterm};
                        if ($Aspect=~/F/ && $GO_MF!~/$GOterm/)
                        {
                            $GO_MF.=",$GOterm";
                        }
                        elsif ($Aspect=~/P/ && $GO_BP!~/$GOterm/)
                        {
                            $GO_BP.=",$GOterm";
                        }
                        elsif ($Aspect=~/C/ && $GO_CC!~/$GOterm/)
                        {
                            $GO_CC.=",$GOterm";
                        }
                    }
                }
            }
            $rfam="$rfam_list";
        }
        if (exists($rnacentral_dict{$chain}))
        {
            my $rnacentral_list=$rnacentral_dict{$chain};
            foreach my $rnacentral(split(/,/,$rnacentral_list))
            {
                if (exists($goa_dict{$rnacentral}))
                {
                    my $GOterm_list=$goa_dict{$rnacentral};
                    foreach my $GOterm(split(/,/,$GOterm_list))
                    {
                        my $Aspect=$aspect_dict{$GOterm};
                        if ($Aspect=~/F/ && $GO_MF!~/$GOterm/)
                        {
                            $GO_MF.=",$GOterm";
                        }
                        elsif ($Aspect=~/P/ && $GO_BP!~/$GOterm/)
                        {
                            $GO_BP.=",$GOterm";
                        }
                        elsif ($Aspect=~/C/ && $GO_CC!~/$GOterm/)
                        {
                            $GO_CC.=",$GOterm";
                        }
                    }
                }
            }
            $rnacentral="$rnacentral_list";
        }
        if (exists($pubmed_dict{$pdbid}))
        {
            $pubmed=$pubmed_dict{$pdbid};
        }
        if (exists($taxon_dict{$chain}))
        {
            $taxon=$taxon_dict{$chain};
        }
        if (exists($cssr_dict{$chain}))
        {
            $cssr=$cssr_dict{$chain};
        }
        if (exists($dssr_dict{$chain}))
        {
            $dssr=$dssr_dict{$chain};
        }
        if (exists($title_dict{$pdbid}))
        {
            $title="$title_dict{$pdbid}";
        }
        if (exists($name_dict{$target}))
        {
            $name=$name_dict{$target};
        }
        
        my $sequence=$fasta_dict{$chain};
        my $L=length $sequence;

        $GO_MF=substr($GO_MF,1) if ($GO_MF=~/^,/);
        $GO_BP=substr($GO_BP,1) if ($GO_BP=~/^,/);
        $GO_CC=substr($GO_CC,1) if ($GO_CC=~/^,/);

        foreach my $GOterm(split(/,/,$GO_MF))
        {
            if (exists($go2ec_dict{$GOterm}))
            {
                my $EC=$go2ec_dict{$GOterm};
                next if ($EC_list=~/$EC/);
                if ($EC=~/(EC:[.\d]+).[-.]+$/)
                {
                    my $prefix="$1";
                    next if ($EC_list=~/$prefix/);
                }
                $EC_list.=",$EC";
            }
        }
        if (length $EC_list==0)
        {
            foreach my $GOterm(split(/,/,$GO_MF))
            {
                if (!exists($go2ec_dict{$GOterm}) && exists($isa_dict{$GOterm}))
                {
                    foreach my $parent(split(/,/,$isa_dict{$GOterm}))
                    {
                        next if (!exists($go2ec_dict{$parent}));
                        my $EC=$go2ec_dict{$parent};
                        next if ($EC_list=~/$EC/);
                        if ($EC=~/(EC:[.\d]+).[-.]+$/)
                        {
                            my $prefix="$1";
                            next if ($EC_list=~/$prefix/);
                        }
                        $EC_list.=",$EC";
                    }
                }
            }
        }
        $EC_list=substr($EC_list,1) if ($EC_list=~/^,/);

        $txt.="$pdbid\t$chainid\t$L\t$resolu\t$rfam\t$rnacentral\t$pubmed\t$EC_list\t$GO_MF\t$GO_BP\t$GO_CC\t$taxon\t$sequence\t$cssr\t$dssr\t$title\t$name\n";
    }
}
open(FP,">$rootdir/data/rna.tsv");
print FP "$txt";
close(FP);

&gzipFile("$rootdir/data/rna.tsv");

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
