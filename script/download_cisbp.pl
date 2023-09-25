#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download cisbp-RNA\n";
system("mkdir -p $rootdir/cisbp/");
my $url="http://cisbp-rna.ccbr.utoronto.ca/tmp/entiredata_2023_09_25_9:33_am.zip";
my $cmd="wget -q '$url' -O $rootdir/cisbp/entiredata.zip";
print  "$cmd\n";
system("$cmd");
system("cd $rootdir/cisbp/; unzip -f entiredata.zip");
if (!-s "$rootdir/cisbp/RBP_Information_all_motifs.txt")
{
    system("cd $rootdir/cisbp/; unzip entiredata.zip");
}

print "$rootdir/data/rna.tsv.gz\n";
my %taxon_dict;
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |grep -v '^#'|cut -f12|sort|uniq`)
{
    chomp($line);
    $taxon_dict{$line}=1;
}
my $taxon_num=scalar %taxon_dict;
print "$taxon_num species\n";


print "$rootdir/data/taxon.tsv.gz\n";
my %species2taxon;
foreach my $line(`zcat $rootdir/data/taxon.tsv.gz`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $taxon=$items[0];
    next if (!exists($taxon_dict{$taxon}));
    my $species=$items[1];
    $species=~s/ /_/g;
    $species2taxon{$species}=$taxon;
}
$taxon_num=scalar %species2taxon;
print "$taxon_num species\n";

my @motif_list;
my %motif_dict;
my $txt;
print "$rootdir/cisbp/RBP_Information_all_motifs.txt\n";
foreach my $line(`cat $rootdir/cisbp/RBP_Information_all_motifs.txt`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $species=$items[7];
    my $Motif_ID=$items[3];
    next if (!exists($species2taxon{$species}) || $Motif_ID eq ".");
    my $taxon=$species2taxon{$species};
    $items[7]=$taxon;
    $line=join("\t",@items);
    if (!exists($motif_dict{$Motif_ID}))
    {
        push(@motif_list,($Motif_ID));
        $motif_dict{$Motif_ID}=$taxon;
    }
    $txt.="$line\n";
}

print "$rootdir/data/RBP_Information_all_motifs.txt\n";
system("head -1 $rootdir/cisbp/RBP_Information_all_motifs.txt > $rootdir/data/RBP_Information_all_motifs.txt");
open(FP,">>$rootdir/data/RBP_Information_all_motifs.txt");
print FP "$txt";
close(FP);

my @taxon_list;
foreach my $taxon(`cat $rootdir/data/RBP_Information_all_motifs.txt |cut -f8|uniq|sort|uniq|grep -v _Species`)
{
    chomp($taxon);
    if (exists($taxon_dict{$taxon}))
    {
        push(@taxon_list,($taxon));
    }
}
$taxon_num=scalar @taxon_list;
print "$taxon_num taxon\n";

my $motif_num=scalar %motif_dict;
print "$motif_num motif\n";

foreach my $taxon(@taxon_list)
{
    if (-s "$rootdir/genomes/$taxon.fasta.gz")
    {
        if (-d "$rootdir/genomes/$taxon")
        {
            system("rm -r $rootdir/genomes/$taxon");
        }
        next;
    }
    system("mkdir -p $rootdir/genomes/$taxon");
    system("cd $rootdir/genomes/$taxon; wget -r --no-parent -A '*_rna.fna.gz' ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/$taxon/");
    foreach my $filename(`find $rootdir/genomes/$taxon/ftp.ncbi.nlm.nih.gov -type f|grep _rna.fna.gz`)
    {
        chomp($filename);
        system("mv $filename  $rootdir/genomes/$taxon/");
    }
    system("rm -rf $rootdir/genomes/$taxon/ftp.ncbi.nlm.nih.gov");
    next if (`ls $rootdir/genomes/$taxon/*_rna.fna.gz|wc -l`+0==0);
    my $filename=`ls -sS $rootdir/genomes/$taxon/*_rna.fna.gz| head -1|grep -ohP "\\S+\$"`;
    chomp($filename);
    system("mv $filename $rootdir/genomes/$taxon.fasta.gz");
    system("rm -rf $rootdir/genomes/$taxon");
}


my $taxon_txt="";
foreach my $taxon(@taxon_list)
{
    if (!-s "$rootdir/genomes/$taxon.fasta.gz")
    {
        $taxon_txt.="$taxon\n";
    }
}
open(FP,">$rootdir/RNAcentral/genomes.list");
print FP $taxon_txt;
close(FP);
if (length $taxon_txt)
{
    system("zcat $rootdir/RNAcentral/rnacentral_species_specific_ids.fasta.gz|$bindir/fasta2taxon $rootdir/RNAcentral/genomes.list - - $rootdir/RNAcentral/genome.fasta");
}
foreach my $taxon(@taxon_list)
{
    if (!-s "$rootdir/genomes/$taxon.fasta.gz")
    {
        system("grep -PA1  --no-group-separator '_$taxon ' RNAcentral/genome.fasta|gzip - > $rootdir/genomes/$taxon.fasta.gz");
        $taxon_txt.="$taxon\n";
    }
}

foreach my $taxon(@taxon_list)
{
    print "$rootdir/genomes/$taxon.freq\n";
    next if (-s "$rootdir/genomes/$taxon.freq" && -s "$rootdir/genomes/$taxon.bg");
    system("zcat $rootdir/genomes/$taxon.fasta.gz | $bindir/fasta2freq - $rootdir/genomes/$taxon.freq");

    system("zcat $rootdir/genomes/$taxon.fasta.gz > $rootdir/genomes/$taxon.fasta");
    system("$bindir/fasta-get-markov -norc $rootdir/genomes/$taxon.fasta $rootdir/genomes/$taxon.bg");
    system("rm $rootdir/genomes/$taxon.fasta");
}


my %meme_dict;
foreach my $taxon(@taxon_list)
{
    $meme_dict{$taxon}=<<EOF
MEME version 4

ALPHABET= ACGT
strands: +

Background letter frequencies
EOF
;
    if (-s "$rootdir/genomes/$taxon.freq")
    {
        $meme_dict{$taxon}.=`cat $rootdir/genomes/$taxon.freq`;
    }
    else
    {
        $meme_dict{$taxon}.="A 0.25 C 0.25 G 0.25 T 0.25\n";
    }
}

foreach my $Motif_ID(@motif_list)
{
    my $taxon=$motif_dict{$Motif_ID};
    print "$Motif_ID => $taxon\n";
    if (!-s "$rootdir/cisbp/pwms/$Motif_ID.txt")
    {
        my $cmd="wget -q http://cisbp-rna.ccbr.utoronto.ca/data/0.6/DataFiles/PWMs/Files/$Motif_ID.txt -O $rootdir/cisbp/pwms/$Motif_ID.txt";
        print "$cmd\n";
        system("$cmd");
        next if (!-s "$rootdir/cisbp/pwms/$Motif_ID.txt");
    }
    $meme_dict{$taxon}.="\n".`$bindir/cisbp2meme $rootdir/cisbp/pwms/$Motif_ID.txt -|grep -A1000 '^MOTIF '`;
}

system("mkdir -p $rootdir/cisbp/fimo/");
my %fasta_dict;
foreach my $taxon(@taxon_list)
{
    open(FP,">$rootdir/cisbp/fimo/$taxon.meme.txt");
    print FP $meme_dict{$taxon};
    close(FP);
    $fasta_dict{$taxon}="";
}
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2,12,13`)
{
    my @items=split(/\t/,$line);
    my $taxon=$items[2];
    my $sequence=$items[3];
    $sequence=~s/u/t/g;
    $sequence=~s/i/a/g;
    $sequence=uc($sequence);
    chomp($sequence);
    $fasta_dict{$taxon}.=">$items[0]_$items[1]\n$sequence\n";
}



foreach my $taxon(@taxon_list)
{
    open(FP,">$rootdir/cisbp/fimo/$taxon.fasta");
    print FP $fasta_dict{$taxon};
    close(FP);
    my $cmd="$bindir/fimo --norc --bfile $rootdir/genomes/$taxon.bg -o $rootdir/cisbp/fimo/$taxon $rootdir/cisbp/fimo/$taxon.meme.txt $rootdir/cisbp/fimo/$taxon.fasta";
    print  "$cmd\n";
    system("$cmd");
    $cmd="cat $rootdir/cisbp/fimo/$taxon/fimo.tsv |grep -vP '\\t\\d+\\t\\d+\\t-\\t' > $rootdir/cisbp/fimo/$taxon.tsv";
    print  "$cmd\n";
    system("$cmd");
    if (-s "$rootdir/cisbp/fimo/$taxon.tsv")
    {
        system("rm -rf $rootdir/cisbp/fimo/$taxon");
    }
}

my %cisbp_dict;
foreach my $line(`cat $rootdir/data/RBP_Information_all_motifs.txt |cut -f4,6,7,20|sort|uniq`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $Motif_ID ="$items[0]";
    my $Gene_id  ="$items[1]";
    my $Gene_name="$items[2]";
    my $Pubmed   ="$items[3]";
    $cisbp_dict{$Motif_ID}="\t$Gene_name\t$Gene_id\t$Pubmed";
}

my $txt="";
foreach my $taxon(@taxon_list)
{
    foreach my $line(`cat $rootdir/cisbp/fimo/$taxon.tsv|cut -f1,3,4,5,9,10|sort -k2|grep -vP '^#'`)
    {
        if ($line=~/^(\S+)\t(\w+)_(\w+)\t(\d+)\t(\d+)\t(\S+)\t(\w+)/)
        {
            my $Motif_ID="$1";
            my $pdbid   ="$2";
            my $chainid ="$3";
            my $start   ="$4";
            my $end     ="$5";
            my $qvalue  ="$6";
            my $sequence="$7";
            $txt.="$pdbid\t$chainid\t$start\t$end\t$qvalue\t$sequence\t$Motif_ID$cisbp_dict{$Motif_ID}\n";
        }
    }
}
open(FP,">$rootdir/data/cisbp_fimo.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/cisbp_fimo.tsv");


system("mkdir -p $rootdir/data/cisbp");
foreach my $Motif_ID(`zcat $rootdir/data/cisbp_fimo.tsv.gz|cut -f7|sort|uniq`)
{
    chomp($Motif_ID);
    my $pngfile="logo$Motif_ID";
    $pngfile=~s/\./_/g;
    $pngfile="$rootdir/data/cisbp/$pngfile.png";
    if (!-s "$pngfile")
    {
        system("cp $rootdir/cisbp/logos/${Motif_ID}_fwd.png $pngfile")
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

