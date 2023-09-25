#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download attract\n";
system("mkdir -p $rootdir/attract/logo/");
system("mkdir -p $rootdir/data/attract/");
system("wget -q https://attract.cnic.es/attract/static/ATtRACT.zip -O $rootdir/attract/ATtRACT.zip");
system("cd $rootdir/attract; unzip -f ATtRACT.zip");
if (!-f "$rootdir/attract/pwm.txt")
{
    system("cd $rootdir/attract; unzip ATtRACT.zip");
}
foreach my $line(`grep '^>' $rootdir/attract/pwm.txt`)
{
    if ($line=~/>(\S+)/)
    {
        my $MotifID="$1";
        print "$MotifID\n";
        my $pngfile="$MotifID";
        $pngfile=~s/\./_/g;
        $pngfile="logo$pngfile.png";
        my $pngfile1="$pngfile";
        $pngfile1=~s/_pfm/.pfm/g;
        if (! -s "$rootdir/data/attract/$pngfile")
        {
            my $cmd="wget https://attract.cnic.es/attract/static/logo/$pngfile1 -O $rootdir/data/attract/$pngfile";
            print  "$cmd\n";
            system("$cmd");
        }
    }
}

my %taxon_dict;
foreach my $taxon(`zcat $rootdir/data/rna.tsv.gz |grep -v '^#'|cut -f12|sort|uniq`)
{
    chomp($taxon);
    $taxon_dict{$taxon}=1;
}

my %species2taxon;
foreach my $line(`zcat $rootdir/data/taxon.tsv.gz`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $taxon="$items[0]";
    next if (!exists($taxon_dict{$taxon}));
    my $species="$items[1]";
    $species=~s/ /_/g;
    $species=~s/\[//g;
    $species=~s/\]//g;
    $species2taxon{$species}=$taxon;
}

my @taxon_list;
my $taxon_txt="";
foreach my $species(`cut -f4 $rootdir/attract/ATtRACT_db.txt|sort|uniq|grep -v "^Organism\$"`)
{
    chomp($species);
    if (!exists($species2taxon{$species}))
    {
        print "$species not in $rootdir/data/taxon.tsv.gz\n";
        next;
    }
    my $taxon=$species2taxon{$species};
    push(@taxon_list,($taxon));

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
    if (`ls $rootdir/genomes/$taxon/*_rna.fna.gz|wc -l`+0==0)
    {
        $taxon_txt.="$taxon\n";
        next;
    }
    my $filename=`ls -sS $rootdir/genomes/$taxon/*_rna.fna.gz| head -1|grep -ohP "\\S+\$"`;
    chomp($filename);
    system("mv $filename $rootdir/genomes/$taxon.fasta.gz");
    system("rm -rf $rootdir/genomes/$taxon");
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
    }

    next if (-s "$rootdir/genomes/$taxon.freq");
    system("zcat $rootdir/genomes/$taxon.fasta.gz | $bindir/fasta2freq - $rootdir/genomes/$taxon.freq");
}


my %Motif2taxon;
foreach my $line(`cat $rootdir/attract/ATtRACT_db.txt |cut -f4,12|sort|uniq`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $species=$items[0];
    my $MotifID=$items[1];
    if (!exists($species2taxon{$species}))
    {
        next;
    }
    $Motif2taxon{$MotifID}=$species2taxon{$species};
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

my $rst="\n";
$rst.=`cat $rootdir/attract/pwm.txt`;
my @blocks=split(/\n>/,$rst);
for (my $b=1; $b<scalar @blocks; $b++)
{
    my $MotifID="";
    my $nsites ="";
    my @lines=split(/\n/,$blocks[$b]);
    if ($lines[0]=~/(\S+)\t(\d+)/)
    {
        $MotifID="$1";
        $nsites ="$2";
    }
    next if (!exists($Motif2taxon{$MotifID}));
    my $taxon=$Motif2taxon{$MotifID};
    my $txt  ="";
    for (my $i=1;$i<scalar @lines;$i++)
    {
        $txt.="$lines[$i]\n";
    }
    $txt=~s/\t/  /g;

    $meme_dict{$taxon}.="\nMOTIF $MotifID\nletter-probability matrix: alength= 4 w= $nsites nsites= $nsites\n$txt";
}

system("mkdir -p $rootdir/attract/fimo/");
my %fasta_dict;
foreach my $taxon(@taxon_list)
{
    open(FP,">$rootdir/attract/fimo/$taxon.meme.txt");
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
    if (!-s "$rootdir/genomes/$taxon.bg")
    {
        system("zcat $rootdir/genomes/$taxon.fasta.gz > $rootdir/genomes/$taxon.fasta");
        system("$bindir/fasta-get-markov -norc $rootdir/genomes/$taxon.fasta $rootdir/genomes/$taxon.bg");
        system("rm $rootdir/genomes/$taxon.fasta");
    }

    open(FP,">$rootdir/attract/fimo/$taxon.fasta");
    print FP $fasta_dict{$taxon};
    close(FP);
    my $cmd="$bindir/fimo --norc --bfile $rootdir/genomes/$taxon.bg -o $rootdir/attract/fimo/$taxon $rootdir/attract/fimo/$taxon.meme.txt $rootdir/attract/fimo/$taxon.fasta";
    print  "$cmd\n";
    system("$cmd");
    $cmd="cat $rootdir/attract/fimo/$taxon/fimo.tsv |grep -vP '\\t\\d+\\t\\d+\\t-\\t' > $rootdir/attract/fimo/$taxon.tsv";
    print  "$cmd\n";
    system("$cmd");
    if (-s "$rootdir/attract/fimo/$taxon.tsv")
    {
        system("rm -rf $rootdir/attract/fimo/$taxon");
    }
}


my %attract_dict;
foreach my $line(`cat $rootdir/attract/ATtRACT_db.txt |cut -f1,2,9,12|sort|uniq|grep -v "^Gene_name"`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $Gene_name="$items[0]";
    my $Gene_id  ="$items[1]";
    my $Pubmed   ="$items[2]";
    my $Motif_ID ="$items[3]";
    $attract_dict{$Motif_ID}="\t$Gene_name\t$Gene_id\t$Pubmed";
}


my $txt="";
foreach my $taxon(@taxon_list)
{
    foreach my $line(`cat $rootdir/attract/fimo/$taxon.tsv|cut -f1,3,4,5,9,10|sort -k2|grep -vP '^#'`)
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
            $txt.="$pdbid\t$chainid\t$start\t$end\t$qvalue\t$sequence\t$Motif_ID$attract_dict{$Motif_ID}\n";
        }
    }
}
open(FP,">$rootdir/data/attract_fimo.tsv");
print FP $txt;
close(FP);
&gzipFile("$rootdir/data/attract_fimo.tsv");

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
