#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "download cisbp\n";
system("mkdir -p $rootdir/cisbp/");
#&download_from_cisbp("TF_Information_all_motifs.txt.zip", "cisbp/TF_Information_all_motifs.txt.zip");
#&download_from_cisbp("PWMs.zip", "cisbp/PWMs.zip");

#system("cd $rootdir/cisbp; unzip -o TF_Information_all_motifs.txt.zip");
#system("cd $rootdir/cisbp; unzip -o PWMs.zip");

print "$rootdir/data/taxon.tsv.gz\n";
my %species2taxon;
foreach my $line(`zcat $rootdir/data/taxon.tsv.gz`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $taxon=$items[0];
    my $species=$items[1];
    $species=~s/ /_/g;
    $species2taxon{$species}=$taxon;
}

my @motif_list;
my %motif_dict;
my $txt;
print "$rootdir/cisbp/TF_Information_all_motifs.txt\n";
foreach my $line(`cat $rootdir/cisbp/TF_Information_all_motifs.txt`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    my $species=$items[7];
    my $Motif_ID=$items[3];
    next if (!exists($species2taxon{$species}) || $Motif_ID eq ".");
    $items[7]=$species2taxon{$species};
    $line=join("\t",@items);
    if (!exists($motif_dict{$Motif_ID}))
    {
        push(@motif_list,($Motif_ID));
        $motif_dict{$Motif_ID}=$line;
    }
    $txt.="$line\n";
}

system("head -1 $rootdir/cisbp/TF_Information_all_motifs.txt > $rootdir/data/TF_Information_all_motifs.txt");
open(FP,">>$rootdir/data/TF_Information_all_motifs.txt");
print FP "$txt";
close(FP);

my %fasta_dict;
my @taxon_list;
foreach my $taxon(`cat $rootdir/data/TF_Information_all_motifs.txt |cut -f8|uniq|sort|uniq|grep -v TF_Species`)
{
    chomp($taxon);
    push(@taxon_list,($taxon));
    $fasta_dict{$taxon}="";
}
my $taxon_num=scalar @taxon_list;
print "$taxon_num taxon\n";

my %Motif2taxon;
foreach my $line(`cat $rootdir/data/TF_Information_all_motifs.txt|cut -f4,8|grep -v TF_Species`)
{
    if ($line=~/^(\S+)\t(\d+)/)
    {
        my $Motif_ID="$1";
        my $taxon="$2";
        $Motif2taxon{$Motif_ID}=$taxon;
    }
}

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
    #system("rm -rf $rootdir/genomes/$taxon/ftp.ncbi.nlm.nih.gov");
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
    next if (-s "$rootdir/genomes/$taxon.freq");
    system("zcat $rootdir/genomes/$taxon.fasta.gz | $bindir/fasta2freq - $rootdir/genomes/$taxon.freq");
}

system("mkdir -p $rootdir/cisbp/meme");
foreach my $Motif_ID(@motif_list)
{
    my $out_file="$rootdir/cisbp/meme/$Motif_ID.txt";
    my $pwm_file="$rootdir/cisbp/pwms/$Motif_ID.txt";
    next if (-s "$out_file" || `cat $pwm_file|wc -l`+0<=1);
    my $cmd="$bindir/cisbp2meme.py --pwm_file $pwm_file --out_file $out_file";
    if (exists($Motif2taxon{$Motif_ID}))
    {
        my $taxon="$Motif2taxon{$Motif_ID}";
        if (-s "$rootdir/genomes/$taxon.freq")
        {
            $cmd.=" --freq_file $rootdir/genomes/$taxon.freq";
        }
    }
    print  "$cmd\n";
    system("$cmd");
}

exit();

sub download_from_cisbp
{
    my ($url_query,$outfile)=@_;
    system("wget -q http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/$url_query -O $rootdir/$outfile") if (!-s "$rootdir/$outfile");
}
