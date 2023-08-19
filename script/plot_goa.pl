#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "get parent GO term\n";
my %isa_dict;
foreach my $line(`cat $rootdir/goa/is_a.csv`)
{
    if ($line=~/^GO:/)
    {
        chomp($line);
        my @items=split(/\t/,$line);
        my $GOterm=$items[0];
        my $parent=$items[2];
        my $indirect=$items[3];
        if (length $indirect)
        {
            $parent.=",$indirect";
        }
        $isa_dict{$GOterm}=$parent;
    }
}

my $txt="";
print "$rootdir/data/rna.tsv.gz\n";
foreach my $line(`zcat $rootdir/data/rna.tsv.gz |cut -f1,2,9,10,11|grep -F "GO:"`)
{
    chomp($line);
    my @items=split(/\t/,$line);

    for (my $i=2;$i<=4;$i++)
    {
        next if (length $items[$i]==0);
        my %GOterm_dict;
        my @GOterm_list;
        foreach my $GOterm(split(/,/,$items[$i]))
        {
            next if (exists($GOterm_dict{$GOterm}));
            $GOterm_dict{$GOterm}=1;
            push(@GOterm_list,($GOterm));
            if (exists($isa_dict{$GOterm}))
            {
                foreach my $parent(split(",",$isa_dict{$GOterm}))
                {
                    next if (exists($GOterm_dict{$parent}));
                    $GOterm_dict{$parent}=1;
                    push(@GOterm_list,($parent));
                }
            }
        }
        my @sort_list = sort @GOterm_list;
        $items[$i]=join(",",@sort_list);
    }
    $txt.="$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]\n";
}
print "$rootdir/data/protein.tsv.gz\n";
foreach my $line(`zcat $rootdir/data/protein.tsv.gz |cut -f1,2,6,7,8|grep -F "GO:"`)
{
    chomp($line);
    my @items=split(/\t/,$line);

    for (my $i=2;$i<=4;$i++)
    {
        next if (length $items[$i]==0);
        my %GOterm_dict;
        my @GOterm_list;
        foreach my $GOterm(split(/,/,$items[$i]))
        {
            next if (exists($GOterm_dict{$GOterm}));
            $GOterm_dict{$GOterm}=1;
            push(@GOterm_list,($GOterm));
            if (exists($isa_dict{$GOterm}))
            {
                foreach my $parent(split(",",$isa_dict{$GOterm}))
                {
                    next if (exists($GOterm_dict{$parent}));
                    $GOterm_dict{$parent}=1;
                    push(@GOterm_list,($parent));
                }
            }
        }
        my @sort_list = sort @GOterm_list;
        $items[$i]=join(",",@sort_list);
    }
    $txt.="$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]\n";
}


open(FP,">$rootdir/data/parent.tsv");
print FP "$txt";
close(FP);

&gzipFile("$rootdir/data/parent.tsv");

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
