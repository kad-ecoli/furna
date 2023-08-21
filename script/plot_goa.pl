#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
my $bindir = dirname(abs_path(__FILE__));
my $rootdir = dirname($bindir);

print "$rootdir/goa/is_a.csv\n";
my %isa_dict;
my %direct_dict;
foreach my $line(`cat $rootdir/goa/is_a.csv`)
{
    if ($line=~/^GO:/)
    {
        chomp($line);
        my @items=split(/\t/,$line);
        my $GOterm=$items[0];
        my $parent=$items[2];
        my $indirect=$items[3];
        $direct_dict{$GOterm}=$parent;
        if (length $indirect)
        {
            $parent.=",$indirect";
        }
        $isa_dict{$GOterm}=$parent;
    }
}

print "$rootdir/data/go2name.tsv.gz\n";
my %go2name_dict;
foreach my $line(`zcat $rootdir/data/go2name.tsv.gz`)
{
    chomp($line);
    my @items=split(/\t/,$line);
    $go2name_dict{$items[0]}=$items[2];
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

print "$rootdir/data/gosvg\n";
system("mkdir -p $rootdir/data/gosvg") if (!-d "$rootdir/data/gosvg");
my $rst=`cat -n $rootdir/data/gosvg/list`;
$txt="";
foreach my $GOterm_line(`zcat $rootdir/data/parent.tsv.gz|cut -f3-|sed 's/\\t/\\n/g'|sort|uniq`)
{
    next if ($GOterm_line!~/GO:\d+/);
    chomp($GOterm_line);
    if ($rst!~/\t$GOterm_line\n/)
    {
        $txt.="$GOterm_line\n";
    }
}
open(FP,">>$rootdir/data/gosvg/list");
print FP $txt;
close($txt);

my $maxwidth=20;
foreach my $line(`cat -n $rootdir/data/gosvg/list`)
{
    if ($line=~/(\d+)\t(\S+)/)
    {
        my $idx="$1";
        my $GOterm_line="$2";
        my $svgfile="$rootdir/data/gosvg/$idx.svg";
        next if (-s "$svgfile");
        my $dotfile="$rootdir/data/gosvg/$idx.dot";
        print "[$idx] $GOterm_line\n";

        my $GVtxt="digraph G{ graph[splines=true,rankdir=\"BT\"];\n";
        foreach my $GOterm(split(/,/,$GOterm_line))
        {
            my $label;
            my $curlen=0;
            foreach my $word(split(/ /,substr("$GOterm $go2name_dict{$GOterm}",0,100)))
            {
                if ($curlen)
                {
                    if (1+(length $word)+$curlen<=$maxwidth)
                    {
                        $label.=" $word";
                        $curlen+=1+length $word;
                    }
                    else
                    {
                        if (length $word<=$maxwidth)
                        {
                            $label.="\n$word";
                            $curlen=length $word;
                        }
                        else
                        {
                            while (length $word)
                            {
                                $label.="\n";
                                $label.=substr($word,0,$maxwidth);
                                $curlen=length substr($word,0,$maxwidth);
                                $word=substr($word,$maxwidth);
                            }
                        }
                    }
                }
                else
                {
                    $label.="$word";
                    $curlen+=length $word;
                }
            }
            
            $GVtxt.="\"$GOterm\"[label=\"$label\" shape=rectangle fillcolor=white style=filled];\n";
            foreach my $parent(split(/,/,$direct_dict{$GOterm}))
            {
                $GVtxt.="\"$GOterm\"->\"$parent\";\n";
            }
        }
        $GVtxt.="}\n";
        open(FP,">$dotfile");
        print FP "$GVtxt";
        close(FP);
        system("$rootdir/graphviz/bin/dot -Tsvg -o$svgfile $dotfile");
        system("rm $dotfile") if (-s "$svgfile");
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
