#! /usr/bin/env perl

use strict;
use Carp;
use Data::Dumper;
# use gjoseqlib;
# use gjocolorlib;
use Getopt::Long;
# use DT;
use JSON;

my $usage = <<"End_of_Usage";

usage: dnadiff_genomes [ optioins ] -f contig-files

       -l list_file      - use tab-delimited file to specify contig files: [ name location abbrev ]
       -o out_dir        - output directory (D = ./out/)
       -g names          - list of genome names
       -j out_json       - json output
       -r                - reuse computed data from prior runs
}

End_of_Usage

my ($help, $list, $outdir, @names, $outjson, $reuse);

GetOptions("h|help"   => \$help,
           "l|list=s" => \$list,
           "o|out=s"  => \$outdir,
           "g=s"      => \@names,
           "j|json=s" => \$outjson,
           "r|reuse"  => \$reuse
          );

$help and die $usage;

my @contig_files;
if ($list && -s $list) {
    @contig_files = map { chomp; my ($name, $file, $abbr) = split /\t/; [ $name, $file, $abbr ] } grep { ! /^#/ }  `cat $list`;
} else {
    @contig_files = map { chomp; my $name = $_; $name =~ s|.*/||; $name =~ s/\.(fa|fasta)//; [ $name, $_ ] } @ARGV;
}

# print STDERR '\@names = '. Dumper(\@names);
# print STDERR '\@files = '. Dumper(\@contig_files);


my @sets;
my $no = 0;
for (@contig_files) {
    my ($name, $file, $abbr) = @$_;
    $abbr ||= ++$no;
    # print join("\t", $name, $file, $abbr) . "\n";
    push @sets, [ $name, $file, $abbr ];
}

$outdir ||= 'out';
run("mkdir -p $outdir");

my @hdrs = ("Sample", 'Ref \\ Qry', map { $_->[2] } @sets);
my @rows;

my @objs;

for my $i (0 .. $#sets) {
    my @c = ($sets[$i]->[0], $sets[$i]->[2]);
    for my $j (0 .. $#sets) {
        my $ii = $i+1;
        my $jj = $j+1;
        my $qname = $names[$i] || $sets[$i]->[0];
        my $rname = $names[$j] || $sets[$j]->[0];
        print STDERR "\rComparing genomes: #$ii vs #$jj (Q=$qname, R=$rname) ...  ";
        my $prefix = "$outdir/$i-$j";
        my $report = "$prefix.report";
        run("dnadiff -p $prefix  $sets[$i]->[1] $sets[$j]->[1] >$prefix.log 2>&1") if ! $reuse || ! -s $report;

        my $diff  = dnadiff_report_to_diff($report);
        my $stat  = dnadiff_report_to_key_stats($report);
        my $hover = dnadiff_report_to_mouseover($report);
        my $text  = dnadiff_report_text($report);
        my $sim   = 1 - $diff;

        # hue, saturation, brightness
        my $mindiff = 0.00009;
        $diff = $mindiff if $diff <= $mindiff;
        my $hue = (1 - $diff ** 0.07) / (1 - $mindiff ** 0.07);
        # my $hue = 1-(log(0.05+$diff) + 3)/3.05;
        # my $hue = 1-(log(0.1+$diff) + 2.31)/2.5;
        # my $hue = 1-(log(0.2+$diff) + 0.19)/1.8;
        # my $color = gjocolorlib::rgb2html( gjocolorlib::hsb2rgb( (1-$hue)*200/360+$diff/4,
        #                                                          $hue/3+(1-$diff)/5+0.36,
        #                                                          1 )
        #                                  );
        # $color = "#fffdd8" if $diff == 0;
        print STDERR "$diff $hue\n";
        print STDERR "$text\n";

        # my $cell = qq(<span style="background-color: $color;"> $stat </span>);

        # $stat = DT::span_mouseover($cell, DT::span_css("strains", "mouseover"), $hover);

        push @c, $stat;

        my $json_obj = { query => $qname, reference => $rname, report => $text, similarity => $sim };
        push @objs, $json_obj;
    }
    print STDERR "\n";
    push @rows, \@c;
}

# DT::print_dynamic_table(\@hdrs, \@rows, { title => 'Strain Comparison' });

if ($outjson) {
    print STDERR $outjson."\n";
    open(F, ">$outjson") or die "Could not open $outjson";
    print F JSON::encode_json(\@objs);
    close(F);
}

print STDERR "done\n";

sub dnadiff_report_to_diff {
    my ($report_file) = @_;
    my ($alignedBases) = `grep "^AlignedBases" $report_file` =~ /\(([0-9.]+)%\)/;
    my ($avgIdentity)  = [split(/\s+/, `grep "^AvgIdentity"  $report_file |head -n 1`)]->[1];

    my $sim = sprintf("%.4f", $alignedBases/100 * $avgIdentity/100);
    return sprintf("%.4f", 1 - $sim);
}

sub dnadiff_report_to_key_stats {
    my ($report_file) = @_;
    my ($alignedBases) = `grep "^AlignedBases" $report_file` =~ /\(([0-9.]+)%\)/;
    my ($avgIdentity)  = [split(/\s+/, `grep "^AvgIdentity"  $report_file |head -n 1`)]->[1];

    my $sim = sprintf("%.4f", $alignedBases/100 * $avgIdentity/100);

    # $alignedBases = sprintf "%5g", $alignedBases;
    # $avgIdentity = sprintf "%5g", $avgIdentity;

    # my $stat =  "$avgIdentity / $alignedBases";
    # $stat =~ s/^\s+//; $stat =~ s/\s+$//;
    # return $stat;
}

sub dnadiff_report_text {
    my ($report_file) = @_;
    return `cat $report_file |head -n41 |tail -n39`;
}

sub dnadiff_report_to_mouseover {
    my ($report_file) = @_;
    my $first = `head -n 1 $report_file`; chomp($first);
    my ($c1, $c2) = split(/\s+/, $first);
    $c1 =~ s|.*/||; $c1 =~ s/\.(fa|fna|fasta)$//; $c1 =~ s/(_filtered|)_scaffolds//; $c1 =~ s/Mycobacterium_[a-z]+_//; $c1 =~ s/Subsample_//;
    $c2 =~ s|.*/||; $c2 =~ s/\.(fa|fna|fasta)$//; $c2 =~ s/(_filtered|)_scaffolds//; $c2 =~ s/Mycobacterium_[a-z]+_//; $c2 =~ s/Subsample_//;

    my @lines = `cat $report_file |head -n41 |tail -n39`;

    my @rows;
    push @rows, [ undef, $c1, $c2 ];
    for (@lines) {
        chomp;
        my @cols = /^\[/ ? ("<b>".$_."<\/b>") : split/\s+/;
        push @rows, [$cols[0], $cols[1], $cols[2]];
    }
    my $hover = "<table class=mouseoverTable>";
    for my $row (@rows) {
        $hover .= " <tr>";
        $hover .= join('', map { "<td>$_</td>" } @$row);
        $hover .= " </tr>";
    }
    # my $info = join("", map { chomp; "$_<br/>" } @lines);
    return $hover;
}


sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
