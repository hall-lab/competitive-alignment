#!/usr/bin/env genome-perl

use strict;
use warnings;
use List::Util qw(max);

my $n = 100;
while (my $line = <>) {
    chomp $line;
    my ($chr, $pos) = split("\t", $line);
    my $start = max(1,$pos-$n);
    my $stop = $pos+$n;
    print "$chr:$start-$stop\n";
}
