#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    chomp;
    if (/>/) {
        s/>//;
        print "$_\t"; # just the seq id
    }
    elsif (/\((-\d+\.\d+)\)$/) {
        my $mfe = $1;
        s/ \($mfe\)//;
        print "$_\t$mfe\n"; # fold + MFE
    }
    else {
        print "$_\t"; # the seq
    }
}
