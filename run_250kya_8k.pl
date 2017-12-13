#!/usr/bin/perl -w
use strict;

# simple script to set values for model with bottleneck 250,000 years ago, with a final population size of ~8000

my $s = 0;
my $Ninit = 10000;
my $r1 = 1;
my $nr1 = 13;
my $bsize = 2;
my $r2 = 0.0;
my $nr2 = 9986;
my $sampsize = -1;
my $mu = 0.5;
my $start_sel = -1;
my $file = "250kya";
my $when_bneck = 0;
my $com = "forward $s $Ninit $r1 $nr1 $bsize $r2 $nr2 $sampsize $mu $start_sel $file $when_bneck";
my $resp = `$com`;
print $resp;
