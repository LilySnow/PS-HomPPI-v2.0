#!/usr/bin/env perl
#===============================================================================
#
#  DESCRIPTION: generate haddock restraint file based on CA-CA file
#       INPUT (ca-ca file. The first 4 columns and the last 3 columns are used):
#           chnID_Qry1  aa_Qry1 chnID_Qry2  aa_Qry2 A_template1 aa1_template1   B_template1 aa2_template1   mean    min max
#           A   65  B   73  A 1806 B 76 15.55 14.80 16.20
#           A   65  B   46  A 1806 B 49 13.92 13.52 14.52
#
#       OUTPUT:
#           assign (name ca and segid A and resi 1806) (name ca and segid B and resi 76) 15.550 1.250 1.150
#           assign (name ca and segid A and resi 1806) (name ca and segid B and resi 49) 13.920 0.900 1.100
#
#      CREATED: 10/12/2015 08:46:24 AM
#===============================================================================

use strict;
use warnings;
use utf8;

my $cacaFL    = shift @ARGV;
my $outputFL  = shift @ARGV;
my $CaDistThr = shift @ARGV;

&header( $outputFL, $CaDistThr );

open( INPUT,  "<$cacaFL" )    or die("Cannot open $cacaFL:$!");
open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");

while (<INPUT>) {
    s/[\n\n]//gm;
    if (/^\w{1}\s+[-\w]+\s+\w{1}\s+/) {

        #           A   65  B   73  A 1806 B 76 15.55 14.80 16.20
        my @a            = split( /\s+/, $_ );
        my $chnID_qry1   = $a[0];
        my $resiNum_qry1 = $a[1];
        my $chnID_qry2   = $a[2];
        my $resiNum_qry2 = $a[3];
        my $max          = pop @a;
        my $min          = pop @a;
        my $mean         = pop @a;

        #--give 0.5 angstroms more flexibility to the restraints

        my $lower = $mean - $min + 0.5;
        my $upper = $max - $mean + 0.5;

        printf OUTPUT
"assign (name ca and segid $chnID_qry1 and resi %s) (name ca and segid $chnID_qry2 and resi %s) %.3f %.3f %.3f\n",
          $resiNum_qry1, $resiNum_qry2, $mean, $lower, $upper;

    }

}
close INPUT;
close OUTPUT;

print "$outputFL generated\n";

#-------------------------

sub header {
    my $outputFL  = shift @_;
    my $CaDistThr = shift @_;

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    print OUTPUT "! generated by $0\n";
    print OUTPUT "! CA-CA cutoff = $CaDistThr\n";
    print OUTPUT
"! 0.5 angstroms are added both sides of the distance restraints derived from templates\n";
    close OUTPUT;

}

