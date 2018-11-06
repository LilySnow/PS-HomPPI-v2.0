#!/usr/bin/env perl
#===============================================================================
#
#         FILE: xue.pl
#
#        USAGE: ./xue.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (),
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 10/31/2018 06:06:13 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use globalVariables;

#my $m = '/data/web_servers/lxue/call_hhpred/MICKSCODE/runquery_abs_report_all.py ../../uploadData/test1/P:Q/P/hhpred_raw_result';
#system($m) ==0 or die ("faile:$!");

our $PYTHON;
our $callHHpredPY ;
my $raw_resultDIR = '../../uploadData/test1/A:B/A/hhpred_raw_result';
my @command =
    ( $PYTHON, $callHHpredPY, $raw_resultDIR );

    #my @command =
    #( $PYTHON, $callHHpredPY, $raw_resultDIR, " >> $logFL_tmp" );
print "CALL HHPRED COMMAND: @command\n";

system(@command) == 0 or die("ERROR: @command failed:$?");


