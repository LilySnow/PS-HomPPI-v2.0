#!/usr/bin/env perl -w
#===============================================================================
#
#        Li Xue, L.Xue@uu.nl
#        11/30/2016 05:50:43 PM
#
#  DESCRIPTION: cgi for PS-HomPPI v2.0
#
#===============================================================================

use strict;
use File::Path;
use diagnostics;
use CGI ":standard";
use CGI::Carp "fatalsToBrowser";

use PSHomPPI_resiPairs;
use myfun_cgi;
use globalVariables;
our $perl5LibDIR;
use lib $perl5LibDIR;
use File::chmod::Recursive;

#--------------   set environment variables -----------#
our $PSHomPPI_path;
$ENV{'PATH'} = "$PSHomPPI_path";

#--------------   global variables -----------#
our $safe_filename_characters;    #  "a-zA-Z0-9_\.-";
our $dataDIR
  ; # "$serverDIR/uploadData";        #store all the user-uploaded data for HOMPPI
our $logDIR;
our $intDef;

#----------------------------------------------------------------#
#----------           Program begins                -------------#
#----------------------------------------------------------------#

my $jobID = shift @ARGV;    # '201704201451-71-T122'
mkdir $logDIR if ( !-d $logDIR );
our $logFL = "../../uploadData/LOGs/$jobID.log";
unlink $logFL if ( -e $logFL );
open LOG,    ">>$logFL" or die("Cannot write to $logFL:$!");
open STDERR, ">>$logFL" or die("Cannot write to $logFL:$!");

#----------predict------
my $jobDIR       = "../../uploadData/$jobID";
my $qryPairLstFL = "$jobDIR/input/qryIDpair.lst";
my $protSeqFL    = "$jobDIR/input/qrySeqs.fasta.txt";
my $delFL        = "$jobDIR/input/delete.lst";
my $seqInt_protFL =
  "$jobDIR/input/seq_int_qryComplexes.lst"
  ; #No need to generate. If this file does not exist, the code will automatically generate it using '?' as actual interface.

our $CaDistThr = 8.5;
our $atomDistThr = 8.5;
my $rasaThr     = 0;     #obsolete
our $safeMode_thr = 0.7;
our $twiMode_thr1 = 0.5;
our $twiMode_thr2 = 0.4;
our $twiMode_thr3 = 0.2;
our $darkMode_thr = 0.1;

#----------3. predict ------#
my $flag_prediction =
  &batch_PSHomPPI( $jobDIR, $qryPairLstFL, $protSeqFL, $delFL, $seqInt_protFL,
    $intDef, $atomDistThr, $rasaThr );

#----------4. send result email ------#
#-- $flag_prediction=1: none of the query pairs have homo-interologs in $FullStatFL
if ( !defined $flag_prediction ) {
    die("flag_prediction not defined:$!");
}

#-- clean up
#chmod_recursive( 'g+w', $jobDIR);
chmod_recursive( 0755, $jobDIR );
chmod 0755, "$jobDIR.tar.gz";

close LOG;
close STDERR;

exit 0;

