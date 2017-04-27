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

#--------------  Input Variables -------------#

my @params = &getHTMLinput;

(my $email,my $jobtitle,my $jobID, my $qryIDpairs,my $qrySeqs,my $pdbFL_name,my $delLst,my $rasaThr,our $atomDistThr,our $safeMode_thr,our $twiMode_thr1,our $twiMode_thr2,our $twiMode_thr3,our $darkMode_thr,our $hhpred_param, my $flag_QryPDBFL, my $flag_rmSameProt) = @params;


#--------           XRAY or NMR       ---------#
#my $methodCode;
#my $method_XRAY=param('method_XRAY')||0;
#my $method_NMR=param('method_NMR')||0;
#$methodCode=$method_XRAY + $method_NMR;#1: XRAY, 2:NMR, 3. XRAY + NMR
#my $resolutionThr = param('resolutionThr') || 3.5;

#--------           validate input variables        ---------#
&validate_input (\@params);


#--------           log file       ---------#
mkdir $logDIR if ( !-d $logDIR );
our $logFL = "$logDIR/$jobID.log";
unlink $logFL if ( -e $logFL );
open LOG,    ">>$logFL" or die("Cannot write to $logFL:$!");
open STDERR, ">>$logFL" or die("Cannot write to $logFL:$!");

#my $submitter_IP = $ENV{REMOTE_ADDR};
#&recordSubmitInfo( $submitter_IP, $email );


#----------------------------------------------------------------#
#----------           Program begins                -------------#
#----------------------------------------------------------------#

#-------- 1. prepare jobDIR       ---------#
my ($qryPairLstFL, $protSeqFL, $delFL, $seqInt_protFL) = &prepareJobDIR($jobDIR,      $qryIDpairs  ,     $qrySeqs  ,     $pdbFL_name  ,     $delLst  );

#-- step 1a:  Print auto-refresh result page
##
#print LOG "\n*step 3: Print auto-refresh result page.\n";
#
#my $indexHtml_template = '../html/result.indexHtml';
#copy($indexHtml_template, "$jobDIR/result.html");
#&autoRefreshResultPage($jobID);
#&writeThankYouHtml( \@params );

#--------- 2. send an email to the user-----#
#&jobSubmitEmail( \@params );


#----------3. predict ------#
#my $flag_prediction =
#  &batch_PSHomPPI( $jobDIR, $qryPairLstFL, $protSeqFL, $delFL, $seqInt_protFL,
#    $intDef, $atomDistThr, $rasaThr );


my $flag_prediction = 0;

#----------4. send result email ------#
#-- $flag_prediction=1: none of the query pairs have homo-interologs in $FullStatFL
if ( !defined $flag_prediction ) {
    die("flag_prediction not defined:$!");
}

&sendResultEmail($flag_prediction, \@params, $logFL);

#-- clean up
#chmod_recursive( 'g+w', $jobDIR);
chmod_recursive( 0755, $jobDIR);
chmod 0755, "$jobDIR.tar.gz";

close LOG;
close STDERR;


exit 0;

