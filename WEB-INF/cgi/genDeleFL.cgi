#!/usr/bin/perl -w

#author: Xue, Li
#date: Dec. 3rd, 2015
#
# To generate deletion files for PS-HomPPI server
#
#
use strict;
use diagnostics;
use CGI ":standard";
use CGI::Carp "fatalsToBrowser";

use lib '/home/lxue/perl5/lib/perl5';
use PSHomPPI_resiPairs;
use myfun_genDelFLcgi;

#------- globle variables

our $SERVER  = 'smtp.psu.edu';    #'mailhub.iastate.edu';
our $FROM    = 'L.Xue@uu.nl';     #'lxue@ist.psu.edu';    #'lixue@iastate.edu';
our $ENVpath = $ENV{'PATH'};

our $serverName = 'PSHOMPPIv1.4';
my $serverDIR = "/data/web_servers/lxue/$serverName";
our $HOMPPIcgibin = "$serverDIR/WEB-INF/cgi";

our $dataDIR =
  "$serverDIR/uploadData/DelFL_data"
  ;                               #store all the user-uploaded data for HOMPPI
our $logDIR = "$dataDIR/LOGs";

our $safe_filename_characters = "a-zA-Z0-9_\.-";

our $mailprog =
  "/usr/sbin/sendmail";    # Location of sendmail; may vary on your system

#-------- Variables
my $email    = param('email')    || 'No input';
my $jobtitle = param('JobTitle') || 'No input';

my $qryIDpairs =
  param('QryIDpairs') || 'No input';    #query ID pairs pasted by the user
my $qrySeqs = param('QrySeqs') || 'No input';    #fasta file pasted by the user
my $seqIdentityThr = param('seqIdentityThr');
my $deletionType   = param('deletionType');

#--------making the JobTitle safe------------

my $jobID = $jobtitle;    #each submission is assigned an unique $jobID

#$jobID=~tr/ /_/;
$jobID =~ s/[^$safe_filename_characters]/_/g;

if ( $jobID =~ /^([$safe_filename_characters]+)$/ ) {
    $jobID = $1;
}
else {
    die "Jobtitle contains invalid characters.\n";
}

#--------generate unique jobID------------
my $t = localtime;
$t =~ s/[\s\:]/_/g;
my $random_number = int( rand(100) );
$jobID = $t . "." . $random_number . '.' . $jobID;

#Directories
my $result_html = $jobID . "_result.html";
my $jobDIR      = "$dataDIR/$jobID";

#chmod  755, $dataDIR;
#rmtree($jobDIR) if(-d $jobDIR);
mkdir( $jobDIR, 0755 ) || die("cannot makedir $jobDIR:$!");

#---------File variables----------------
my $protSeqFL = "$jobDIR/qrySeqs.fasta.txt";

#my $proteinFL1="$jobDIR/qryProt1.fasta.txt"; #save uploaded file or pasted seq as $proteinFL
#my $proteinFL2="$jobDIR/qryProt2.fasta.txt"; #save uploaded file or pasted seq as $proteinFL
my $qryPairLstFL = "$jobDIR/qryIDpair.lst";
my $delFL        = "$jobDIR/delete.txt";

#--------           log file       ---------#
mkdir $logDIR if ( !-d $logDIR );

our $logFL = "$logDIR/$jobID.log";

unlink $logFL if ( -e $logFL );
open LOG,    ">>$logFL" or die("Cannot write to $logFL:$!");
open STDERR, ">>$logFL" or die("Cannot write to $logFL:$!");

#-------- check the number of input ---------#

if (   $email eq 'No input'
    || $qryIDpairs eq 'No input'
    || $qrySeqs eq 'No input'
    || $seqIdentityThr eq 'No input' )
{

    print header, start_html('warning'),
'<font size="5">Warning: <font color="red">Please input all the required fields. And try again.</font>',
      p,
"User input the following: email = $email; qry ID pairs = $qryIDpairs ; $seqIdentityThr ",
      p,
      "Qry sequence fasta file: ", p,
      $qrySeqs,
      end_html();
    exit;
}

#-------- check email ---------#
use Email::Valid;

unless ( Email::Valid->address($email) ) {

    print header;
    print
      '<font size="5">Warning: You supplied an invalid email address.</font>',
      p,
      end_html();
    exit;
}

#-------- check algorithm parameter ranges ---------#
if ( $seqIdentityThr =~ /[^\d\.]+/ ) {
    print header,
      start_html('warning'),
'<font size="5">Warning: Please enter a positive real number for sequence identity threshold. And try again.</font>',
      p,
      end_html();
    exit;
}

#----------------------------------------------------------------#
#----------           Program begins                -------------#
#----------------------------------------------------------------#

#----------check the pasted sequences or the uploaded sequences------

my $FLAG = 1;    #flag the type of input file

if ( $qrySeqs ne 'No input' ) {
    &saveFastaFile( $protSeqFL, $qrySeqs );
}

if ( $qryIDpairs ne 'No input' ) {
    &saveFastaFile( $qryPairLstFL, $qryIDpairs );
}

#------check the pasted/uploaded file is in FASTA format---------
if ( $FLAG == 1 || $FLAG == 2 ) {

    #user input/uploaded a file of protein seqs

    &checkfasta($protSeqFL);
}

#----------print 'THANK YOU' html page--
#$\ = "\n";

my @params = ( $jobtitle, $jobID, $email, $deletionType, $seqIdentityThr );

&writeThankYouHtml( \@params );

#----------send an email to the user-----
&jobSubmitEmail( \@params );

#----------predict------

if ( $deletionType =~ /^strict$/i ) {

    print LOG
"deletionType = $deletionType. Call GenDeleteFL_PS_strict to generate $delFL ...\n\n";

    &GenDeleteFL_PS_strict( $qryPairLstFL, $protSeqFL, $seqIdentityThr,
        $delFL );
}
elsif ( $deletionType =~ /^loose$/i ) {

    print LOG
"deletionType = $deletionType. Call GenDeleteFL_PS_loose to generate $delFL ...\n\n";

    &GenDeleteFL_PS_loose( $qryPairLstFL, $protSeqFL, $seqIdentityThr,
        $delFL );
}

#--------- send result email to the user
&resultEmail( \@params, $logFL );

close LOG;
close STDERR;
exit 0;



