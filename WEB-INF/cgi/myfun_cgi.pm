use strict;
use globalVariables;

sub getHTMLinput {

    #--------------   global variables for callhhpred -----------#
    our $flag_hhpredrun;
    our $addssPL;
    our $uniprotDB_dir;
    our $pdb70HMM_dir;
    our $uniprotDB_name;
    our $pdb70DB_name;
    our $queryName;
    our $numItr_uniprotDB;
    our $numItr_PDB70;
    our $HHR1OUTPUTLOCATION;
    our $lookuptable_dir;
    our $PPIDBseqFL;
    our $evalThr_uniprot;
    our $miniCov_hhpredPDB70;
    our $maxSID_hhpredPDB70;
    our $EvalThr_hhpredPDB70;
    our $PvalThr_hhpredPDB70;
    our $minResNum_in_alignment;

    my $email    = param('email')    || 'no_input';
    my $jobtitle = param('JobTitle') || 'no_input';
#    my $querySeqFL1 = param('querySeqFL1') || 'no_input';    #fasta file uploaded by the user
#    my $querySeqFL2 = param('querySeqFL2') || 'no_input';    #fasta file uploaded by the user
    my $qryIDpairs =
      param('QryIDpairs') || 'no_input';     #query ID pairs pasted by the user
    my $qrySeqs = param('QrySeqs') || 'no_input'; #fasta file pasted by the user
    my $pdbFL_name = param('pdbFL') || 'no_input';    # query pdb files
    my $delLst = param('delFL')
      || 'no_input'
      ; #a list of PDB ID and chain IDs that the user wants to exclude from the homolog list
    my $rasaThr = 0;    #-- this parameter is obsolete
    our $atomDistThr = param('atomDistThr')
      || 8.5;           # -- PS-HomPPI only uses CA-CA to define interface.
     #Ca-Ca distance threshold used for mapping template interating residue distances to query

    #-- prediction mode parameters:
    our $safeMode_thr = param('safeMode') || 0.7;
    our $twiMode_thr1 = param('twiMode1') || 0.5;
    our $twiMode_thr2 = param('twiMode2') || 0.4;
    our $twiMode_thr3 = param('twiMode3') || 0.2;
    our $darkMode_thr = param('darkMode') || 0.1;

    our %hhpred_param = (
        'RUNNING' => $flag_hhpredrun
        ,    #if false, abort on error. if true, ignore and attempt to continue
        'ADDSSSCRIPTLOCATION' => $addssPL,    #location of addss.pl file
        'UNIPROTDATABASELOCATION' =>
          $uniprotDB_dir,                     # location of UNIPROT DB folder
        'PDB70DATABASELOCATION' => $pdb70HMM_dir
        , #/home/mwalter/project/v0.0.02/projectfiles/databases/pdb70HMMs #location of PDB70 DB folder
        'UNIPROTDATABASENAME' =>
          $uniprotDB_name,    #uniprot20_2016_02 #name of UNIPROT DB file
        'PDB70DATABASENAME' => $pdb70DB_name,   #='pdb70' #name of PDB70 DB file
        'QUERY1NAME'        => $queryName
        , # $query1.fa #for now unsafe to modify, need to also modify HHR1OUTPUTLOCATION
        'NUMBERITERATIONSUNIPROT' => param('numItr_uniprotDB')
          || $numItr_uniprotDB
        ,    #=1 #number of hhsearch iterations against the UNIPROT DB
        'NUMBERITERATIONSPDB' => param('numItr_PDB70')
          || $numItr_PDB70
        ,    #=1 #number of hhsearch iterations against the PDB70 DB

#HHR1OUTPUTLOCATION = "%s/%s.hhr" #contains QUERY1NAME and jobFolder, DO NOT MODIFY/UNCOMMENT
#S2CLOCATION = /home/mwalter/project/v0.0.02/projectfiles/databases/s2c #location of S2C DB: currently unused
        'LOOKUPDATABASELOCATION' => $lookuptable_dir
        , # location of the folder with the lookup table created by startscript.py
        'PDBSEQRESFILE' => $PPIDBseqFL
        , #=>/home/mwalter/project/v0.0.02/projectfiles/databases/PDBallFASTA/pdb_seqres.txt # location of PDB all sequences fasta file
        'HHRUNIPROTEVALUE' => param('evalThr_uniprot')
          || $evalThr_uniprot
        , #0.0001 ##both the -e -E options in the HHsearch against UNIPROT, E-value cutoff

        #HHRPDBMINSID => 5 #-qid, dont use for now
        'HHRPDBMINCOVERAGE' => param('miniCov_hhpredPDB70')
          || $miniCov_hhpredPDB70
        ,    #5 # the -cov option in the HHsearch against pdb70, min coverage
        'HHRPDBMAXID' => param('maxSID_hhpredPDB70')
          || $maxSID_hhpredPDB70
        , #100 #the -id option in the HHsearch against pdb70, max SID cutoff (for less redundancy)
        'HHRPDBEVALUE' => param('evalThr_hhpredPDB70')
          || $EvalThr_hhpredPDB70
        , #1 #both the -e -E options in the HHsearch against pdb70, E-value cutoff
        'HHRPDBPVALUE' => param('PvalThr_hhpredPDB70')
          || $PvalThr_hhpredPDB70
        ,    #=>40 #min P-value cutoff for the HHsearch against PDB70
        'MINIMUMALIGNMENTLENGTH' => param('minResNum_in_alignment')
          || $minResNum_in_alignment
        , #=>1 #min number of residues in alignment. must be >0 or results include start => stop residues
    );

    #------------------------------------------------------------#
    my $flag_QryPDBFL;
    if ( $pdbFL_name =~ /no_input/i ) {
        $flag_QryPDBFL = 'No';
    }
    else {
        $flag_QryPDBFL = 'Yes';
    }

    #------------------------------------------------------------#
    my $flag_rmSameProt = 1;
    if ( $delLst eq 'no_input' || $delLst =~ /^\s{0,}$/ ) {
        $flag_rmSameProt = 0;
    }

    #------------------------------------------------------------#
    return (
        $email,        $jobtitle,
        $qryIDpairs,   $qrySeqs,      $pdbFL_name,   $delLst,
        $rasaThr,      $atomDistThr,  $safeMode_thr, $twiMode_thr1,
        $twiMode_thr2, $twiMode_thr3, $darkMode_thr, \%hhpred_param,
        $flag_QryPDBFL, $flag_rmSameProt
    );

}

sub genJobID {

    my $jobtitle = shift @_;

    #--------making the JobTitle safe------------

    my $jobID = $jobtitle;    #each submission is assigned an unique $jobID
    $jobID = &makeJobIDsafe($jobID);

    ##--------generate unique jobID------------
    $jobID = &makeJobIDunique($jobID);

#    my $jobID =      '2017-2-8-10-28.40.test';    #each submission is assigned an unique $jobID

    return $jobID;
}

sub prepareJobDIR {

    #-- 0. save query pair ID file
    #-- 1. save the fasta files
    #-- 2. save the query pdb file, if provided
    #-- 3. save the delete file, if provided
    #-- 4. check fasta file format

    my $jobDIR     = shift @_;
    my $qryIDpairs = shift @_;
    my $qrySeqs    = shift @_;
    my $pdbFL_name = shift @_;
    my $delLst     = shift @_;

    if (!-d $jobDIR){
        mkdir $jobDIR;
    }

    my $inputDIR = "$jobDIR/input";
    mkdir($inputDIR) if ( !-d $inputDIR );
    my $protSeqFL    = "$inputDIR/qrySeqs.fasta.txt";
    my $qryPairLstFL = "$inputDIR/qryIDpair.lst";
    my $ori_pdbFL =
      "$inputDIR/pdb/$pdbFL_name";    #save uploaded pdb files as $ori_pdbFL
    my $delFL = "$inputDIR/delete.lst";
    my $seqInt_protFL =
      "$inputDIR/seq_int_qryComplexes.lst"
      ; #No need to generate. If this file does not exist, the code will automatically generate it using '?' as actual interface.


    #-- 0. save query pair ID file
    if ( $qryIDpairs ne 'no_input' ) {
        &saveFastaFile( $qryPairLstFL, $qryIDpairs );
    }

    #-- 1. save the fasta files
    if ( $qrySeqs ne 'no_input' ) {
        &saveFastaFile( $protSeqFL, $qrySeqs );
    }

    #-- 2. save the query pdb file, if provided
    if ( $pdbFL_name ne 'no_input' ) {

        #----------save the uploaded pdb file of queries------

        if ( !-d "$inputDIR/pdb" ) {
            mkdir("$inputDIR/pdb") || dir("Cannot make dir $inputDIR/pdb:$!");
        }

        &uploadZipFile( $ori_pdbFL, 'pdbFL' )
          ;    #save uploaded file as $ori_dockFL;

        my $outdir =
          &unzipFL($ori_pdbFL);    #$outdir: the files are extracted to $outdir

    }

    #-- 3. save the delete file, if provided
    if ( $delLst ne 'no_input' &&  $delLst !~ /^\s{0,}$/ ) {
        &saveFastaFile( $delFL, $delLst );
    }

    #-- 4. check fasta file format
    &checkfasta($protSeqFL);

   return  ($qryPairLstFL, $protSeqFL, $delFL, $seqInt_protFL);
}

sub autoRefreshResultPage {
    use CGI ":standard";
    use CGI::Carp qw(fatalsToBrowser warningsToBrowser);

    our $SERVERurl;    # 'http://ailab1.ist.psu.edu/PSHOMPPIv2.0/';

    my $jobID = shift @_;
    my $url   = "$SERVERurl/$jobID/result.html";
    print LOG "Refreshing result html: $url\n";
    print header( -Refresh => "1;url=$url ", -type => 'text/html' );

    print '<html>';
    print '<body>';
    print "<p></p>";
    print " <b>Status: RUNNING</b><br>";
    print " <br>                      ";
    print " </body>                   ";
    print '</html>';

}

sub makeJobIDsafe {
    our $safe_filename_characters;

    my $jobID = shift @_;

    $jobID =~ s/[^$safe_filename_characters]/_/g;

    if ( $jobID =~ /^([$safe_filename_characters]+)$/ ) {
        $jobID = $1;
    }
    else {
        die "Jobtitle contains invalid characters:$!";
    }

    return $jobID;
}

sub makeJobIDunique {

    use POSIX qw(strftime);
    my $jobID         = shift @_;
    my $random_number = int( rand(100) );
    my $datestring    = strftime "%F%H%M", localtime $^T;   #2016-03-01 11:34:46
    $datestring =~ s/-//g;
    $datestring =~ s/\s/_/g;
    $jobID = $datestring . "_" . $random_number . "_" . $jobID;

    return $jobID;
}

sub validate_input {

    use CGI ":standard";
    use CGI::Carp "fatalsToBrowser";

    (my $email,my $jobtitle,my $qryIDpairs,my $qrySeqs,my $pdbFL_name,my $delLst,my $rasaThr,my $atomDistThr,my $safeMode_thr,my $twiMode_thr1,my $twiMode_thr2,my $twiMode_thr3,my $darkMode_thr,my $hhpred_param, my $flag_QryPDBFL, my $flag_rmSameProt) = @{shift @_};


    #-------- check query ID pairs: only allow \w ---------#
    if ($qryIDpairs !~ /^[\w\n\r:]+$/){
         print header, start_html('warning'),
'<font size="5">Warning: <font color="red">Please only use numbers (0-9) or letters (a-z, A-Z) or underscore (_) for query ID.</font>',
          p, end_html();
          exit;
    }

    #-------- check the number of input ---------#
    my $flag_paste1   = 0;
    my $flag_paste2   = 0;
    my $flag_paste    = 0;
    my $flag_QrySeqFL = 0;

#    if ( $querySeqFL1 ne 'No input' && $querySeqFL2 ne 'No input' ) {
#        $flag_QrySeqFL = 1;
#    }

    if ( $pdbFL_name !~ 'no_input' && $pdbFL_name !~ /(zip|tar\.gz)$/ ) {

        #if the user uploaded the wrong format, then exit the program!
        print LOG "User upload pdb files for queries:  $pdbFL_name\n";
        print header, start_html('warning'),
'<font size="5">Warning: <font color="red">Please upload PDB formated files of your docked models in the format of *.tar.gz OR *.zip. And try again.</font>',
          p, end_html();
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

    if ( $atomDistThr =~ /[^\d\.]+/ ) {
        print header,
          start_html('warning'),
'<font size="5">Warning: Please enter a positive number for Atom distance Threshold. And try again.</font>',
          p,
          end_html();
        exit;
    }

    if ( $rasaThr =~ /[^\d\.]+/ ) {
        print header,
          start_html('warning'),
'<font size="5">Warning: Please enter a positive number for Atom distance Threshold. And try again.</font>',
          p,
          end_html();
        exit;
    }
    elsif ( $rasaThr < 0 || $rasaThr > 100 ) {
        print header,
          start_html('warning'),
'<font size="5">Warning: Please enter a number between 0 and 100 for The relative solvent accessible area (RASA) Threshold. And try again.</font>',
          p,
          end_html();
        exit;

    }
#
#if($method_XRAY ==1 && $resolutionThr !~ /^\d+[\.\d]{0,}[\s\t]{0,}$/){
#	print header,
#		start_html('warning'),
#				"<font size=\"5\">Warning: Please enter a positive number for X-RAY resolution Threshold. And try again.</font>",p,
#		end_html();
#		exit;
#}

}

#------------
sub saveFastaFile {

#author: Xue, Li
#date: Apr., 2009
#This script is part of the script for the RNA-protein interface prediction webserver. It is called by RPI_predict.cgi
#
#To save user input sequences into fasta file to /uploadData/HOMPPI/$jobID/$taskID/$taskID.fasta.txt.
#$taskID is unique for each query
#
#Usage: perl ./saveFile1.pl $ProteinToBePredictedDIR/$fastaFL  $QrySeq
#e.g.: perl ./saveFile1.pl /uploadData_RPI/test/xue.txt 'UOOIENINIONONINONIONO'

    use strict;
    use diagnostics;

    use CGI ":standard";
    use CGI::Carp qw(fatalsToBrowser);
    use File::Basename;

    #my $LOG;
    #BEGIN {
    #use CGI::Carp qw(carpout);
    #open $LOG, '>>', '/usr/local/apaches/logs/carperror.log' or die
    #"Cannot open file: $!\n";
    #carpout($LOG);
    #}

    #Variables
    my $proteinFL = shift @_;

    #Directories
    my $jobDIR = dirname($proteinFL);
    my $qrySeq = shift @_;              #fasta file pasted by the user

    #File variales

#---------------------------------------------------------------------------------#
#
#              program begins
#

    #mkdir( $jobDIR, 0777 ) || die("cannot makedir $jobDIR\n");
    unlink $proteinFL if ( -e $proteinFL );

    #--------save file---------
    open( INPUT, ">>$proteinFL" ) || fail($proteinFL);
    print INPUT "$qrySeq";
    close(INPUT);

    #---------

    sub fail {
        my $FL = shift @_;
        print "<title>Error</title>", "<p>Error: cannot open $FL:$!</p>";
        exit;
    }
}

sub uploadFile {

    #save uploaded FASTA file as $proteinFL = $ProteinToBePredictedDIR/$fastaFL
    #Author: Xue, Li
    #date: May, 2009
    #
    #This script is part of RNA-protein interface predictioin
    #
    #This script is to upload FASTA query file into $ProteinToBePredictedDIR
    #The maximum upload filesize is 100K

    use strict;

    #use CGI ":standard";
    my $masFLSize = 100;    # maximum upload filesize is 100K
    $CGI::POST_MAX = 1024 * $masFLSize;    # maximum upload filesize is 100K

    my $query          = new CGI;
    my $proteinFL      = shift @_;
    my $uploadFLhandle = shift @_;
    my $jobDIR         = dirname($proteinFL);

    #----- Look for uploads that exceed $CGI::POST_MAX-------

    if ( !$query->param('querySeqFL') && $query->cgi_error() ) {
        print $query->cgi_error();
        print p,
"The file you are attempting to upload exceeds the maximum allowable file size $masFLSize K.",
          p, 'Please refer to your system administrator', p;
        print $query->hr, $query->end_html;
        exit 0;
    }

    #----- Upload file-------

    my $upload_filehandle = $query->upload($uploadFLhandle);

    #mkdir( $jobDIR, 0777 ) || die("Cannot make dir\n");
    unlink $proteinFL if ( -e $proteinFL );

    open UPLOADFILE, ">$proteinFL" || die("Cannot open $proteinFL!");

    while (<$upload_filehandle>) {
        s/[\n\r]//mg;
        print UPLOADFILE "$_\n";
    }

    close UPLOADFILE;

}

#-------------

sub checkfasta {

    #Xue, Li
    #May 2009
    #check input file format is FASTA format or not

    my $FL     = shift @_;
    my $lineNo = 0;
    open( INPUT, "<$FL" ) || die("Cannot open $FL!\n");
    foreach my $line (<INPUT>) {
        $line =~ s/[\n\r]//mg;    #remove \r and \n
        $lineNo = $lineNo + 1;

        if ( $line =~ /^#/ ) {
            next;
        }
        elsif ($line =~ /^\>.+$/
            || $line =~ /^[a-zA-Z]+$/
            || $line =~ /^[\s\t]{0,}$/ )
        {

            next;
        }
        else    #this line is not in FASTA format
        {
            print header, start_html('warning'),
"<font size=\"5\">Warning: Please check the input file format.Start from line $lineNo.</font>",
              p, "<a href=\"SEQexample.txt\"> Click here for an example<\/a>",
              p, end_html();
            exit 1;
        }

    }
}

#------------

sub checkPDBChainIDFL {

    my $FL     = shift @_;
    my $lineNo = 0;

    open( INPUT, "<$FL" ) || die("Cannot open $FL!\n");
    foreach (<INPUT>) {
        $lineNo++;

        s/[\n\r]//mg;

        if (/^\#/) {

            #this is a comment line
            next;
        }

        elsif (/^[\da-zA-Z]{4}[\s\t]{0,},[\s\t]{0,}[a-zA-Z\d]{1}$/
            || /^[\s\t]{0,}$/ )
        {

            #this is a normal pdb/chain ID line or an empty line
            next;
        }

        else {

            #this line is not in the format of HOMPPI server
            print header, start_html('Check the input file format'), p,
              "Please check the input file format. Start from line $lineNo.", p,
"<a href=\"PDBIDChainIDexample.txt\"> Click here for an example<\/a>",
              p, end_html();
            exit 1;
        }

    }
    close(INPUT);

}

#------------
sub convertFasta2onelineFormat {

    #perl convertFasta2onelineFormat.pl ../data/dockingTestSet.fasta
    use strict;

    my $inputFL = shift @_;
    our $dataDIR;
    my $tempFL = "$dataDIR/xue.tmp";
    my $seq;

    unlink $tempFL if ( -e $tempFL );
    open( OUTPUT, ">>$tempFL" ) || die("Cannot open $tempFL\n");

    open( INPUT, "<$inputFL" ) || die("Cannot open $inputFL\n");

    foreach (<INPUT>) {
        s/[\n\r]+//mg;

        if (/>/) {
            if ($seq) {
                print OUTPUT "\n";
            }
            print OUTPUT "$_\n";
            next;
        }

        if (/^([a-zA-Z]+)$/) {
            $seq = $1;
            print OUTPUT $seq;

        }

    }
    close(INPUT);
    close(OUTPUT);

    unlink($inputFL);
    rename( $tempFL, $inputFL );

}

sub writeThankYouHtml {

    #insert input parameters to $basicISUhtml
    use strict;
    use CGI ":standard";

    our $serverName;

    my (
        $jobtitle, $jobID, $email,
        $intDef, $atomDistThr, $k, $flag_rmSameProt, $flag_QryPDBFL, $pdbFL_name

    ) = @{ shift @_ };

    #	my $method = &getMethod($methodCode);

    my $basicISUhtml = 'basicISU.html';

    open( INPUT, " < $basicISUhtml " )
      || die(" Cannot open $basicISUhtml:$!");

    print header;

    while (<INPUT>) {
        if (/<!-- INSERT your content here -->/) {
            print "$_ ";
            print "Thank you for using $serverName.</h3>";
            print " Your job <i>$jobtitle </i> is currently running ......";
            print '<p></p> ';
            print
" The zipped results will be downlaodable from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID.tar.gz\">here</a>.<p></p>";
            print
' You will be notified by email after the prediction is finished. ';
            print ' <p> </p>';
            print "Please check <font color=\"red\">$email </font>is correct.";
            print ' <p> </p>';
            print ' <h3> Query Parameters </h3>';
            print ' <TABLE border = 0> ';
            print
"<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  &#197</TD></TR>";

        #print "TR> <TD><FONT COLOR=\"Blue\"> <hr /></FONT></TD><TD></TD></TR>";
            print
"<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>";
            print
"<TR> <TD><FONT COLOR=\"Blue\"> User provided query PDB files: </FONT></TD><TD> $flag_QryPDBFL</TD></TR>";

#			print
#"<TR> <TD><FONT COLOR=\"Blue\"> Types of Homologs to be used: </FONT></TD><TD> $method </TD></TR>";
#
#			if ( $method =~ /XRAY/ ) {
#				print
#"<TR> <TD><FONT COLOR=\"Blue\"> X-RAY resolution threshold </FONT></TD><TD> $resolutionThr </TD></TR>";
#			}

            print ' </TABLE>';
        }
        else {
            print "$_";
        }
    }
    close(INPUT);

}

#------------

sub getMethod {
    my $methodCode = shift @_;
    my $method;

    if ( $methodCode == 1 ) {
        $method = 'XRAY';

    }
    elsif ( $methodCode == 2 ) {
        $method = 'NMR';
    }
    elsif ( $methodCode == 3 ) {
        $method = 'XRAY and NMR';
    }
    else {
        die(
"Please provide valide method code: 1 for XRAY, 2 for NMR, 3 for XRAY and NMR:$!"
        );
    }
}

#------------

sub jobSubmitEmail_todo {

    use strict;
    use warnings;

    # use MIME::Lite::TT::HTML;

    our $serverName;
    our $version;
    our $SERVER;
    our $FROM;

    my (
        $jobtitle, $jobID, $email,
        $intDef, $atomDistThr,  $flag_rmSameProt, $flag_QryPDBFL, $pdbFL_name

    ) = @{ shift @_ };

    #	use Mail::Mailer;
    use Email::Valid;
    use MIME::Lite;

    #	my $method = &getMethod($methodCode);

    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.";
        exit;
    }

#	my $EmailContent = "<h3>Thank you for using PS-HomPPI.</h3>
#Your job <i>$jobtitle</i> is currently running ......<p></p>
#You will be notified by email after the prediction is finished.<p></p>
#<h3>Query Parameters</h3>
#<TABLE border=0>
#<TR> <TD><FONT COLOR=\"Blue\"> The relative solvent accessible area (RASA) Threshold</FONT></TD><TD> $rasaThr  \%</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> Atom distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> <hr /></FONT></TD><TD></TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> BLAST Positive Score Threshold for homologs</FONT></TD><TD> $PositiveS \%</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> Expected value threshold for BLAST</FONT></TD><TD> $Evalue</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Threshold of FracQH </FONT></TD><TD> $FracQH_Thr</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Threshold of FracPHP </FONT></TD><TD> $FracPHP_Thr</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Types of Homologs to be used: </FONT></TD><TD> $method</TD></TR>";

    my $EmailContent;
    if ( $flag_QryPDBFL eq 'Yes' ) {
        $EmailContent = "<h3>Thank you for using PS-HomPPI $version.</h3>
        Your job <i>$jobtitle</i> is currently running ......<p></p>
        You will be notified by email after the prediction is finished.<p></p>

        <p></p>
        <p></p>
        The fasta files of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        The PDB files of the query proteins are <a  href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/pdb/$pdbFL_name\">here</a>.<p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    }
    else {
        $EmailContent = "<h3>Thank you for using PS-HomPPI $version.</h3>
        Your job <i>$jobtitle</i> is currently running ......<p></p>
        You will be notified by email after the prediction is finished.<p></p>

        <p></p>
        <p></p>
        The fasta files of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>";
    }

    my $delFL_line =
"<p></p>The delete file is <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/delete.lst\">here</a>.<p></p>";

    if ( $flag_rmSameProt == 1 ) {
        $EmailContent = $EmailContent . $delFL_line . '<p></p>';
    }
    my $paramTable = "<h3>Query Parameters</h3>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> User provided query PDB files: </FONT></TD><TD> $flag_QryPDBFL</TD></TR>
";

    $EmailContent = $EmailContent . $paramTable;

#	if ( $method =~ /XRAY/ ) {
#		$EmailContent =
#		  $EmailContent
#		  . "<TR> <TD><FONT COLOR=\"Blue\"> X-RAY resolution threshold </FONT></TD><TD> $resolutionThr </TD></TR>";
#	}
#
    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . " <p></p><p></p> If you have questions regarding <a href = \"http://ailab1.ist.psu.edu/$serverName/\">$serverName</a>, please visit the corresponding web page(s) or write to Li Xue "
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName Submission  Details",
        From    => $FROM,
        To      => $email,
        Type    => 'text/html',
        Data    => $EmailContent
    );

    $msg->send( 'smtp', $SERVER, Timeout => 120 );

    print LOG "A submission confirmation email is sent to $email. \n";

}

sub resultEmail {

    use strict;
    use warnings;

    #	use Mail::Mailer;
    use Email::Valid;
    use MIME::Lite;

    # use MIME::Lite::TT::HTML;
    our $serverName;
    our $SERVER;    #email sever
    our $FROM;
    our $dataDIR;
    my $params_ref = shift @_;

    #	my $jobID      = shift @_;
    my $logFL = shift @_;

    #my $jobDIR     = "$dataDIR/$jobID";
    my (
        $jobtitle,        $jobID,         $email,
        $intDef,          $atomDistThr,   $k,
        $flag_rmSameProt, $flag_QryPDBFL, $pdbFL_name
    ) = @$params_ref;

    #print header,
    #      start_html('warning'),
    #"flag_qryPDB=$flag_QryPDBFL, pdbFL_name = $pdbFL_name",
    #      p,
    #      end_html();
    #      exit;

    #get the type of homolog: X-RAY or NMR
    #	my $method = &getMethod($methodCode);

    #

    #check the validation of email address
    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.\n";
        exit;
    }

    my $statFLname;
    if ( $flag_rmSameProt == 1 ) {
        $statFLname = 'statistics_qryIDpair_wo_sameProt.txt';
    }
    else {
        $statFLname = 'statistics_qryIDpair.txt';
    }

    my $EmailContent;
    if ( $flag_QryPDBFL eq 'Yes' ) {
        $EmailContent = "<h3>Thank you for using $serverName.</h3>
        The prediction results for your job <i>$jobtitle</i> are now available.<p></p>

        <h4>Results</h4>

        Predicted CA-CA distances for the interface residue pairs can be downloaded from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/Ca_Ca_distances.tar.gz\">here</a>.<p></p>

        The statistics of sequence similarity of query and its homo-interologs can be downloaded <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/$statFLname\">here</a>. Please note that not all the homo-interologs listed in this file was used by PS-HomPPI. <p></p>

        The complete run can be downloaded as a zipped tar file from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID.tar.gz\">here</a>.<p></p>



        <h4>User input Data</h4>

        <p></p>
        The fasta file of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        <p></p>
        The PDB files of the query proteins are <a  href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/pdb/$pdbFL_name\">here</a>.<p></p>
        <p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    }
    else {
        $EmailContent = "<h3>Thank you for using $serverName.</h3>
        The prediction results for your job <i>$jobtitle</i> are now available.<p></p>

        <h4>Results</h4>

        Predicted CA-CA distances for the interface residue pairs can be downloaded from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/Ca_Ca_distances.tar.gz\">here</a>.<p></p>

        The statistics of sequence similarity of query and its homo-interologs can be downloaded <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/$statFLname\">here</a>. Please note that not all the homo-interologs listed in this file was used by PS-HomPPI. <p></p>

        The complete run can be downloaded as a zipped tar file from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID.tar.gz\">here</a>.<p></p>



        <h4>User input Data</h4>

        <p></p>
        The fasta file of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        <p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    }

    my $delFL_line =
"<p></p>The delete file is <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/delete.lst\">here</a>.<p></p>";

    if ( $flag_rmSameProt == 1 ) {

        $EmailContent = $EmailContent . $delFL_line . '<p></p><p></p>';
    }
    my $paramTable = "<h4>Query Parameters</h4>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> User provided query PDB files: </FONT></TD><TD> $flag_QryPDBFL</TD></TR>";

    $EmailContent = $EmailContent . $paramTable;

#	if ( $method =~ /XRAY/ ) {
#		$EmailContent =
#		  $EmailContent
#		  . "<TR> <TD><FONT COLOR=\"Blue\"> X-RAY resolution threshold </FONT></TD><TD> $resolutionThr </TD></TR>";
#	}

    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . ' <p></p><p></p> If you have questions regarding <a href = "http://ailab1.ist.psu.edu/'
      . $serverName . '/">'
      . $serverName
      . '</a>, please visit the corresponding web page(s) or write to Li Xue '
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName Prediction Results",
        From    => $FROM,
        To      => $email,
        Type    => 'text/html',
        Data    => $EmailContent
    );

    $msg->send( 'smtp', $SERVER, Timeout => 60 );

    #	$msg->send();

    open( LOG, ">>$logFL" ) || die("Cannot open $logFL.\n");
    print LOG "\n\nA prediction result email is sent to $email. \n";
    close LOG;

}

#------------

sub resultEmail_old {

    use strict;
    use warnings;

    #	use Mail::Mailer;
    use Email::Valid;
    use MIME::Lite;

    # use MIME::Lite::TT::HTML;
    our $serverName;
    our $SERVER;    #email sever
    our $FROM;
    our $dataDIR;
    my $params_ref = shift @_;

    #	my $jobID      = shift @_;
    my $logFL = shift @_;

    #my $jobDIR     = "$dataDIR/$jobID";
    my ( $jobtitle, $jobID, $email,
        $intDef, $atomDistThr, $k, $flag_rmSameProt ) = @$params_ref;

    #get the type of homolog: X-RAY or NMR
    #	my $method = &getMethod($methodCode);

    #

    #check the validation of email address
    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.\n";
        exit;
    }

#	my $EmailContent = "<h3>Thank you for using $serverName.</h3>
#The prediction results for your job <i>$jobtitle</i> are now available.<p></p>
#You may download the results <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/predictionsFor$jobID.txt\">here</a>.<p></p>
#
#<p></p>
#
#	<p></p>
#You may also download the list of homo-interologs <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/homo-interologs.txt\">here</a>.<p></p>
#
#
#<p></p>
#The fasta file of the query protein 1 is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryProt1.fasta.txt\">here</a>.<p></p>
#
#<p></p>
#The fasta file of the query protein 2 is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryProt2.fasta.txt\">here</a>.<p></p>
#
#<p></p>
#The BLASTp result file for query protein 1 is <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryProt1_qryProt2/qryProt1/qryProt1.blast\">here</a>.<p></p>
#The BLASTp result file for query protein 2 is <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryProt1_qryProt2/qryProt2/qryProt2.blast\">here</a>.<p></p>
#<p></p>
#
#<h3>Query Parameters</h3>
#<TABLE border=0>
#<TR> <TD><FONT COLOR=\"Blue\"> The relative solvent accessible area (RASA) Threshold</FONT></TD><TD> $rasaThr  \%</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> <hr /></FONT></TD><TD></TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> BLAST Positive Score Threshold for homologs</FONT></TD><TD> $PositiveS \%</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> Expected value threshold for BLAST</FONT></TD><TD> $Evalue</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Threshold of FracQH </FONT></TD><TD> $FracQH_Thr</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Threshold of FracPHP </FONT></TD><TD> $FracPHP_Thr</TD></TR>
#<TR> <TD><FONT COLOR=\"Blue\">  Types of Homologs to be used: </FONT></TD><TD> $method</TD></TR>";
#

    my $statFLname;
    if ( $flag_rmSameProt == 1 ) {
        $statFLname = 'statistics_qryIDpair_wo_sameProt.txt';
    }
    else {
        $statFLname = 'statistics_qryIDpair.txt';
    }

    my $EmailContent = "<h3>Thank you for using $serverName.</h3>
The prediction results for your job <i>$jobtitle</i> are now available.<p></p>

<h4>Results</h4>

You may download the interface prediction result file
<a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/InterfacePredictions.final.txt\">here</a>.<p></p>

A machine-friendly format of the interface prediction result file <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/InterfacePredictions.internal.txt\">here</a>.<p></p>

Predicted CA-CA distances for the interface residue pairs can be downloaded from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/resiPairs/predicted_resiPair_distance.tar.gz\">here</a>.<p></p>

The statistics of sequence similarity of query and its homo-interologs can be downloaded <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/$statFLname\">here</a>. Please note that not all the homo-interologs listed in this file was used by $serverName.

<h4>User input Data</h4>

<p></p>
The fasta files of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>

<p></p>
The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    my $delFL_line =
"<p></p>The delete file is <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/delete.lst\">here</a>.<p></p>";

    if ( $flag_rmSameProt == 1 ) {

        $EmailContent = $EmailContent . $delFL_line . '<p></p><p></p>';
    }
    my $paramTable = "<h4>Query Parameters</h4>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>";

    $EmailContent = $EmailContent . $paramTable;

#	if ( $method =~ /XRAY/ ) {
#		$EmailContent =
#		  $EmailContent
#		  . "<TR> <TD><FONT COLOR=\"Blue\"> X-RAY resolution threshold </FONT></TD><TD> $resolutionThr </TD></TR>";
#	}

    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . ' <p></p><p></p> If you have questions regarding <a href = "http://ailab1.ist.psu.edu/'
      . $serverName . '/">'
      . $serverName
      . '</a>, please visit the corresponding web page(s) or write to Li Xue '
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName Prediction Results",
        From    => $FROM,
        To      => $email,
        Type    => 'text/html',
        Data    => $EmailContent
    );

    $msg->send( 'smtp', $SERVER, Timeout => 60 );

    #	$msg->send();

    open( LOG, ">>$logFL" ) || die("Cannot open $logFL.\n");
    print LOG "\n\nA prediction result email is sent to $email. \n";
    close LOG;

}

#------------

sub resultFailEmail {

    use strict;
    use warnings;

    #	use Mail::Mailer;
    use Email::Valid;
    use MIME::Lite;

    # use MIME::Lite::TT::HTML;
    our $serverName;
    our $SERVER;
    our $FROM;
    our $dataDIR;
    my $params_ref = shift @_;

    #	my $jobID      = shift @_;
    my $logFL = shift @_;

    #my $jobDIR     = "$dataDIR/$jobID";

    my (
        $jobtitle,        $jobID,         $email,
        $intDef,          $atomDistThr,   $k,
        $flag_rmSameProt, $flag_QryPDBFL, $pdbFL_name
    ) = @$params_ref;

    #	my $method = &getMethod($methodCode);

    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.\n";
        exit;
    }

    my $statFLname;
    if ( $flag_rmSameProt == 1 ) {
        $statFLname = 'statistics_qryIDpair_wo_sameProt.txt';
    }
    else {
        $statFLname = 'statistics_qryIDpair.txt';
    }

    my $EmailContent = "<h3>Thank you for using $serverName.</h3>
The prediction results for your job <i>$jobtitle</i> are now available.<p></p>


Unfortunately, no homo-interologs are found by PS-HomPPI using the thresholds that user specified.<p></p>
User may loosen the thresholds to include remote sequence homologs into prediction. For parameters of PS-HomPPI of Safe/Twilight/Dark zones, please refer to <p></p>

<ul>Li C. Xue, Drena Dobbs, Vasant Honavar: HomPPI: A Class of Sequence Homology Based Protein-Protein Interface Prediction Methods.
BMC Bioinformatics 2011, 12:244.</ul>


<p></p>

<p></p>
The fasta files of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
The PDB files of the query proteins are <a  href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/pdb/$pdbFL_name\">here</a>.<p></p>
<p></p>
The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/qryIDpair.lst\">here</a>.<p></p>

	The statistics of sequence similarity of query and its homo-interologs can be downloaded <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/$jobID/$statFLname\">here</a>. Please note that not all the homo-interologs listed in this file was used by $serverName.
<p></p>
<p></p>




<p></p>
<h3>Query Parameters</h3>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> CA-CA distance Threshold</FONT></TD><TD> $atomDistThr  Å</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> K in K-nearest homologs</FONT></TD><TD> $k</TD></TR>";

#<TR> <TD><FONT COLOR=\"Blue\">  Types of Homologs to be used: </FONT></TD><TD> $method</TD></TR>";

#	if ( $method =~ /XRAY/ ) {
#		$EmailContent =
#		  $EmailContent
#		  . "<TR> <TD><FONT COLOR=\"Blue\"> X-RAY resolution threshold </FONT></TD><TD> $resolutionThr </TD></TR>";
#	}
#
    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . ' <p></p><p></p> If you have questions regarding <a href = "http://ailab1.ist.psu.edu/$serverName/">'
      . $serverName
      . '</a>, please visit the corresponding web page(s) or write to Li Xue '
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName Prediction Results - No Homo-interolog found",
        From    => $FROM,
        To      => $email,
        Type    => 'text/html',
        Data    => $EmailContent
    );

    $msg->send( 'smtp', $SERVER, Timeout => 60 );

    #	$msg->send();

    open( LOG, ">>$logFL" ) || die("Cannot open $logFL:$!");
    print LOG "\n\nA prediction result email is sent to $email. \n";
    close LOG;

}

sub sendResultEmail {

    my $flag_prediction = shift @_;
    my @params          = @{ shift @_ };
    my $logFL           = shift @_;

    if ( $flag_prediction == 0 ) {

#         $flag_prediction=0: at least one of the protein pairs have homo-interologs in $FullStatFL
#send prediction result email

        #&resultEmail( \@params, $logFL );

    }

    if ( $flag_prediction == 1 ) {

        #no prediction can be made. Send prediction result email.
        #&resultFailEmail( \@params, $logFL );
    }
    return;
}

#------------
1;
