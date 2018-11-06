#
#===============================================================================
#
#         FILE: myfun_genDelFLcgi.pm
#
#  DESCRIPTION:
#
#      CREATED: 12/03/2015 05:08:32 AM
#===============================================================================

use strict;
use warnings;

sub writeThankYouHtml {

    #insert input parameters to $basicISUhtml
    use strict;
    use CGI ":standard";

    our $serverName;
    my ( $jobtitle, $jobID, $email, $deletionType, $seqIdentityThr ) =
      @{ shift @_ };

    my $basicISUhtml = 'basicISU_genDelFL.html';

    open( INPUT, " < $basicISUhtml " )
      || die(" Cannot open $basicISUhtml:$!");

    print header;

    while (<INPUT>) {
        if (/<!-- INSERT your content here -->/) {
            print "$_ ";
            print
"Thank you for using $serverName to generate deletion files.</h3>";
            print " Your job <i>$jobtitle </i> is currently running ......";
            print '<p></p> ';
            print
" The result will be downlaodable from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/delete.txt\">here</a>.<p></p>";
            print
' You will be notified by email after the prediction is finished. ';
            print ' <p> </p>';
            print "Please check <font color=\"red\">$email </font>is correct.";
            print ' <p> </p>';
            print ' <h3> Query Parameters </h3>';
            print ' <TABLE border = 0> ';
            print
"<TR> <TD><FONT COLOR=\"Blue\"> Sequence identity threshold </FONT></TD><TD> $seqIdentityThr  %</TD></TR>";
            print
"<TR> <TD><FONT COLOR=\"Blue\"> Deletion type </FONT></TD><TD> $deletionType </TD></TR>";

            print ' </TABLE>';
        }
        else {
            print "$_";
        }
    }
    close(INPUT);

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
    open( INPUT, ">>$proteinFL" ) || fail();
    print INPUT "$qrySeq";
    close(INPUT);

    #---------

    sub fail {
        print "<title>Error</title>", "<p>Error: cannot open proteinFL!</p>";
        exit;
    }
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

sub resultEmail {

    use strict;
    use warnings;

    use Email::Valid;
    use MIME::Lite;

    our $serverName;
    our $SERVER;    #email sever
    our $FROM;
    our $dataDIR;

    my ( $jobtitle, $jobID, $email, $deletionType, $seqIdentityThr ) =
      @{ shift @_ };

    my $logFL = shift @_;

    #check the validation of email address
    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.\n";
        exit;
    }

    my $EmailContent =
      "<h3>Thank you for using $serverName to generate deletion files.</h3>
        The  results for your job <i>$jobtitle</i> are now available.<p></p>

        <h4>Results</h4>

        The deletion file can be downloaded from  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/delete.txt\">here</a>.<p></p>

        <h4>User input Data</h4>

        <p></p>
        The fasta file of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        <p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    my $paramTable = "<h4>Query Parameters</h4>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> Sequence identity threshold </FONT></TD><TD> $seqIdentityThr %</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> Deletion type </FONT></TD><TD> $deletionType</TD></TR>";

    $EmailContent = $EmailContent . $paramTable;
    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . ' <p></p><p></p> If you have questions regarding <a href = "http://ailab1.ist.psu.edu/'
      . $serverName . '/">'
      . $serverName
      . '</a>, please visit the corresponding web page(s) or write to Li Xue '
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName: Deletion File Generated",
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

sub jobSubmitEmail {

    use strict;
    use warnings;

    # use MIME::Lite::TT::HTML;

    our $serverName;
    our $SERVER;
    our $FROM;

    my ( $jobtitle, $jobID, $email, $deletionType, $seqIdentityThr ) =
      @{ shift @_ };

    use Email::Valid;
    use MIME::Lite;

    unless ( Email::Valid->address($email) ) {
        print header;
        print "You supplied an invalid email address.";
        exit;
    }

    my $EmailContent =
      "<h3>Thank you for using $serverName to generate deletion files.</h3>
        The  results for your job <i>$jobtitle</i> are now available.<p></p>

        <h4>Results</h4>

         The result will be downloadable from <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/delete.txt\">here</a>.<p></p>

        <h4>User input Data</h4>

        <p></p>
        The fasta file of the query proteins is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/qrySeqs.fasta.txt\">here</a>.<p></p>
        <p></p>
        The query ID pair file is  <a href=\"http://ailab1.ist.psu.edu/$serverName/uploadData/DelFL_data/$jobID/qryIDpair.lst\">here</a>.<p></p>";

    my $paramTable = "<h4>Query Parameters</h4>
<TABLE border=0>
<TR> <TD><FONT COLOR=\"Blue\"> Sequence identity threshold </FONT></TD><TD> $seqIdentityThr %</TD></TR>
<TR> <TD><FONT COLOR=\"Blue\"> Deletion type </FONT></TD><TD> $deletionType</TD></TR>";

    $EmailContent = $EmailContent . $paramTable;
    $EmailContent = $EmailContent . '</TABLE>';

    $EmailContent =
        $EmailContent
      . ' <p></p><p></p> If you have questions regarding <a href = "http://ailab1.ist.psu.edu/'
      . $serverName . '/">'
      . $serverName
      . '</a>, please visit the corresponding web page(s) or write to Li Xue '
      . $FROM . '.';

    my $msg = MIME::Lite->new(
        Subject => "$serverName: Generate Deletion File Submission Details",
        From    => $FROM,
        To      => $email,
        Type    => 'text/html',
        Data    => $EmailContent
    );

    $msg->send( 'smtp', $SERVER, Timeout => 120 );

    print LOG "A submission confirmation email is sent to $email. \n";

}

sub GenDeleteFL_PS_strict {

#
#Li Xue
#Dec 16th, 2012
#
#GenDeleteFL_PS: generate deletion file with constraints on both side (strict condiction):
#A':B' is removed if A and A' share >=30% sequence identity OR B and B' share >=30% sequence identity.
#
#
#
#
#Given a list of protein pairs, generate a list of homo-interologs that share high sequence similarities witht the query pair.
#For A:B, A':B' is written into the output file, if identity(A,A')>=90% OR identity(B,B')>=90%
#
#perl GenDeleteFL.pl prot_pairFL protSeqFL
#perl GenDeleteFL.pl ../data/test/test_protPair.lst ../data/test/test_protSeq.txt
#perl GenDeleteFL.pl ../data/BM3/BM3_pairedChains.txt  ../data/BM3/BM3_pairedChains.protSeq.txt
#perl GenDeleteFL.pl ../data/trans212/trans212pairedChains_onlyTransInterfaces.lst  ../data/trans212/trans212pairedChains_onlyTransInterfaces.protSeq.txt

    use strict;
    use File::Basename;

    #	use myfun;

    # ---------
    our $BlastDIR;    # "$homeDIR/ncbi-blast-2.2.25+";
    our $BlastDataset
      ; #	  "$BlastDIR/data/nr_pdbaa/nr_pdbaa" :non-redundant protein dataset from BLAST FTP database, resolution info added.

    # ---------
    my $protPairFL       = shift @_;
    my $protSeqFL        = shift @_;
    my $IdentityScoreThr = shift @_;    #90
    my $outputFL         = shift @_;    # "$dataDIR/delete.lst";
    my $EvalThr          = 10;

    my $dataDIR = dirname($protPairFL);

    #---------------
    print LOG
"\n\nGenerating deletion file with Identity Threshold: $IdentityScoreThr ...\n\n";

    my @QRYpairs = @{ &readPairFL($protPairFL) };
    my %protSeqs = %{ &readProtSeqFL($protSeqFL) };

    # search for homo-interologys for each test protein chain pair

    foreach my $QRYpair (@QRYpairs) {

        #	$QRYpair='d1h2ka_:d1h2ks_';

        print LOG "\n\n----------------------------------\n\n";
        print LOG "\tNow search homo-interologs for $QRYpair ";  #d1flee_:d1fei_
        print LOG "\n\n----------------------------------\n\n";

        my ( $prot1, $prot2 ) = split( /:/, $QRYpair );
        my @prots = ( $prot1, $prot2 );

        #dir

        #my $pairDIR  = "$dataDIR/$prot1\_$prot2";
        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $dataDIR );
        my $prot1DIR = "$pairDIR/$prot1";
        my $prot2DIR = "$pairDIR/$prot2";

        mkdir($pairDIR)  if ( !-d $pairDIR );
        mkdir($prot1DIR) if ( !-d $prot1DIR );
        mkdir($prot2DIR) if ( !-d $prot2DIR );

        #files
        my $hom_complexesLst       = "$pairDIR/hom-complexes.lst";
        my $seq_int_homComplexesFL = "$pairDIR/seq_int_hom-complexes.lst";
        my $homologFL1 =
          "$pairDIR/$prot1/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl
        my $homologFL2 =
          "$pairDIR/$prot2/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl

        foreach my $prot (@prots) {

            my $QryFstFL =
              "$pairDIR/$prot/$prot.fasta.txt"
              ;    #"$pairDIR/$query1/$pdbID\_$chainID1.fasta.txt";
            my $QryBlastFL = "$pairDIR/$prot/$prot.blast";
            my $homologFL =
              "$pairDIR/$prot/homologsOfquerySeq.lst"
              ;    #the output file of parse_blast.pl
            my $protSeq = $protSeqs{$prot};
            if ( !defined $protSeq ) {

                die(
                    "seq of protein $prot is not available. Check $protSeqFL:$!"
                );
            }
            &writeFastaFL( $prot, $protSeq, $QryFstFL );

            #call Blast and get homolog lists for two queries
            #&callpsiBlast($BlastDataset, $QryFstFL, $QryBlastFL, $EvalThr );
            &callpsiBlast_old( $BlastDataset, $QryFstFL, $QryBlastFL,
                $EvalThr );

            #find hom-complexes from the homolog lists

            &parse_psiblast_userProtDB_globalThr( $QryBlastFL,
                $IdentityScoreThr );

        }

    }

    #--collect homo-interologs into one file

    &collectInterologs_strict( $dataDIR, \@QRYpairs, $outputFL,
        $IdentityScoreThr );

    #--delete other files

    &cleanUp( \@QRYpairs, $dataDIR );

    #--
    print LOG "$outputFL generated.\n";

}

sub GenDeleteFL_PS_loose {

    #
    #Li Xue
    #Dec 16th, 2012
    #

#GenDeleteFL_PS_loose: generate deletion file with constraints on only one side (loose condition):
#A':B' is removed if A and A' share >=30% sequence identity AND B and B' share >=30% sequence identity.

#Modified from PS-HomDDI/GenDeleteFL.pl
#
#
#Given a list of protein pairs, generate a list of homo-interologs that share high sequence similarities witht the query pair.
#For A:B, A':B' is written into the output file, if identity(A,A')>=90% AND identity(B,B')>=90%
#
#perl GenDeleteFL.pl prot_pairFL protSeqFL
#perl GenDeleteFL.pl ../data/trans212/trans212pairedChains_onlyTransInterfaces.lst  ../data/trans212/trans212pairedChains_onlyTransInterfaces.protSeq.txt

    use strict;
    use File::Basename;

    #	use myfun;

    # ---------
    our $BlastDIR;    # "$homeDIR/ncbi-blast-2.2.25+";
    our $BlastDataset
      ; #	  "$BlastDIR/data/nr_pdbaa/nr_pdbaa" :non-redundant protein dataset from BLAST FTP database, resolution info added.

    # ---------
    my $protPairFL       = shift @_;
    my $protSeqFL        = shift @_;
    my $IdentityScoreThr = shift @_;    #90
    my $outputFL         = shift @_;    # "$dataDIR/delete.lst";
    my $EvalThr          = 10;

    my $dataDIR = dirname($protPairFL);

    #---------------
    print LOG
"\n\nGenerating deletion file with Identity Threshold: $IdentityScoreThr ...\n\n";

    my @QRYpairs = @{ &readPairFL($protPairFL) };
    my %protSeqs = %{ &readProtSeqFL($protSeqFL) };

    # search for homo-interologys for each test protein chain pair

    foreach my $QRYpair (@QRYpairs) {

        #	$QRYpair='d1h2ka_:d1h2ks_';

        print LOG "\n\n----------------------------------\n\n";
        print LOG "\tNow search homo-interologs for $QRYpair ";  #d1flee_:d1fei_
        print LOG "\n\n----------------------------------\n\n";

        my ( $prot1, $prot2 ) = split( /:/, $QRYpair );
        my @prots = ( $prot1, $prot2 );

        #dir

        my $pairDIR  = "$dataDIR/$prot1\_$prot2";
        my $prot1DIR = "$pairDIR/$prot1";
        my $prot2DIR = "$pairDIR/$prot2";

        mkdir($pairDIR)  if ( !-d $pairDIR );
        mkdir($prot1DIR) if ( !-d $prot1DIR );
        mkdir($prot2DIR) if ( !-d $prot2DIR );

        #files
        my $hom_complexesLst = "$pairDIR/hom-complexes.lst";

     #        my $seq_int_homComplexesFL = "$pairDIR/seq_int_hom-complexes.lst";
        my $homologFL1 =
          "$pairDIR/$prot1/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl
        my $homologFL2 =
          "$pairDIR/$prot2/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl

        foreach my $prot (@prots) {

            my $QryFstFL =
              "$pairDIR/$prot/$prot.fasta.txt"
              ;    #"$pairDIR/$query1/$pdbID\_$chainID1.fasta.txt";
            my $QryBlastFL = "$pairDIR/$prot/$prot.blast";
            my $homologFL =
              "$pairDIR/$prot/homologsOfquerySeq.lst"
              ;    #the output file of parse_blast.pl
            my $protSeq = $protSeqs{$prot};
            if ( !defined $protSeq ) {

                die(
                    "seq of protein $prot is not available. Check $protSeqFL:$!"
                );
            }
            &writeFastaFL( $prot, $protSeq, $QryFstFL );

            #call Blast and get homolog lists for two queries

            &callpsiBlast_old( $BlastDataset, $QryFstFL, $QryBlastFL,
                $EvalThr );

            #find hom-complexes from the homolog lists

            &parse_psiblast_userProtDB_globalThr( $QryBlastFL,
                $IdentityScoreThr );

        }

        #compare two homolog lists
        #if they have the same pdbID and different chainID
        #then this pdb complex is written into the output file

        &compareHomlogLsts( $homologFL1, $homologFL2 )
          ;    #outputFL: ../data/hom-complexes.lst

    }

    #--collect homo-interologs into one file

    &collectInterologs_loose( $dataDIR, \@QRYpairs, $outputFL,
        $IdentityScoreThr );

    #--delete other files

    &cleanUp( \@QRYpairs, $dataDIR );

}

sub GenDeleteFL_PS_loose_species {

    #
    #Li Xue
    #Dec 16th, 2012
    #

#GenDeleteFL_PS_loose: generate deletion file with constraints on only one side (loose condition):
#A':B' is removed if A and A' share >=30% sequence identity AND B and B' share >=30% sequence identity.

#Modified from PS-HomDDI/GenDeleteFL.pl
#
#
#Given a list of protein pairs, generate a list of homo-interologs that share high sequence similarities witht the query pair.
#For A:B, A':B' is written into the output file, if [identity(A,A')>=90% and isSameSpecies(A,A')] AND [identity(B,B')>=90% and isSameSpecies(B,B')]
#
#perl GenDeleteFL.pl prot_pairFL protSeqFL
#perl GenDeleteFL.pl ../data/test/test_protPair.lst ../data/test/test_protSeq.txt
#perl GenDeleteFL.pl ../data/BM3/BM3_pairedChains.txt  ../data/BM3/BM3_pairedChains.protSeq.txt
#perl GenDeleteFL.pl ../data/trans212/trans212pairedChains_onlyTransInterfaces.lst  ../data/trans212/trans212pairedChains_onlyTransInterfaces.protSeq.txt

    use strict;
    use File::Basename;

    #	use myfun;

    # ---------
    our $BlastDIR;    # "$homeDIR/ncbi-blast-2.2.25+";
    our $BlastDataset
      ; #	  "$BlastDIR/data/nr_pdbaa/nr_pdbaa" :non-redundant protein dataset from BLAST FTP database, resolution info added.

    # ---------
    my $protPairFL       = shift @_;
    my $protSeqFL        = shift @_;
    my $IdentityScoreThr = shift @_;    #90
    my $taxonomyMap = shift @_;    #my $taxonomyMap = &readP_chain_taxonomyLST;
    my $outputFL    = shift @_;    # "$dataDIR/delete.lst";

    my $EvalThr = 10;

    my $dataDIR =
      dirname($protPairFL)
      ; #....\predicted_int_PSHomPPI_5angstroms_SeqIdenThr90_strict_newPSHomPPI\1BKD\predicted_int\R_S

    #---------------
    print LOG
"\n\nGenerating deletion file with Identity Threshold: $IdentityScoreThr ...\n\n";

    my @QRYpairs = @{ &readPairFL($protPairFL) };
    my %protSeqs = %{ &readProtSeqFL($protSeqFL) };

#--qry pdb ID chn ID map file. Use to determine the species of the query protein
    my $caseDIR = dirname($dataDIR);

    my $qryPDBIDchnIDmapFL = "$caseDIR/uniqueID2pdbIDchnID.map";

    my $pdbIDmap = &readQryPDBIDchnIDmapFL($qryPDBIDchnIDmapFL);

    # search for homo-interologys for each test protein chain pair

    foreach my $QRYpair (@QRYpairs) {

        #	$QRYpair='d1h2ka_:d1h2ks_';

        print LOG "\n\n----------------------------------\n\n";
        print LOG "\tNow search homo-interologs for $QRYpair ";  #d1flee_:d1fei_
        print LOG "\n\n----------------------------------\n\n";

        my ( $prot1, $prot2 ) = split( /:/, $QRYpair );
        my @prots = ( $prot1, $prot2 );

        #dir

        my $pairDIR  = "$dataDIR/$prot1\_$prot2";
        my $prot1DIR = "$pairDIR/$prot1";
        my $prot2DIR = "$pairDIR/$prot2";

        mkdir($pairDIR)  if ( !-d $pairDIR );
        mkdir($prot1DIR) if ( !-d $prot1DIR );
        mkdir($prot2DIR) if ( !-d $prot2DIR );

        #files
        my $hom_complexesLst       = "$pairDIR/hom-complexes.lst";
        my $seq_int_homComplexesFL = "$pairDIR/seq_int_hom-complexes.lst";
        my $homologFL1 =
          "$pairDIR/$prot1/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl
        my $homologFL2 =
          "$pairDIR/$prot2/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl

        foreach my $prot (@prots) {

            my $QryFstFL =
              "$pairDIR/$prot/$prot.fasta.txt"
              ;    #"$pairDIR/$query1/$pdbID\_$chainID1.fasta.txt";
            my $QryBlastFL = "$pairDIR/$prot/$prot.blast";
            my $homologFL =
              "$pairDIR/$prot/homologsOfquerySeq.lst"
              ;    #the output file of parse_blast.pl
            my $protSeq = $protSeqs{$prot};
            if ( !defined $protSeq ) {

                die(
                    "seq of protein $prot is not available. Check $protSeqFL:$!"
                );
            }
            &writeFastaFL( $prot, $protSeq, $QryFstFL );

            #call Blast and get homolog lists for two queries

            &callpsiBlast_old( $BlastDataset, $QryFstFL, $QryBlastFL,
                $EvalThr );

            #find hom-complexes from the homolog lists

            my @homologs = @{
                &parse_psiblast_userProtDB_globalThr( $QryBlastFL,
                    $IdentityScoreThr )
            };

      #-check whether each homolog is from the same species as the query protein
            my $qryPDBIDchnID = $pdbIDmap->{$prot};

            if ( !defined $qryPDBIDchnID ) {

                die("Qry $prot is not defined in $qryPDBIDchnIDmapFL:$!");
            }
            my @sameSpecieshomologs =
              &getSameSpeciesHomologs( $qryPDBIDchnID, \@homologs,
                $taxonomyMap );

            &writeArrayIntoFL( \@sameSpecieshomologs, $homologFL );

        }

        #compare two homolog lists
        #if they have the same pdbID and different chainID
        #then this pdb complex is written into the output file

        &compareHomlogLsts( $homologFL1, $homologFL2 )
          ;    #outputFL: ../data/hom-complexes.lst

    }

    #--collect homo-interologs into one file

    &collectInterologs_loose_species( $dataDIR, \@QRYpairs, $outputFL,
        $IdentityScoreThr );

    #--delete other files

    &cleanUp( \@QRYpairs, $dataDIR );

}

sub getSameSpeciesHomologs {

    my $qryPDBIDchnID = shift @_;
    my @homologs      = @{ shift @_ };
    my $taxonomyMap   = shift @_;

    my $qryPDBID = substr( $qryPDBIDchnID, 0, 4 );
    my $qrychnID = substr( $qryPDBIDchnID, 4, 1 );

    my @sameSpeciesHomologs;

    foreach (@homologs) {
        my $pdbID_homolog = substr( $_, 0, 4 );
        my $chnID_homolog = substr( $_, 4, 1 );

        my $ans1 = &isSamespecies(
            $taxonomyMap,   $qryPDBID, $qrychnID,
            $pdbID_homolog, $chnID_homolog
        );

        if ( ( $ans1 == 1 || $ans1 == 2 ) ) {

            print LOG
              "qry $qryPDBIDchnID and homolog $_ are from the same species.\n";

            push @sameSpeciesHomologs, $_;
        }
    }

    return \@sameSpeciesHomologs;

}

sub readPairFL {
    my $protPairFL = shift @_;
    my @pairs;

    open( INPUT, "<$protPairFL" ) || die("Cannot open $protPairFL:$!");

    while (<INPUT>) {

        s/[\n\r]//mg;
        if (/^\w+:\w+/) {
            push @pairs, $_;

        }
    }
    close INPUT;
    if ( !@pairs ) {
        die("$protPairFL does not have pairs:$!");
    }

    return \@pairs;
}

sub readProtSeqFL {

    #	>Actin
    #	MDSEVAALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGIMVGMGQKDSYVGDEAQS
    #	KRGILTLRYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPMNPKSNREKMT
    #	QIMFETFNVPAFYVSIQAVLSLYSSGRTTGIVLDSGDGVTHVVPIYAGFSLPHAILRIDL
    #	AGRDLTDYLMKILSERGYSFSTTAEREIVRDIKEKLCYVALDFEQEMQTAAQSSSIEKSY
    #	ELPDGQVITIGNERFRAPEALFHPSVLGLESAGIDQTTYNSIMKCDVDVRKELYGNIVMS
    #	GGTTMFPGIAERMQKEITALAPSSMKVKIIAPPERKYSVWIGGSILASLTTFQQMWISKQ
    #	EYDESGPSIVHHKCF

    my $protSeqFL = shift @_;
    my %seqs;
    my $ID;
    my $seq;

    open( INPUT, "<$protSeqFL" ) || die("Cannot open $protSeqFL:$!");
    foreach (<INPUT>) {
        s/[\n\r]//mg;
        if (/^\s*$/) {
            next;
        }
        if (/^>(\w+)/) {

            if ( defined $ID && defined $seq ) {
                $seqs{$ID} = $seq;

                #                print LOG "$ID:$seq\n";

            }

            $seq = '';
            $ID  = $1;
            next;
        }
        if (/^[a-zA-Z]+/) {
            s/\s+//;
            $seq = $seq . $_;

        }

    }

    #read the last protein
    if ( defined $ID && defined $seq ) {
        $seqs{$ID} = $seq;

        #        print LOG "$ID:$seq\n";
    }

    close INPUT;

    if ( !%seqs ) {
        die("Check $protSeqFL:$!");
    }

    return \%seqs;

}

sub writeFastaFL {
    my $ID       = shift @_;
    my $seq      = shift @_;
    my $outputFL = shift @_;

    unlink $outputFL if ( -e $outputFL );

    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");

    print OUTPUT ">$ID\n";
    print OUTPUT "$seq\n";
    close OUTPUT;

    print LOG "$outputFL generated.\n";

}

sub parse_psiblast_userProtDB_globalThr {

#
#This script is part of Protein-Protein interface prediction using customized training dataset
#parse pdb ID and chain No. from the BLAST result file, and write them into files, "homologsOfquerySeq.lst"
#user can specify the sequence similarity threshold (Positive Score threshold)
#
#note: This script deals with the format when Blast program blasts against user-defined proteinProtein database.
#
#
#1. Input file format (BLAST output file):
#
#   HHPRED&PBLAST
#   # Iteration: 1
#   # Query: Q
#   # Database: /home/mwalter/project/v0.0.02/projectfiles/databases/pdb70HMMs
#   # Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, HHpredprob, query acc.
#   # 85 hits found
#
#   gi|530537729|pdb|4LQW|A164163VNPTVFFDIAVDGEPLGRVSFELFADKVPKTAENFRALSTGEKGFGYKGSCFHRIIPGFMCQGGDFTRHNGTGGKSIYGEKFEDENFILKHTGPGILSMANAGPNTNGSQFFICTAKTEWLDGKHVVFGKVKEGMNIVEAMERFGSRNGKTSKKITIADCGQLTNPVVFFDVCADGEPLGRITMELFSNIVPRTAENFRALCTGEKGFGFKNSIFHRVIPDFVCQGGDITKHDGTGGQSIYGDKFEDENFDVKHTGPGLLSMANQGQNTNNSQFVITLKKAEHLDFKHVVFGFVKDGMDTVKKIESFGSPKGSVCRRITITECGQI016314177285.21.3E-4716366.25NaN100.0Q
#
#
#2. When chain ID is in lowercase, Blast will output double uppercase chainID
# For example: pdb|3KIX|LL  Chain l,  >pdb|2QNH|MM Chain m,  >pdb|3KIS|LL Chain ...   259    3e-70
# So we need to extract chain ID from "Chain l" part not from "pdb|3KIX|LL".
#
#

#
# Usage: perl parse_psiblast_userProtDB.pl filename.blast PositiveSThreshold IdentityThreshold
# e.g., perl parse_psiblast_userProtDB.pl ../data/test/test.blast 0  90

    use strict;
    use File::Basename;

    #Variables
    my $inputfile    = shift @_;    # BLAST output file: Q.blast
    my $identitySThr = shift @_;
    my $taskID = basename( dirname($inputfile) );
    my $inputfilename = basename( $inputfile, '.blast' );
    my $PDBID;
    my $chainID;
    my $roundNum;                   #the number of iteration of PSI-BLAST
    my $homologRef;
    my $homologID;
    my @homlog
      ; #to store the PDBID and chainID of homologs of query sequence returned Blastp
    my $similarity;
    my $num_seq = 0;    #number of the homologs for each query sequence
    my $positiveSThr = 0;

    #directories
    print LOG
"\n\nParse BLAST output file ($inputfile) with the following thresholds:\n";
    print LOG "\tGLOBAL Positive Score Thr: $positiveSThr\n";
    print LOG "\tGLOBAL Identity Score Thr: $identitySThr\n\n\n";

    my $outputDIR   = dirname($inputfile);
    my $outfilename = "homologsOfquerySeq.lst";    #define output file name
    my $outfile     = "$outputDIR/$outfilename";
    unlink $outfile if ( -e $outfile );

    #------------- program begins --------------------

#read the blast output file. Only last round of PSI-blast results are returned in $blastResults.
    my ( $qryID, $qryLength, $blastResults ) =
      &readBlastOUTPUTFL_userProtDB($inputfile);

    if ( !defined $qryID && !defined $qryLength && !defined $blastResults ) {
        print LOG "No hit found in $inputfile\n";
        @$homologRef = ();
        &writeHomologLstfile( $outfile, $homologRef, $positiveSThr );
        print LOG "parse_psiblast_userDB.pl finish!  Outputfile: $outfile\n\n";
        return;
    }

    #--get homologs that have PositiveS and IdentityS >= thresholds
    #print LOG "\nPassed homologs:\n";
    #    my @homologs = keys(%$blastResults);

    foreach my $homolog ( keys %{$blastResults} ) {

        foreach my $num ( keys %{ $blastResults->{$homolog} } ) {
            my $identityS_local =
              $blastResults->{$homolog}->{$num}->{'identityS'};    #[0-100]
            my $homolog_Len = $blastResults->{$homolog}->{$num}->{'slen'};
            my $LAL         = $blastResults->{$homolog}->{$num}->{'LAL'};

            my $num_identicalResi = $identityS_local * $LAL * 0.01;

            print LOG
"$homolog: Num_identicalResi/Qry_Len = $num_identicalResi/$qryLength and Num_identicalResi/Homolog_Len = $num_identicalResi/$homolog_Len\n";

            if (   $num_identicalResi / $qryLength * 100 >= $identitySThr
                && $num_identicalResi / $homolog_Len * 100 >= $identitySThr )
            {    #check the sequece similarity
                push @$homologRef, $homolog;
            }

        }
    }

#write the homologs that are in the last iteration of PSI-BLAST and pass the positiveS threshold
#to the output file

    if ( !$homologRef ) {
        print LOG
          "No homolog is returned by BLAST pass the Identity thresholds.\n";
        @$homologRef = ();
    }

    my @homologs_final = &unique(@$homologRef);
    &writeHomologLstfile( $outfile, $homologRef, $positiveSThr, $identitySThr );

    print LOG "parse_psiblast_userProtDB{} finish!  Outputfile: $outfile\n\n";
}

#--------------------------

sub parse_psiblast_userProtDB_globalThr_old {

#Author: Xue, Li
#Date: 12/18/2012
#
#
#parse_psiblast_userProtDB_globalThr: use GLOBAL sequence identity threshold
#parse_psiblast_userProtDB: use LOCAL sequence identity threshold
#
#
#This script is part of Protein-Protein interface prediction using customized training dataset
#parse pdb ID and chain No. from the BLAST result file, and write them into files, "homologsOfquerySeq.lst"
#user can specify the sequence similarity threshold (Positive Score threshold)
#
#note: This script deals with the format when Blast program blasts against user-defined proteinProtein database.
#
#
#1. Input file format (BLAST output file):
#
#	 # PSIBLAST 2.2.25+
#	# Iteration: 1
#	# Query: 1avxA
#	# Database: C:\Users\lixue/ncbi-blast-2.2.25+/data/scop_v175_prot/scop_v175_prot
#	# Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, query acc.
#	# 485 hits found
#	pdb|1tx6|D	223	223	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	1	223	1	223	 413	2e-116	223	100.00	100.00	1avxA
#	pdb|1eja|A	223	223	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	1	223	1	223	 413	2e-116	223	100.00	100.00	1avxA
#	pdb|1tx6|A	223	223	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	1	223	1	223	 413	2e-116	223	100.00	100.00	1avxA
#
#
#
#
#2. When chain ID is in lowercase, Blast will output double uppercase chainID
# For example: pdb|3KIX|LL  Chain l,  >pdb|2QNH|MM Chain m,  >pdb|3KIS|LL Chain ...   259    3e-70
# So we need to extract chain ID from "Chain l" part not from "pdb|3KIX|LL".
#
#

#
# Usage: perl parse_psiblast_userProtDB.pl filename.blast PositiveSThreshold IdentityThreshold
# e.g., perl parse_psiblast_userProtDB.pl ../data/test/test.blast 0  90

    use strict;
    use File::Basename;

    #Variables
    my $inputfile = shift @_;

    #	my $positiveSThr  = shift @_;
    my $identitySThr  = shift @_;
    my $taskID        = basename( dirname($inputfile) );
    my $inputfilename = basename( $inputfile, '.blast' );
    my $PDBID;
    my $chainID;
    my $roundNum;    #the number of iteration of PSI-BLAST
    my $homologRef;
    my $homologID;
    my @homlog
      ; #to store the PDBID and chainID of homologs of query sequence returned Blastp
    my $similarity;
    my $num_seq = 0;    #number of the homologs for each query sequence

    #directories
    #use Cwd;
    #my $DIR     = cwd;

    my $positiveSThr = 0;

    print LOG
"\n\nParse BLAST output file ($inputfile) with the following thresholds:\n";
    print LOG "\tGLOBAL Positive Score Thr: $positiveSThr\n";
    print LOG "\tGLOBAL Identity Score Thr: $identitySThr\n\n\n";

    my $outputDIR   = dirname($inputfile);
    my $outfilename = "homologsOfquerySeq.lst";    #define output file name
    my $outfile     = "$outputDIR/$outfilename";
    unlink $outfile if ( -e $outfile );

    #program begins

#read the blast output file. Only last round of PSI-blast results are returned in $blastResults.
    my ( $qryID, $qryLength, $blastResults ) =
      &readBlastOUTPUTFL_userProtDB($inputfile);

    if ( !defined $qryID && !defined $qryLength && !defined $blastResults ) {
        print LOG "No hit found in $inputfile\n";
        @$homologRef = ();
        &writeHomologLstfile( $outfile, $homologRef, $positiveSThr,
            $identitySThr );

        print LOG
"parse_psiblast_userProtDB_globalThr() finish!  Outputfile: $outfile\n\n";
        return;

    }

    #get homologs that have IdentityS >= thresholds
    my @homologs = keys(%$blastResults);

    foreach my $homolog (@homologs) {
        my $identityS_local = $blastResults->{$homolog}->{'identityS'}; #[0-100]
        my $homolog_Len     = $blastResults->{$homolog}->{'slen'};
        my $LAL             = $blastResults->{$homolog}->{'LAL'};

        my $num_identicalResi = $identityS_local * $LAL * 0.01;

        print LOG
"$homolog: Num_identicalResi/Qry_Len = $num_identicalResi/$qryLength and Num_identicalResi/Homolog_Len = $num_identicalResi/$homolog_Len\n";

        #xuexuexuexue

        if (   $num_identicalResi / $qryLength * 100 >= $identitySThr
            && $num_identicalResi / $homolog_Len * 100 >=
            $identitySThr )    #check the sequece similarity
        {
            push @$homologRef, $homolog;

            #print LOG "$homolog passed.\n"; #xuexuexuexue
        }
        else {

            #print LOG "$homolog NOT passed.\n";#xuexuexuexue
        }

    }

#write the homologs that are in the last iteration of PSI-BLAST and pass the positiveS threshold
#to the output file

    if ( !defined $homologRef ) {
        print LOG
          "No homolog is returned by BLAST pass the Identity thresholds.\n";
        @$homologRef = ();
    }

    &writeHomologLstfile( $outfile, $homologRef, $positiveSThr, $identitySThr );

    print LOG "parse_psiblast_userDB.pl finish!  Outputfile: $outfile\n\n";

    return $homologRef;
}

sub collectInterologs_strict {

    #collect homo-interologs into one file

    #Input files (homologsOfquerySeq.lst for Qry A and B):

    #		3areD
    #		3mnzA
    #		4d9lJ
    #		1q9lA
    #		3mlzL
    #		3f7vB
    #		3q6fK
    #		3tygL
    #		1rzgB

    #Output file (delete.lst):

    #		A:C => 1jpsL:*,*:1ahwC,1ahwA:1ahwF,1uj3A:*,1ahwD:1ahwC,1ahwD:1ahwF
    #		B:C => 1uj3B:1uj3C,*:1jpsT,1ahwE:1ahwF,*:1ahwC,1ahwE:1ahwC,1ahwB:1ahwF

#3i29A:* meaning any homologous protein complex with 3i29A on the left side should be deleted.

    my $dataDIR        = shift @_;
    my @QRYpairs       = @{ shift @_ };
    my $outputFL       = shift @_;
    my $SeqIdentityThr = shift @_;

    &headerDelFL_strict( $outputFL, $SeqIdentityThr );

    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    foreach my $pair (@QRYpairs) {

        my ( $prot1, $prot2 ) = split( /:/, $pair );
        my @prots = ( $prot1, $prot2 );

        #--
        #		my $pdbID1 = lc( substr( $prot1, 0, 4 ) );
        #		my $chnID1 = substr( $prot1, 4, 1 );
        #		my $pdbID2 = lc( substr( $prot2, 0, 4 ) );
        #		my $chnID2 = substr( $prot2, 4, 1 );
        #		$pair = "$pdbID1$chnID1:$pdbID2$chnID2";

        #dir

        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $pair, $dataDIR );

        if ( !-d $pairDIR ) {
            die("$pairDIR does not exist:$!");
        }

        my $homologsOfquerySeqLST1 = "$pairDIR/$prot1/homologsOfquerySeq.lst";
        my $homologsOfquerySeqLST2 = "$pairDIR/$prot2/homologsOfquerySeq.lst";

        my @homologs1 = @{ &readLstFL($homologsOfquerySeqLST1) };
        my @homologs2 = @{ &readLstFL($homologsOfquerySeqLST2) };

        my $homointerolog;
        if ( !@homologs1 && !@homologs2 ) {

            print LOG "\n\n**** No homo-interologs found for $pair.**** \n\n";
            $homointerolog = '';
        }

        my @homointerologs_new;    #('3i29A:*','2qyiC:*')

        if (@homologs1) {

            foreach my $homolog (@homologs1) {

                #				$homolog= '3fclB'

                push @homointerologs_new, "$homolog:*";

            }
        }

        if (@homologs2) {

            foreach my $homolog (@homologs2) {

                push @homointerologs_new, "*:$homolog";

            }
        }

        #--

        $homointerolog = join( ',', @homointerologs_new );

        print OUTPUT "$pair => $homointerolog\n";

    }
    close OUTPUT;

    print LOG "\n$outputFL generated.\n";

}

sub collectInterologs_loose {

    #collect homo-interologs into one file

#output:
#Highly similar homo-interologs.
#For qry A:B, A':B' is a highly similar homo-interolog if A and A' share >= 90% identity AND B and B' share >= 90% identity.
#A:B => 1acbE:1acbI

    my $dataDIR        = shift @_;
    my @QRYpairs       = @{ shift @_ };
    my $outputFL       = shift @_;
    my $SeqIdentityThr = shift @_;

    &headerDelFL_loose( $outputFL, $SeqIdentityThr );

    open( OUTPUT, ">>$outputFL" );

    foreach my $pair (@QRYpairs) {

        my ( $prot1, $prot2 ) = split( /:/, $pair );
        my @prots = ( $prot1, $prot2 );

        #--
        #		my $pdbID1 = lc( substr( $prot1, 0, 4 ) );
        #		my $chnID1 = substr( $prot1, 4, 1 );
        #		my $pdbID2 = lc( substr( $prot2, 0, 4 ) );
        #		my $chnID2 = substr( $prot2, 4, 1 );
        #		$pair = "$pdbID1$chnID1:$pdbID2$chnID2";

        #dir

        my $pairDIR = "$dataDIR/$prot1\_$prot2";

        my $homInterologFL = "$pairDIR/templates.lst";

        print LOG "Read $homInterologFL ..\n";
        my $homointerolog;
        if ( !-e $homInterologFL ) {
            print LOG "\n\n**** $homInterologFL does not exist! No homo-interologs found for $pair.**** \n\n";
            $homointerolog = '';
        }
        else {
            my @homointerologs =
              @{ &readHomComplexFL($homInterologFL) }
              ;    #('3i29    A:B','2qyi    C:D')

            #--
            my @homointerologs_new;    #('3i29A:3i29B','2qyiC:2qyiD')
            foreach my $a (@homointerologs) {

                #$a ='2qyi    C:D'
                my ( $pdbID, $chn1, $chn2 ) = split( /[\s\t:]+/, $a );

                push @homointerologs_new, "$pdbID$chn1:$pdbID$chn2";

            }

            #--

            $homointerolog = join( ',', @homointerologs_new );
        }

        print OUTPUT "$pair => $homointerolog\n";

    }
    close OUTPUT;

    print LOG "\n$outputFL generated.\n";

}

sub collectInterologs_loose_species {

    #collect homo-interologs into one file

#output:
#Highly similar homo-interologs.
#For qry A:B, A':B' is a highly similar homo-interolog if A and A' share >= 90% identity AND B and B' share >= 90% identity.
#A:B => 1acbE:1acbI

    my $dataDIR  = shift @_;
    my @QRYpairs = @{ shift @_ };
    my $outputFL = shift @_;

    our $SeqIdentityThr;

    &headerDelFL_loose_species( $outputFL, $SeqIdentityThr );

    open( OUTPUT, ">>$outputFL" );

    foreach my $pair (@QRYpairs) {

        my ( $prot1, $prot2 ) = split( /:/, $pair );
        my @prots = ( $prot1, $prot2 );

        #--
        #		my $pdbID1 = lc( substr( $prot1, 0, 4 ) );
        #		my $chnID1 = substr( $prot1, 4, 1 );
        #		my $pdbID2 = lc( substr( $prot2, 0, 4 ) );
        #		my $chnID2 = substr( $prot2, 4, 1 );
        #		$pair = "$pdbID1$chnID1:$pdbID2$chnID2";

        #dir

        my $pairDIR = "$dataDIR/$prot1\_$prot2";

        my $homInterologFL = "$pairDIR/templates.lst";

        print LOG "Read $homInterologFL ..\n";
        my $homointerolog;
        if ( !-e $homInterologFL ) {
            print LOG "\n\n**** $homInterologFL does not exist! No homo-interologs found for $pair.**** \n\n";
            $homointerolog = '';
        }
        else {
            my @homointerologs =
              @{ &readHomComplexFL($homInterologFL) }
              ;    #('3i29    A:B','2qyi    C:D')

            #--
            my @homointerologs_new;    #('3i29A:3i29B','2qyiC:2qyiD')
            foreach my $a (@homointerologs) {

                #$a ='2qyi    C:D'
                my ( $pdbID, $chn1, $chn2 ) = split( /[\s\t:]+/, $a );

                push @homointerologs_new, "$pdbID$chn1:$pdbID$chn2";

            }

            #--

            $homointerolog = join( ',', @homointerologs_new );
        }

        print OUTPUT "$pair => $homointerolog\n";

    }
    close OUTPUT;

    print LOG "\n$outputFL generated.\n";

}

sub headerDelFL_loose {

    my $delFL          = shift @_;
    my $SeqIdentityThr = shift @_;

    unlink $delFL if ( -e $delFL );
    open( OUTPUT, ">>$delFL " ) || die("Cannot open $delFL:$!");
    print OUTPUT "#GenDeleteFL_PS_loose() in PSHomPPI_resiPairs.pm\n";
    print OUTPUT "#Highly similar homo-interologs. \n";
    print OUTPUT
"#For qry A:B, A':B' is a highly similar homo-interolog if A and A' share >= $SeqIdentityThr% identity AND B and B' share >= $SeqIdentityThr% identity. \n";
    close OUTPUT;

}

sub headerDelFL_loose_species {

    my $delFL          = shift @_;
    my $SeqIdentityThr = shift @_;

    unlink $delFL if ( -e $delFL );
    open( OUTPUT, ">>$delFL " ) || die("Cannot open $delFL:$!");
    print OUTPUT "#GenDeleteFL_PS_loose() in PSHomPPI_resiPairs.pm\n";
    print OUTPUT "#Highly similar homo-interologs. \n";
    print OUTPUT
"#For qry A:B, A':B' is a highly similar homo-interolog if [A and A' share >= $SeqIdentityThr% identity and isSameSpecies(A,A')] AND [B and B' share >= $SeqIdentityThr% identity and isSmesSpecies(B,B')] . \n";
    close OUTPUT;

}

sub headerDelFL_strict {

    my $delFL          = shift @_;
    my $SeqIdentityThr = shift @_;

    unlink $delFL if ( -e $delFL );
    open( OUTPUT, ">>$delFL " ) || die("Cannot open $delFL:$!");
    print OUTPUT "#GenDeleteFL_PS_strict() in PSHomPPI_resiPairs.pm\n";
    print OUTPUT "#Highly similar homo-interologs. \n";
    print OUTPUT
"#For qry A:B, A':B' is a highly similar homo-interolog if A and A' share >= $SeqIdentityThr% identity OR B and B' share >= $SeqIdentityThr% identity. \n";
    close OUTPUT;

}

sub cleanUp {

    #-- delete the by-product files generated by GenDeleteFL_PS().
    use strict;

    print LOG "\n *** Now cleaning up ... ***\n";

    my @QRYpairs = @{ shift @_ };
    my $dataDIR  = shift @_;

    foreach my $QRYpair (@QRYpairs) {

        #	$QRYpair='d1h2ka_:d1h2ks_';

        my ( $prot1, $prot2 ) = split( /:/, $QRYpair );
        my @prots = ( $prot1, $prot2 );

        #dir

        my $pairDIR  = "$dataDIR/$prot1\_$prot2";
        my $prot1DIR = "$pairDIR/$prot1";
        my $prot2DIR = "$pairDIR/$prot2";

        #files
        #		my $QryBlastFL1            = "$pairDIR/$prot1/$prot1.blast";
        #		my $QryBlastFL2            = "$pairDIR/$prot2/$prot2.blast";
        my $hom_complexesLst       = "$pairDIR/hom-complexes.lst";
        my $seq_int_homComplexesFL = "$pairDIR/seq_int_hom-complexes.lst";
        my $homologFL1 =
          "$pairDIR/$prot1/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl
        my $homologFL2 =
          "$pairDIR/$prot2/homologsOfquerySeq.lst"
          ;    #the output file of parse_blast.pl

        my $homologs_in_homComplexesFL1 =
          "$pairDIR/$prot1/homologs_in_homComplexes.lst";
        my $homologs_in_homComplexesFL2 =
          "$pairDIR/$prot2/homologs_in_homComplexes.lst";

        my $homologs_globalIdentity1 =
          "$pairDIR/$prot1/homologsOfquerySeq_globalIdentity.lst";
        my $homologs_globalIdentity2 =
          "$pairDIR/$prot2/homologsOfquerySeq_globalIdentity.lst";
        my $homologs_localIdentity1 =
          "$pairDIR/$prot1/homologsOfquerySeq_localIdentity.lst";
        my $homologs_localIdentity2 =
          "$pairDIR/$prot2/homologsOfquerySeq_localIdentity.lst";

        my @filesToBeDel = (
            $homologFL1,                  $homologFL2,
            $hom_complexesLst,            $seq_int_homComplexesFL,
            $homologs_in_homComplexesFL1, $homologs_in_homComplexesFL2,
            $homologs_globalIdentity1,    $homologs_globalIdentity2,
            $homologs_localIdentity1,     $homologs_localIdentity2
        );

        foreach (@filesToBeDel) {

            unlink $_ if ( -e $_ );
            print LOG "$_ deleted.\n";
        }

    }

    print LOG
      "the by-product files generated by GenDeleteFL_PS() are removed.\n";

}

sub readHomComplexFL {

    #Input FL:
    #4b2c    A:D
    #4b2c    A:B
    #3rdz    B:C
    #3rdz    B:D

    my $FL = shift @_;
    my @elements;

    open( INPUT, "<$FL" ) || die("Cannot open $FL:$!");

    foreach (<INPUT>) {
        s/[\n\r]//mg;
        if (/^\w+/) {
            push @elements, $_;
        }
    }
    close INPUT;

    if ( !@elements ) {
        print LOG "$FL is empty.\n";
        @elements = ();
    }

    return \@elements;
}


1;
