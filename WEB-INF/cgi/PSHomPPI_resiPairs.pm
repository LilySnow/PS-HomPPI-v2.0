
#===============================================================================
#
#        Li Xue ( me.lixue@gmail.com )
#        12/02/2016 10:16:40 AM
#        Utrecht University
#
#  DESCRIPTION: functions called by PSHomPPI.cgi
#
#===============================================================================

use strict;
use warnings;

use globalVariables;
our $searchNewPPIDBDIR;           # 'searchNewPPIDB';
our $safe_filename_characters;    # "a-zA-Z0-9_\.-";

sub batch_PSHomPPI {

# Li Xue
# Jul. 13th, 2012
# modified: Feb. 6th 2015
#
# batch.pl : predict proteins using safe/twi/dark homologs
# 1. predict qrys with safe zone homologs
# 2. if qrys do not have safe zone homologs, predict qrys with twilight zone homologs
# 3. if qrys do not have safe and twilight zone homologs,predict qrys with dark zone homologs

# collect all the homologs' info into a file of statistics.txt.
# Apply similarity thresholds to decided which homologs in this file will be used for prediction.
#
#
#Input files:
#
#1. Protein ID pair list
#	prot1:prot2
#	prot3:prot4
#2. seq file of query proteins
#	>prot1
#	seq_prot1
#	>prot2
#	seq_prot2
#3. a list of homo-interologs to be removed (generate by PS-HomDDI/GenDeleteFL.pl)
#	1a2kA:1a2kC =>
#	1a14H:1a14N => 1NMB_H:1NMB_N,1A14_H:1A14_N,1NMA_H:1NMA_N
#4. seq_int file of query prot pairs (optional)
#   >prot1|prot1:prot2
#   seq_prot1
#   interface of prot1 with prot2
#
#
#
# predictionMode: SafeMode, TwilightMode, DarkMode
# SafeMode: only use Safe-Zone homologs
# TwilightMode: use Safe and Twilight-Zone homologs
# DarkMode: use Safe, Twilight and Dark-Zone homologs
# IntDef: atomDistances  vanDerWaalsDistances
# Usage: perl batch.pl FLAG inputFL $k IntDef atomDistThr rasaThr scoreThr

    use strict;
    use File::Basename;
    use File::Copy;
    use File::Path;

    #------------
    #global variables
    our $perlDIR;              # = '/usr/bin/perl';
    our $scoreThr;             #prediction score cutoff value
    our $searchNewPPIDBDIR;    #    = "$NPSHomPPIDIR/searchNewPPIDB";

    #------------
    #input variables

    my $jobDIR        = shift @_;
    my $qryPairFL     = shift @_;    # "$jobDIR/qryIDpair.lst";
    my $protSeqFL     = shift @_;
    my $delLstFL      = shift @_;
    my $seqInt_protFL = shift @_;
    my $intDef        = shift @_;    #'vanDerWaalsDistances';
    my $atomDistThr   = shift @_;    #8.5
    my $rasaThr       = shift @_;    #0

    #------------
    #local variables

    my $flag_prediction = 0;
    ; #$flag_prediction=0: at least one of the protein pairs have homo-interologs in $FullStatFL
      #$flag_prediction=1: none of the query pairs have homo-interologs in $FullStatFL

    my $rmSameProtFlag;
    if ( -e $delLstFL ) {

        $rmSameProtFlag = 'Yes';
    }
    else {
        $rmSameProtFlag = 'No';
    }

    my $jobID = basename( $qryPairFL, ( '.txt', '.lst' ) );

    #------------
    #file variables

    my $seq_int_qryComplexesFL = "$jobDIR/seq_int_qryComplexes.lst";
    my $statisticsFL           = "$jobDIR/statistics.txt";

    #    my $QryWoHomComplexLST =
    #      "$jobDIR/qryWoTemplates.lst"
    #      ;    #keep a list of Qry pairs that cannot find homo-complexes
    #    unlink $QryWoHomComplexLST if ( -e $QryWoHomComplexLST );

    my $FullStatFL;
    if ( $rmSameProtFlag eq 'Yes' ) {

        $FullStatFL =
          "$jobDIR/statistics_wo_sameProt.txt"
          ;    #../data/benchmark180/statistics_qry_wo_sameProteins.txt
    }
    else {
        $FullStatFL =
          "$jobDIR/statistics.txt";    #../data/benchmark180/statistics_qry.txt
    }

    #-------------------------
    print LOG
"\nNote: The full statistics FL is $FullStatFL (contains all possible qry-homolog pairs).\n\n";
    my $startTime  = localtime();
    my $startTime2 = time();
    print LOG "\n\n**Interface prediction starts at $startTime.**\n\n";

    #-------------------------
    # ----- collect all the possible homologs' info into statistics.txt
    my $maxEvalue    = 10;
    my $minPositiveS = 0;

    ( $flag_prediction, my $QRYpairs_ref ) = &conservAnalysis_PSHomPPI(
        $jobDIR,    $qryPairFL,      $protSeqFL, $seqInt_protFL,
        $maxEvalue, $minPositiveS,   $intDef,    $atomDistThr,
        $rasaThr,   $rmSameProtFlag, $delLstFL
    );

    print LOG "\n\n\nCONSERVATION ANALYSIS DONE. \n\n\n";


    if ( $flag_prediction == 1 ) {
        print LOG
"\nNo homo-interologs can be identified for the query protein pairs.\n";

        &writePredictionFL_noHomologFound($jobDIR);
    }

    #------------------

    my $flag_finalPrediction = $flag_prediction;
    if ( $flag_prediction == 0 ) {

        if ( !-e $FullStatFL ) {
            die("$FullStatFL  does not exist. Check:$!");
        }

        #--------
        #predict interface residues

        &safeTwiDark( $jobDIR, $jobID, $qryPairFL, $FullStatFL );

#--------
# evaluate templates in the $FullStatFL file, and label each with safe, twilight, dark zones
# &evaluateTemplates($FullStatFL);
#--------
#predict interacting residue PAIRS based on clustered templates
        $flag_finalPrediction = &cluster_PS_resiPair( $jobDIR, $qryPairFL )
          ; #0: at least 1 qry pair has Ca-Ca prediction; 1: No Ca-Ca predictions made for any qry pair

    }

    #------
    my $endTime  = localtime();
    my $endTime2 = time();
    my $timeUsed = $endTime2 - $startTime2;

#---- collect CA-CA files for all queries into one folder: $jobDIR/Ca_Ca_distances
    my @QRYpairs = @{ &readPairFL($qryPairFL) };    #(1a2kA:1a2kC, 1a14H:1a14N)
    my $outputCACA_DIR = "$jobDIR/Ca_Ca_distances";
    if ( $flag_finalPrediction == 0 ) {

#---- collect all the prediction results (binary prediction and pred_scores) for the queries into one file
        my $outputFL1 =
          "$jobDIR/InterfacePredictions.final.txt";    #output file for usr
        my $outputFL2 = "$jobDIR/InterfacePredictions.internal.txt"
          ; #output file for machine reading, e.g., distOf_pInt_dockInt() in DockRank
        &collectPredictionFLsIntoOneFL1(
            $jobDIR,    \@QRYpairs, $outputFL1,
            $startTime, $endTime,   $timeUsed
        );    #write usr-friendly format of final predictions
        &collectPredictionFLsIntoOneFL2(
            $jobDIR,    \@QRYpairs, $outputFL2,
            $startTime, $endTime,   $timeUsed
        );    #write machine-friendly format of final predictions
        &collectPredictedResiPairFLsIntoOneFolder( $jobDIR, \@QRYpairs,
            $outputCACA_DIR, $startTime, $endTime, $timeUsed )
          ;    # collect all CA-CA files into one folder

# ----- collect all the prediction results (TP,TN,CC...) for the queries into one file
#    my $finalStatFL = "$jobDIR/InterfacePredictions.stat";
#    &collectStat_KNN_prot( $outputFL2, $finalStatFL, $scoreThr );
    }

    #-------
    print LOG
"\n\n\nNote: The full statistics FL is $FullStatFL (contains all possible qry-homolog pairs).\n\n";

    print LOG
"\n\n**Interface prediction is finished at $endTime. $timeUsed seconds used. **\n\n\n";

    return $flag_finalPrediction;

}

sub conservAnalysis_PSHomPPI {

    #author: Xue, Li
    #date: Jul 13th,2012
    #

#Input files:
#
#1. Protein ID pair list ($protPairFL)
#	prot1:prot2
#	prot3:prot4
#2. seq file of query proteins ($protSeqFL)
#	>prot1
#	seq_prot1
#	>prot2
#	seq_prot2
#3. a list of homo-interologs to be removed ($delLstFL, generate by PS-HomDDI/GenDeleteFL.pl)
#	1a2kA:1a2kC =>
#	1a14H:1a14N => 1NMB_H:1NMB_N,1A14_H:1A14_N,1NMA_H:1NMA_N
#4. seq_int file of query prot pairs (optional)
#   >prot1|prot1:prot2
#   seq_prot1
#   interface of prot1 with prot2
#
#

#
#Output:  "..\data\statistics.txt"
#output format:
#>query1 pdbID chainID:
#PDBID+CHAINID Bit_score Positive_Score IdentityScore aligLen_Query aligLen_Homolog CC TP TN FP FN
#...
#>query2
#...
#

# 3. pdbIDchnIDmap
#
# Format:
#	Input protein ID, pdbID_chainID
#
#	A, 1avxA
#	B, 1avxB
#	C, 1qquA
#
#
#perl ./conservAnalysis.pl $protPairFL, $protSeqFL maxEval minPositiveSThr  atomDistThr RASA_thr $rmSameProtFlag, $delLstFL

# 		perl ./conservAnalysis.pl  ../data/trans208/trans208ChainPair_dimer.lst 10 0  atomDistances 4 5
# 		perl ./conservAnalysis.pl  ../data/oblig115/oblig115ChainPair_dimer.lst 10 0  atomDistances 4 5

    #		( $flag_prediction, my $QRYpairs_ref ) = &conservAnalysis_PSHomPPI(
    #			$protPairFL, $protSeqFL, $seqInt_protFL,
    #			$maxEvalue,      $minPositiveS,
    #			$intDef,  $atomDistThr,    $rasaThr,
    #			$rmSameProtFlag, $delLstFL
    #		);

    use strict;
    use File::Basename;
    use File::Copy;
    use File::Path;

    #-------------------
    #global variables
    our $homeDIR;
    our $searchNewPPIDBDIR;    #= "$NPSHomPPIDIR/searchNewPPIDB";
    our $pdbIDchnIDmap
      ; #a map between chain IDs in $inputFL and pdb ID chain ID. It is used to remove homo-interologs that are highly similar to query pairs
    our $intNum_cutoff;    # 4;
    our $PvalThr_hhpredPDB70;

    #------------
    #input variables

    my $jobDIR           = shift @_;
    my $protPairFL       = shift @_;
    my $protSeqFL        = shift @_;
    my $seqInt_protFL    = shift @_;
    my $EvalThr          = shift @_;  #1e-10 , 10
    my $positiveScoreThr = shift @_;  #[0-1] #-- not used for the hhpred version
    my $intDef           = shift @_;  #'vanDerWaalsDistances';
    my $atomDistThr      = shift @_;  #8.5
    my $rasaThr          = shift @_;  #0
    my $rmSameProtFlag   = shift @_;  #'Yes'  or 'No'
    my $delLstFL         = shift @_;
    my $flag             = 1
      ; #$flag=0: at least one of the protein pairs have homo-interologs in $FullStatFL
        #$flag=1: none of the query pairs have homo-interologs in $FullStatFL

    my $jobID    = basename( $protPairFL, ( '.fasta', '.lst', '.txt' ) );
    my $statFL   = "$jobDIR/statistics\_$jobID.txt";
    my $outputFL = $statFL;

    #files
    my $basename = basename( $protPairFL, '.lst' );
    my $statisticsFL = "$jobDIR/statistics.txt";
    my $QryWoHomComplexLST =
      "$jobDIR/qryWoTemplate.lst"
      ;    #keep a list of Qry pairs that cannot find homo-complexes

    #--------------------------------------------------------------#
    #program begins
    my @qryWoTemplate;

    #generate folders for each query pair
    if ( !-e $seqInt_protFL ) {

        #use '?' as interface
        &genNullSeqIntFL( $protPairFL, $protSeqFL, $seqInt_protFL );
    }
    my @QRYpairs =
      @{ &foldersGen1( $jobDIR, $protPairFL, $protSeqFL, $seqInt_protFL ) };

    #conservation analysis on each test chain pair
    foreach my $QRYpair (@QRYpairs) {

        #	A:B
        print LOG "\n\n----------------------------------\n\n";
        print LOG "\tNow conservation analysis on $QRYpair ";    #'A:B'
        print LOG "\n\n----------------------------------\n\n";

        my ( $chain1, $chain2 ) = split( /:/, $QRYpair );

        my $query1 = $chain1;
        my $query2 = $chain2;

        #DIRs

        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $jobDIR );
        my $qryDIR1 = "$pairDIR/$query1";
        my $qryDIR2 = "$pairDIR/$query2";

        #files

        my $Qry1FstFL =
          "$qryDIR1/$chain1.fasta.txt"
          ;    #"$pairDIR/$query1/$pdbID\_$chainID1.fasta.txt";
        my $Qry2FstFL = "$qryDIR2/$chain2.fasta.txt";
        my $Qry1hhpredFL = "$qryDIR1/$chain1.hhpred";  #-- output of Mick's code
        my $Qry2hhpredFL = "$qryDIR2/$chain2.hhpred";  #-- output of Mick's code
        my $homologLisFL1 =
          "$qryDIR1/homologsOfquerySeq.lst";  #the output file of parse_blast.pl
        my $homologLisFL2 =
          "$qryDIR2/homologsOfquerySeq.lst";  #the output file of parse_blast.pl
        my $hom_complexesLst       = "$pairDIR/templates.lst";
        my $seq_int_homComplexesFL = "$pairDIR/seq_int_templates.lst";

        #call Blast and get homolog lists for two queries
        if ( !-e $Qry1FstFL ) {
            print LOG "$Qry1FstFL does not exist. Next pair.\n";
            next;
        }
        if ( !-e $Qry2FstFL ) {
            print LOG "$Qry2FstFL does not exist. Next pair.\n";
            next;
        }

        if ( !-e $Qry1hhpredFL ) {

      #&callpsiBlast( $BlastDataset, $Qry1FstFL, $Qry1hhpredFL, $EvalThr )
      #  ;    #generate PSSM first and then Blast the query against PDB database

            #&callpsiBlast_old( $BlastDataset, $Qry1FstFL, $Qry1hhpredFL,
            #    $EvalThr );    #Blast the query directly against PDB database

            my @param_hhpred = ( $Qry1hhpredFL, $EvalThr );
            &call_hhpred( $Qry1FstFL, \@param_hhpred );    #-- call Mick's code
        }
        else {
            print LOG "$Qry1hhpredFL pre-exist. No need to run callhhpred.\n";
        }

        if ( !-e $Qry2hhpredFL ) {

      #&callpsiBlast( $BlastDataset, $Qry2FstFL, $Qry2hhpredFL, $EvalThr )
      #  ;    #generate PSSM first and then Blast the query against PDB database

            #&callpsiBlast_old( $BlastDataset, $Qry2FstFL, $Qry2hhpredFL,
            #    $EvalThr );    #Blast the query against PDB database
            my @param_hhpred = ( $Qry2hhpredFL, $EvalThr );
            &call_hhpred( $Qry2FstFL, \@param_hhpred );    #-- call Mick's code
        }
        else {

            print LOG "$Qry2hhpredFL pre-exist. No need to run callhhpred.\n";
        }

        #find hom-complexes from the homolog lists
        my $hhpredProbThr = $PvalThr_hhpredPDB70;
        &parse_hhpred_userProtDB( $Qry1hhpredFL, $hhpredProbThr )
          ; #parse the homologs' PDBID from the resulting blast file, and write them into a list as the input of PPIDB
        &parse_hhpred_userProtDB( $Qry2hhpredFL, $hhpredProbThr )
          ;    #parse the homologs' PDBID from the resulting blast

#--note: remove bound complexes from homolog list. This step is for DockRank paper.
        my $caseDIR = dirname($jobDIR);
        my $boundPDBIDsFL =
          "$caseDIR/complexes_to_be_del_from_homcomplexes.lst";
        if ( -e $boundPDBIDsFL ) {
            print LOG "Remove bound PDB ID from homolog list.\n";

            my @boundPDBIDs = @{ &readLstFL($boundPDBIDsFL) };

            foreach my $boundPDBID (@boundPDBIDs) {
                &rmBoundFromHomologFL( $homologLisFL1, $boundPDBID );
                &rmBoundFromHomologFL( $homologLisFL2, $boundPDBID );
            }
        }
        else {
            print LOG
"\nwarning: bound PDB ID file is not available. Bound pdb ID is NOT removed from homologs!!!\n\n";
        }

#compare two homolog lists
#if they have the same pdbID and different chainID and they interact with at least 5 residues
#then this pdb complex is used to infer the interfaces of query1 and query2

        &compareHomlogLsts( $homologLisFL1, $homologLisFL2 )
          ;    #outputFL: ../data/templates.lst

        if ( !-e $hom_complexesLst ) {
            print LOG "\n$QRYpair\tcannot find homo-complexes.\n\n";
            push @qryWoTemplate, $QRYpair;
            next;
        }
        else {

            $flag = 0
              ; #$flag=0: at least one of the protein pairs have homo-interologs in $FullStatFL
                #$flag=1: none of the query pairs have homo-interologs in $FullStatFL

            #search new PPIDB for the interfaces of the hom-complexes of query1:query2
            &callPPIDB_PS( $hom_complexesLst, $intDef, $atomDistThr, $rasaThr );


            #rename the output of callPPIDB_PS to a better name $seq_int_homComplexesFL
            $basename = basename( $hom_complexesLst, ( '.lst', '.txt' ) );
            unlink $seq_int_homComplexesFL if ( -e $seq_int_homComplexesFL );
            move( "$pairDIR/seq_int_qryComplexes_$basename.lst",
                $seq_int_homComplexesFL );

            #--------------
            #remove homologs with interfaces less than $intNum_cutoff residues from seq_int_qryComplexes_templates.lst
            &rmChainsWithLessThan5intAA_PS( $intNum_cutoff,
                $seq_int_homComplexesFL );

            #generate local alignment of interfaces file

            &localAlignmentGen_PS( $pairDIR, $query1, $seq_int_homComplexesFL );
            &localAlignmentGen_PS( $pairDIR, $query2, $seq_int_homComplexesFL );

        }

    }

    if (@qryWoTemplate) {
        &writeArrayIntoFL( \@qryWoTemplate, $QryWoHomComplexLST );
    }

    if (scalar @QRYpairs == scalar @qryWoTemplate){
    #-- none of the queries have dimer templates
        print LOG "\nUnfortunately, none of the query pairs have dimer templates.\n";
        return ( $flag, \@QRYpairs );
    }

    #collect the alignment statistics for each homo-interolog

    my $ori_statFL = "$jobDIR/statistics.txt";
    &collectStat( \@QRYpairs, $jobDIR, $ori_statFL );

    #- predict CC for each template and add the pred_CC to $ori_statFL
    &addPredCC($ori_statFL);

    #- exclude user-specified templates
    my $final_statFL;

    if ( $flag == 0 && $rmSameProtFlag eq 'Yes' ) {
        $flag = &rmSameProtPair_deleteFL( $ori_statFL, $delLstFL );
        my $dirname = dirname($ori_statFL);
        $basename = basename( $ori_statFL, ( '.txt', '.lst' ) );
        $final_statFL = "$dirname/$basename\_wo_sameProt.txt";    #
    }
    else {
        $final_statFL = $ori_statFL;
    }

    #
    print LOG "output1: $QryWoHomComplexLST (Qry without hom-interologs)\n";
    print LOG "output2: $final_statFL\n";

    return ( $flag, \@QRYpairs );

 #$flag=0: at least one of the protein pairs have homo-interologs in $FullStatFL
 #$flag=1: none of the query pairs have homo-interologs in $FullStatFL
}

sub genNullSeqIntFL {

#
# Given query protein pairs and sequence files, generate the null PS-seq_int file, which means interfaces are '?'.
#
# Input:
#1. Protein ID pair list
#	prot1:prot2
#	prot3:prot4
#2. seq file of query proteins
#	>prot1
#	seq_prot1
#	>prot2
#	seq_prot2
# Output:
# seq_int file of query prot pairs
#   >prot1|prot1:prot2
#   seq_prot1
#   int_prot1 with prot2 (all question marks)

    my ( $protPairFL, $protSeqFL, $seqInt_protFL ) = @_;

    my @protPairs = @{ &readPairFL($protPairFL) };   #(1a2kA:1a2kC, 1a14H:1a14N)
    my %protSeqs  = %{ &readProtSeqFL($protSeqFL) };
    my @prots     = keys %protSeqs;

    print LOG
      "\n\nUser did not provide the interface file of their query proteins.\n";
    print LOG
"Now  generate Null seq_int file for query proteins with question marks denoting the interfaces.\n";

    unlink $seqInt_protFL if ($seqInt_protFL);
    open( OUTPUT, ">>$seqInt_protFL" ) || die("Cannot open $seqInt_protFL:$!");
    print OUTPUT
      "# Null seq_int file generated by genNullSeqIntFL() in PSHomPPI.pm.\n";

    foreach my $pair (@protPairs) {
        my ( $protA, $protB ) = split( /:/, $pair );

        if ( !defined $protSeqs{$protA} ) {
            die("No seq of $protA is read from  $protSeqFL:$!");
        }
        if ( !defined $protSeqs{$protB} ) {
            die("No seq of $protB is read from  $protSeqFL:$!");
        }

        #---
        my $seqA  = $protSeqs{$protA};
        my $seqB  = $protSeqs{$protB};
        my $int_A = '?' x length($seqA);
        my $int_B = '?' x length($seqB);

        print OUTPUT ">$protA|$protA:$protB\n";
        print OUTPUT "$seqA\n";
        print OUTPUT "$int_A\n";

        print OUTPUT ">$protB|$protB:$protA\n";
        print OUTPUT "$seqB\n";
        print OUTPUT "$int_B\n";

    }
    close OUTPUT;

    print LOG "The null seqInt file $seqInt_protFL is generated.\n";

}

sub foldersGen1 {

    #author: Xue, Li
    #date: 7/17/2012

#This script is part of the pipeline of partner-spcific HomoPPI. It is different from folderGen.pl for HOMPPI

    #Generate folders from seq_int_file.

# 1.  generate a folder for each interacting pair in $protPairFL
#		under which there are two subfolders. For example,
#		1a2kA_1a2kC/1a2kA
#		1a2kA_1a2kC/1a2kC
# 2.  generate fasta file, eg, 1a2kA.fasta.txt, using $seqFL_Qry file.
# 3.  generate $taskID.int file using $seqIntFL_Qry file. If $taskID does not exist in $seqIntFL_Qry, its int are '?'s.
#
#
#Input seq_int file example:
#
#	>1a2kA|1a2kA:1a2kC
#	KPIWEQIGSSFIQHYYQLFDNDRTQLGAIYIDASCLTWEGQQFQGKAAIVEKLSSLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKADEDPIMGFHQMFLLKNINDAWVCTNDMFRLALHNFG
#	0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#	>1a2kC|1a2kC:1a2kA
#	QVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVATLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFVAMPALAPPEVVMDPALAAQYEHDLEV
#	0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

#
#note: remove other folders before running this script.
#
#perl ./foldersGen1.pl  ../data/seq_int_qryComplexes_test.lst
#perl ./foldersGen1.pl  ../data/seq_int_qryComplexes_perm147pairedChains.lst
#perl ./foldersGen1.pl  ../data/trans212/seq_int_qryComplexes_trans212pairedChains.lst

    use strict;
    use Cwd;
    use File::Basename;
    use File::Path;
    use File::Copy;

    my $outputDIR    = shift @_;                 # = $jobDIR
    my $protPairFL   = shift @_;
    my $seqFL_Qry    = shift @_;
    my $seqIntFL_Qry = shift @_;
    my $jobDIR       = dirname($seqIntFL_Qry);
    my $taskID;               #1t9gS
    my $pdbID;
    my $chain1;
    my $chain2;
    my $homologPartnerSpc;    #>1t9gS|S:R
    my $interfaceFL;
    my $fastaFL;
    my $listFL;               #for the use of download pdb file
    my $seq;
    my $flag = 0;
    my $pairDIR;

    my %QRYpairs_new
      ; #it is the intersection between qry pairs listed in $qrypairFL and $seqIntFL_Qry
    my @qrys_not_in_seqIntFL;
    my @duplicatedQrypairs = ();

    my @QRYpairs_all =
      @{ &readPairFL($protPairFL) };    #(1a2kA:1a2kC, 1a14H:1a14N)
    my $num_qrypairs_inPairFL = scalar @QRYpairs_all;
    @QRYpairs_all = @{ &uniquePair( \@QRYpairs_all ) }
      ;    #-- remove redundant pairs. P:Q and Q:P are redundant

    my $num_pair = 0;

    my $seqs = &readProtSeqFL($seqFL_Qry);
    my ( $seqs_tmp, $ints ) =
      &readQrySeqIntFL($seqIntFL_Qry)
      ;    #$seqs->{2b42A}='SETETENON';  $ints->{2b42A|2b42A:2b42B}='011000000'

    foreach my $pair (@QRYpairs_all) {

        #'chnA:chnB'
        #'1fevA:1fevB'

        ( $pairDIR, my $prot1, my $prot2 ) = &getPairDIR( $pair, $outputDIR );

        #- $pairDIR = '/home/lixue/PSHOMPPIv1.3/uploadData/test/P:Q'

        #-------- creat $pairDIR folders

        print LOG "\n\nGenerating a folder: $pairDIR\n\n";

        mkdir($pairDIR) if ( !-d $pairDIR );

        $num_pair++;
        my $pair_new = basename($pairDIR)
          ;    #-- the protein pair is alphabatically reordered: B:A => A:B
        $QRYpairs_new{$pair_new} = 1;

        my @prots = ( $prot1, $prot2 );
        for ( my $i = 0 ; $i < 2 ; $i++ ) {

            my $prot = $prots[$i];

            #-------- creat $pairDIR/$prot1 subfolders
            print LOG "Generating a folder: $pairDIR/$prot\n";
            if ( !-d "$pairDIR/$prot" ) {
                mkdir("$pairDIR/$prot")
                  || die("Cannot make folder: $pairDIR/$prot: $!");
            }

            #write $prot.fasta.txt
            $fastaFL = "$pairDIR/$prot/$prot.fasta.txt";
            my $seq = $seqs->{$prot};
            &writeFastaFL( $prot, $seq, $fastaFL );

            #---write int file
            if ( !defined $ints->{"$prot1|$prot1:$prot2"} ) {

                print LOG
"\nWarning: The interfaces of $prot1 interacting with $prot2 are not defined in $seqIntFL_Qry. Its int is assigned as question marks.\n\n";
                push @qrys_not_in_seqIntFL, $pair;

                $ints->{"$prot1|$prot1:$prot2"} =
                  '?' x length( $seqs->{$prot1} );
            }

            if ( !defined $ints->{"$prot2|$prot2:$prot1"} ) {

                print LOG
"\nWarning: The interfaces of $prot2 interacting with $prot1 are not defined in $seqIntFL_Qry. Its int is assigned as question marks.\n\n";
                push @qrys_not_in_seqIntFL, $pair;

                $ints->{"$prot2|$prot2:$prot1"} =
                  '?' x length( $seqs->{$prot2} );
            }

            #write $prot.int file (interface file)
            $interfaceFL = "$pairDIR/$prot/$prot.int";

            my $prot_prtner;
            if ( $i == 0 ) {

                $prot_prtner = $prots[1];
            }
            else {
                $prot_prtner = $prots[0];
            }
            my $int = $ints->{"$prot|$prot:$prot_prtner"};
            &writeIntFL( $interfaceFL, $prot, $seq, $int );
        }

    }

    my @QRYpairs_final     = keys %QRYpairs_new;
    my $num_qrypairs_final = scalar @QRYpairs_final;

    print LOG
"\n$num_pair interacting-pair folders are generated! In each interacting-pair folder, there are two subfolders for the binding partners. In each subfolder, there is a fasta file and an interface file.\n\n";

    print LOG "$num_qrypairs_inPairFL qry pairs listed in $protPairFL. There are $num_qrypairs_final unique query pairs.\n";

    return \@QRYpairs_final;

}

sub writeIntFL {
    my $interfaceFL = shift @_;
    my $header      = shift @_;
    my $seq         = shift @_;
    my $int         = shift @_;
    unlink $interfaceFL if ( -e $interfaceFL );
    open( INTFL, ">> $interfaceFL" )
      || die("Cannot creat $interfaceFL:$!");
    print INTFL ">$header\n";
    print INTFL "$seq\n";
    print INTFL "$int\n";
    close(INTFL);

    print LOG "$interfaceFL generated.\n";
}

sub readQrySeqIntFL {

#Input example:
#		>1avxA|1avxA:1avxB
#		IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN
#		0000000000000000000000011000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000001000000000000000000110000000000000000000000000000000

    my $SeqIntFL = shift @_;
    my $header;
    my $qryProtID;
    my %seqs;
    my %int;
    my $length_seq;
    my $length_int;

    open( seq_int_FL, "<$SeqIntFL" )
      || die("Cannot open $SeqIntFL:$!");
    while (<seq_int_FL>) {
        s/[\n\r]//mg;

        if (/^>(([\-\.\w]+)\|[\-\.\w]+:[\-\.\w]+)/) {

            #>1avxA|1avxA:1avxB

            $header    = $1;    #1avxA|1avxA:1avxB
            $qryProtID = $2;    #1avxA
            next;
        }

        if (/^[A-Za-z]+$/) {
            my $seq = $_;

            if (! defined $qryProtID ){
                die ("qryProtID not defined in $SeqIntFL. Check the format:$!");
            }

            $seqs{$qryProtID} = $seq;
            $length_seq       = length($seq);
            $qryProtID        = '';
            next;
        }
        if (/^[01\?\-]+$/) {
            my $int = $_;
            $int{$header} = $int;
            $length_int = length($int);

            if ( $length_seq ne $length_int ) {
                die(
"seq length and int length are not equal. Check $SeqIntFL:$!"
                );
            }
            $header     = '';
            $length_seq = 0;
            $length_int = 0;
            next;
        }

    }
    close(seq_int_FL);

    if ( !%seqs || !%int ) {
        die("Check $SeqIntFL:$!");
    }

    return ( \%seqs, \%int );

}

#---------------------------
sub getPairDIR2 {
    my $pair    = shift @_;    #'A:B'
    my $dataDIR = shift @_;
    my $pairDIR;

    my ( $chain1, $chain2 ) = split( /:/, $pair );

    if ( !defined $chain1 || !defined $chain2 ) {

        die("Check $pair:$!");
    }
    my $pairDIR1 = "$dataDIR/$chain1:$chain2";    #possible pairDIR
    my $pairDIR2 = "$dataDIR/$chain2:$chain1";    #another possible pairDIR

    if ( $chain1 lt $chain2 ) {
        $pairDIR = $pairDIR1;
    }
    else {
        $pairDIR = $pairDIR2;
    }

    return ( $pairDIR, $chain1, $chain2 );

}

sub call_hhpred {
    use strict;
    use File::Spec;

    our $callHHpredPY;
    our $logFL;
    our $PYTHON ; # '/usr/bin/python';
    our $hhpred_param; #-- hash ref

    my $qryFstFL = shift @_;    #P.fasta.txt
    my $outputFL = shift @_;
    my $hhpredParam = shift @_; # a hash ref

    #--- prepare hhpred_raw_result folder
    #my $dataDIR       = File::Spec->rel2abs( dirname($qryFstFL) );
    my $dataDIR = dirname($qryFstFL);
    my $raw_resultDIR = "$dataDIR/hhpred_raw_result";
    mkdir $raw_resultDIR if ( !-d $raw_resultDIR );
    copy( $qryFstFL, "$raw_resultDIR/query.fa" );

    #-- write settings.config file
    my $hhpred_configFL = "$raw_resultDIR/settings.config";
    &write_hhpredConfigFL($hhpred_configFL, $hhpred_param);

    #-- call hhpred

    my $logFL_tmp = "$logFL.hhpredtmp";
    #my @command =
    #  ( $PYTHON, $callHHpredPY, $raw_resultDIR );

    my $command =
    " $PYTHON $callHHpredPY $raw_resultDIR > $logFL_tmp ";
    print LOG "CALL HHPRED COMMAND: $command\n";

    system($command) == 0 or die("ERROR: $command failed:$?"); # there is a bug in addss.pl (see TODO)

    #-- copy final result files to $dataDIR
    my $qryID = basename( $qryFstFL, ( '.fa', '.fasta', '.fasta.txt' ) );
    copy( "$raw_resultDIR/query.fa.outputCOMPUTER", "$dataDIR/$qryID.hhpred" );

    #--- attach $logFL.tmp to $logFL
#    system("cat $logFL_tmp >> $logFL ") ==0 or die ("Cannot attach $logFL_tmp to $logFL:$!");
#    unlink "$logFL_tmp" if (-e "$logFL_tmp");

}


sub call_hhpred_old {
    use strict;
    use File::Spec;

    our $callHHpredPY;
    our $logFL;
    our $PYTHON;    # '/usr/bin/python';
    our $hhpred_param;
    our $queryName;    # 'query.fa'

    my $qryFstFL    = shift @_;    #P.fasta.txt
    my $outputFL    = shift @_;
    my $hhpredParam = shift @_;    # a hash ref

    #--- prepare hhpred_raw_result folder
    #    my $dataDIR       = File::Spec->rel2abs( dirname($qryFstFL) );
    my $dataDIR       = dirname($qryFstFL);
    my $raw_resultDIR = "$dataDIR/hhpred_raw_result";
    mkdir $raw_resultDIR if ( !-d $raw_resultDIR );
    copy( $qryFstFL, "$raw_resultDIR/$queryName" );

    #-- write settings.config file
    my $hhpred_configFL = "$raw_resultDIR/settings.config";
    &write_hhpredConfigFL( $hhpred_configFL, $hhpred_param );

    #-- call hhpred
    my @command = ( $PYTHON, $callHHpredPY, $raw_resultDIR, " >> $logFL" );

    print LOG "CALL HHPRED COMMAND: @command\n";
    system(@command) == 0 or die("ERROR: @command failed:$!");

    #-- copy final result files to $dataDIR
    my $qryID = basename( $qryFstFL, ( '.fa', '.fasta', '.fasta.txt' ) );
    copy( "$raw_resultDIR/$queryName.outputCOMPUTER",
        "$dataDIR/$qryID.hhpred" );

}

#--------------------------

sub parse_hhpred_userProtDB {

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
    my $inputfile     = shift @_;    # Mick's hhpred output file: Q.hhpred
    my $hhpredProbThr = shift @_;
    my $taskID = basename( dirname($inputfile) );
    my $inputfilename = basename( $inputfile, '.hhpred' );
    my $PDBID;
    my $chainID;
    my $roundNum;                    #the number of iteration of PSI-BLAST
    my $homologRef;
    my $homologID;
    my @homlog
      ; #to store the PDBID and chainID of homologs of query sequence returned Blastp
    my $similarity;
    my $num_seq = 0;    #number of the homologs for each query sequence

    #directories

    print LOG
"\n\nParse hhpred output file ($inputfile) with the following thresholds:\n";
    print LOG "\tHHpred prob. Thr: $hhpredProbThr\n\n\n";

    my $outputDIR   = dirname($inputfile);
    my $outfilename = "homologsOfquerySeq.lst";    #define output file name
    my $outfile     = "$outputDIR/$outfilename";
    unlink $outfile if ( -e $outfile );

    #------------- program begins --------------------

#read the blast output file. Only last round of PSI-blast results are returned in $blastResults.
    my ( $qryID, $qryLength, $blastResults ) =
      &readHHpredOUTPUTFL_userProtDB($inputfile);

    if ( !defined $qryID && !defined $qryLength && !defined $blastResults ) {
        print LOG "No hit found in $inputfile\n";
        @$homologRef = ();
        &writeHomologLstfile( $outfile, $homologRef, $hhpredProbThr );
        print LOG "parse_psiblast_userDB.pl finish!  Outputfile: $outfile\n\n";
        return;
    }

    #--get homologs that have PositiveS and IdentityS >= thresholds
    #print LOG "\nPassed homologs:\n";
    #    my @homologs = keys(%$blastResults);

    foreach my $homolog ( keys %{$blastResults} ) {

        foreach my $num ( keys %{ $blastResults->{$homolog} } ) {
            my $hhpredProb =
              $blastResults->{$homolog}->{$num}->{'hhpredProb'};    #[0-100]

            if ( !defined $hhpredProb ) {
                die(
"homolog $homolog ($num) does not have hhpred prob. defined in $inputfile:$!"
                );
            }

            if ( $hhpredProb >= $hhpredProbThr ) {
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
    &writeHomologLstfile( $outfile, \@homologs_final, $hhpredProbThr );

    print LOG "parse_hhpred_userProtDB{} finish!  Outputfile: $outfile\n\n";
}

sub readLstFL {

    #Input:
    #4b2cA
    #3rdzB

    my $FL = shift @_;
    my @elements;

    open( INPUT, "<$FL" ) || die("Cannot open $FL:$!");

    foreach (<INPUT>) {
        s/[\n\r]//mg;
        if (/^([\w\.:]+)[\s\t]{0,}$/) {
            push @elements, $1;
        }
    }
    close INPUT;

    if ( !@elements ) {
        print LOG "$FL is empty.\n";
        @elements = ();
    }

    my $num = scalar @elements;

    #print "$num elements read from $FL.\n";

    return \@elements;

}

#---------
sub rmBoundFromHomologFL {

    #Input FL:
    #1ql4C
    #3nwvC
    #3nbtA
    #2b11B

    my $homologFL  = shift @_;
    my $BoundPdbID = shift @_;

    print LOG "\nRemove Bound pdb ID $BoundPdbID from $homologFL.\n";

    my $r        = rand();
    my $dirname  = dirname($homologFL);
    my $tempFL   = "$dirname/homolog.$r.tmp";
    my @homologs = @{ &readLstFL($homologFL) };

    unlink $tempFL if ( -e $tempFL );
    open( OUTPUT, ">>$tempFL" ) || die("Cannot open $tempFL:$!");
    foreach my $homolog (@homologs) {
        if ( $homolog =~ /^$BoundPdbID/i ) {
            print LOG "$homolog of Bound complex is removed from $homologFL\n";
            next;
        }
        print OUTPUT "$homolog\n";

    }
    close OUTPUT;

    unlink $homologFL;
    move( $tempFL, $homologFL )
      || die("Cannot rename $tempFL to $homologFL:$!");

}

#---------------------
sub compareHomlogLsts {

    #!usr/bin/perl -w
    #Li Xue
    #3/27/2010
    #
    #slow!!

#compare two homolog lists
#if they have the same pdbID and different chainID and they interact with at least 5 residues
#then this pdb complex is used to infer the interfaces of query1 and query2
# perl ./compareHomlogLsts.pl HomologLst1 HomologLst2
#
#INPUT: two lists of homologs of two binding partners
#    1d5nB
#    1tc6A
#    1vewC
#
#OUTPUT:
#1. a list of homo-complexes "../data/templates.lst"
#2. a list of homologs that are used to form hom-complexes for query chain 1
#3. a list of homologs that are used to form hom-complexes for query chain 2
#
#note: if one or two binding partners cannot find homologs, no output file will be generated.
#
#Usage example:
# perl ./compareHomlogLsts.pl ../data/1auiA_A_B/homologsOfquerySeq.lst ../data/1auiB_B_A/homologsOfquerySeq.lst
# perl ./compareHomlogLsts.pl ../data/1dkg_A_D/1dkgD_D_A/homologsOfquerySeq.lst ../data/1dkg_A_D/1dkgA_A_D/homologsOfquerySeq.lst
# perl ./compareHomlogLsts.pl ../data/trans208/2prg/2prg_B_C/2prgB_B_C/homologsOfquerySeq.lst ../data/trans208/2prg/2prg_B_C/2prgC_C_B/homologsOfquerySeq.lst
    use strict;
    use File::Basename;

    my $inputFL1 = shift @_;
    my $inputFL2 = shift @_;

    my $query1    = basename( dirname($inputFL1) );
    my $query2    = basename( dirname($inputFL2) );
    my ($pairDIR) = $inputFL1 =~ /(.+)\/$query1\/homologsOfquerySeq.lst/;

    my @Qry1homologs_in_homComplexes
      ;    #homologs of Qry chain 1 that forms homComplexes
    my @Qry2homologs_in_homComplexes;
    my @homComplexes;
    my $homComplexesFL = "$pairDIR/templates.lst";
    my $homologs_in_homComplexesFL1 =
      "$pairDIR/$query1/homologs_in_homComplexes.lst";
    my $homologs_in_homComplexesFL2 =
      "$pairDIR/$query2/homologs_in_homComplexes.lst";

    #----------------

    print LOG
"\n\nCompare two homolog lists ($inputFL1 and $inputFL2) and write homo-interolog list....\n\n";

    my $homologLst1 =
      &readPDBHomologLstFL($inputFL1);    #hash ref: $list1->{$pdbID} = @chainID
    my $homologLst2 = &readPDBHomologLstFL($inputFL2);    #%hash ref

    unlink($homComplexesFL) if ( -e $homComplexesFL );
    unlink($homologs_in_homComplexesFL1)
      if ( -e $homologs_in_homComplexesFL1 );
    unlink($homologs_in_homComplexesFL2)
      if ( -e $homologs_in_homComplexesFL2 );

    &headerHomcomplexFL($homComplexesFL);

    if ( defined $homologLst1 && defined $homologLst2 ) {

        #if both binding partners have homologs
        my @pdbIDlst1 = keys(%$homologLst1);

        #        print LOG "\nPDB IDs in $inputFL1: @pdbIDlst1\n";

        my @pdbIDlst2 = keys(%$homologLst2);

        #        print LOG "\nPDB IDs in $inputFL2: @pdbIDlst2\n";

        my %seen;
        my @commonPDBIDs;
        if ( scalar @pdbIDlst1 >= scalar @pdbIDlst2 ) {
            @seen{@pdbIDlst1} = (1) x scalar @pdbIDlst1;
            @commonPDBIDs = grep { $seen{$_} } @pdbIDlst2;
        }
        else {
            @seen{@pdbIDlst2} = (1) x scalar @pdbIDlst2;
            @commonPDBIDs = grep { $seen{$_} } @pdbIDlst1;

        }

        if (@commonPDBIDs) {
            print LOG "\nCommon PDB IDs are: @commonPDBIDs\n";

            foreach my $pdbid (<@commonPDBIDs>) {
                my @homologChainIDs_Qry1 =
                  &unique( @{ $homologLst1->{$pdbid} } );
                my @homologChainIDs_Qry2 =
                  &unique( @{ $homologLst2->{$pdbid} } );

                foreach my $chainid1 (@homologChainIDs_Qry1) {
                    foreach my $chainid2 (@homologChainIDs_Qry2) {
                        if ( $chainid1 ne $chainid2 ) {

                            push @homComplexes,
                              "$pdbid$chainid1:$pdbid$chainid2";

                            push @Qry1homologs_in_homComplexes,
                              "$pdbid$chainid1\|$chainid1:$chainid2";
                            push @Qry2homologs_in_homComplexes,
                              "$pdbid$chainid2\|$chainid2:$chainid1";

                        }
                    }
                }

            }

            @Qry1homologs_in_homComplexes =
              &unique(@Qry1homologs_in_homComplexes);
            @Qry2homologs_in_homComplexes =
              &unique(@Qry2homologs_in_homComplexes);

#-- do not use uniquePair(). uniquePair()will not keep the order of chain IDs in one pair
#@homComplexes = @{ &uniquePair( \@homComplexes ) };
            @homComplexes = &unique(@homComplexes);

            #2greB:2greO -> 2gre B:O
            @homComplexes = @{ &reformHomComplexPairs( \@homComplexes ) };

            &writeArrayIntoFL( \@Qry1homologs_in_homComplexes,
                $homologs_in_homComplexesFL1 );
            &writeArrayIntoFL( \@Qry2homologs_in_homComplexes,
                $homologs_in_homComplexesFL2 );
            &writeArrayIntoFL( \@homComplexes, $homComplexesFL );

            print LOG
"compareHomlogLsts.pl is finished. $homComplexesFL is generated. $homologs_in_homComplexesFL1 and $homologs_in_homComplexesFL2 are also generated. \n\n";

        }
        else {
            my $qryComplex = basename( dirname( dirname($inputFL1) ) );
            print LOG
"compareHomlogLsts.pl is finished. $qryComplex does not have homoComplexes.\n";
            unlink $homComplexesFL if ( -e $homComplexesFL );

        }

    }
    else {

        #at least one of the binding partners cannot find homologs
        if ( !$homologLst1 ) {

            #		../data/1i4fA_A_C/homologsOfquerySeq.lst
            my ($partner) = $inputFL1 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
            print LOG
"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";
            unlink $homComplexesFL if ( -e $homComplexesFL );

        }
        if ( !$homologLst2 ) {

            #		../data/1i4fC_C_A/homologsOfquerySeq.lst
            my ($partner) = $inputFL2 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
            print LOG
"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";
            unlink $homComplexesFL if ( -e $homComplexesFL );

        }
    }

    print LOG "\n";

}

sub callPPIDB_PS {

    #!/usr/bin/perl -w

#Author: Xue, Li
#Date: Nov, 2010
#
#This script is part of homolog based Protein-Protein interface prediction
#
#Check the PPIDB. If the proteins in the input file exists in PPIDB, output their sequence and interface information. Input file is the homologs of the query seq, which is generated by BLASTp.
#Input: the homologs of query sequence, e.g, ./data/1lukA/PDBIDchainID_homologsOfquerySeq.lst
#Output: the sequence and interface info of the homologs, e.g., seq_interface_homologs.lst
#
#Usage: perl ./callPPIDB_PS.pl PDBIDchainID.lst $intDef $atomDistThr $rasaThr
#e.g., perl ./callPPIDB_PS.pl ../data/3IKHA/homologsOfquerySeq.lst vanDerWaalsDistances 6 0

    use strict;
    use File::Basename;

    our $logFL;
    our $searchNewPPIDBDIR;
    our $searchNewPPIDB_PS_PL;

    my $inputfile   = shift @_;
    my $intDef      = shift @_;
    my $atomDistThr = shift @_;
    my $rasaThr     = shift @_;

    print LOG "\n\n\nSearching New PPIDB:\n";
    print LOG "\tInt Def: $intDef\n";
    print LOG "\tDist Thr: $atomDistThr\n";
    print LOG "\tRasa Thr: $rasaThr\n\n\n";

    #Directories
    use Cwd;
    my $codeDIR = cwd;

    #	my ($cgiDIR) = $codeDIR =~ /(.+)\/code$/;
    #	my $dirname = dirname($inputfile);
    #	($dirname) = $dirname =~ /^[\.]+\/(.+)/;
    #
    #	my $taskDIR = "$cgiDIR/$dirname";
    my $taskDIR = dirname($inputfile);

    #
    #my $taskID =basename(dirname($inputfile));
    #my $taskDIR="$dataDIR/$taskID";

   #
   #my $filename = basename($inputfile);
   #$inputfile="$dataDIR/$taskID/$filename";#get the absolute dir for $inputfile

#modify config file to designate that the output dir should be $dataDIR/$taskID
#To generate a tmp file, make the modification, and then copy the tmp file to the original file
#chdir("$PPIDBclient")||die("Cannot enter $PPIDBclient\n");

    my $r            = rand(1);
    my $basename     = basename( $inputfile, ( '.lst', '.txt' ) );
    my $new_configFL = "$searchNewPPIDBDIR/conf/$basename\_$r.config";
    my $old          = "$searchNewPPIDBDIR/conf/config_PS";

    open( oldFL, "<$old" ) || die("Cannot open $old:$!");
    unlink $new_configFL if ( -e $new_configFL );
    open( newFL, ">>$new_configFL" ) or die("Cannot write to $new_configFL:$!");
    while (<oldFL>) {
        if (/^interfaceDefinition/) {
            print newFL "interfaceDefinition = $intDef\n";
        }

        elsif (/^distanceThreshold/) {
            print newFL "distanceThreshold = $atomDistThr\n";
        }
        elsif (/^rasaThr/) {
            print newFL "rasaThr = $rasaThr\n";
        }
        elsif (/^outDir/) {
            print newFL "outDir = $taskDIR\n";
        }
        else {
            print newFL "$_";
        }
    }
    close(oldFL);
    close(newFL);

    print LOG
      "configure file for searchNewPPIDB is generated: $new_configFL.\n";

    #search PPIDB
    #	our $ENVpath=$ENV{'PATH'};
    #	our $perlDIR;
    my $dataDIR = dirname($inputfile);
    $r = rand(1);
    my $tmplogFL = "$dataDIR/temp.$r.log";

    my $command = "perl $searchNewPPIDB_PS_PL  $inputfile $new_configFL";
    $command =~ /(.+)/;

    print LOG "\nCOMMAND: $command\n";
    system("$command > $tmplogFL") == 0
      or die("FAILED: $command:$!");    #search new PPIDB
    &appendFL( $tmplogFL, $logFL );     #"cat $tmplogFL  >> $logFL"; #xue
    unlink($tmplogFL) if ( -e $tmplogFL );    #xue

    if ( -e $new_configFL ) {
        unlink $new_configFL;
        print LOG "$new_configFL removed.\n";
    }


}

sub localAlignmentGen_PS {

    #!usr/bin/perl -w
    #Li Xue
    #4/12/2010

#INPUT FL:
#1. homologs list of query of one query: ../data/1efvB/homologsOfquerySeq.lst
#2. blast output file: e.g., 1efvB.blast
#3. seq_int_homologComplexes file: e.g., $searchNewPPIDBDIR/out/seq_interface_queryComplexes_templates.lst
#
#OUTPUT:
#1. the file of seq_int for homologs of Qry: e.g., ../data/1efvB/seq_interface_homologs.txt
#2. "seq_interface_homologs.localAlign", a fasta file of homologs of query Sequence, in which only the locally aligned part of the homologs are included.
#The output files are put under the same directory with the input files.
#
#perl ./localAlignmentGen.pl ../data/1foe_G_H  1foeG_G_H ../data/1foe_G_H/seq_int_homComplexes_test.lst

    use strict;
    use File::Basename;

    #variables
    my $pairDIR = shift @_;
    my $query   = shift @_;
    my $seq_int_homComplexesFL =
      shift @_;    #  "$searchNewPPIDBDIR/out/seq_interface_templates.lst";

    #my $homologsOfquerySeqLST = "$pairDIR/$query/homologsOfquerySeq.lst";
    my $homologs_in_homComplexesFL =
      "$pairDIR/$query/homologs_in_homComplexes.lst";
    my $QryHHpredFL = "$pairDIR/$query/$query.hhpred";
    my $QrySeqFL    = "$pairDIR/$query/$query.fasta.txt";
    my $seq_int_homologFL =
      "$pairDIR/$query/seq_interface_homologs.txt";    #output 1

    #read the file of homolog list of query
    die("$homologs_in_homComplexesFL does not exist:$!")
      if ( !-e $homologs_in_homComplexesFL );
    my @homologs =
      &readHomologLstFL($homologs_in_homComplexesFL);    #1t9gS|S:A #-- xue-fix

    #read the seq and int for homolog complexes of query
    die("$seq_int_homComplexesFL does not exist:$!")
      if ( !-e $seq_int_homComplexesFL );
    my ($comment, $seqs_complexes, $int_complexes ) =
      &readSeqIntFL($seq_int_homComplexesFL);

    #compare @homologs and @complexes, get the seq and int for homologs of query
    my @complexes = keys( %{$seqs_complexes} );          #1o94F|F:E
    &writeSeqIntHomologs( $seq_int_homologFL, \@homologs, $seqs_complexes,
        $int_complexes );

    die("$QryHHpredFL does not exist:$!") if ( !-e $QryHHpredFL );
    &localAlign_sub_PS( $seq_int_homologFL, $QryHHpredFL, $QrySeqFL );

}

sub collectStat {
#
# Xue, Li
# Oct. 14, 2016
#
# sub:readLocalAlignFL. homologs with interfaces that are all ?s will not be used in prediction.
#
# To collect prediction performance for each test pair if we use one hom-complex interface as the prediction
#
# Input format:
#       Q:P
#
# output: ../data/statistics.txt
# output format:
# >query1 pdbID chainID:
# PDBID+CHAINID Num_Int num_residue1    num_residue2    num_residue1_P Bit_score Positive_Score IdentityScore aligLen_Query aligLen_Homolog CC TP TN FP FN
# ...
# >query2
# ...
#
# >1foeG_G_H      28      377     361     177
# 3bjiA|A:D       24      378     372     177      75.5    3e-14   1e-102  51
#
# perl ./collectStat.pl ../data/seq_int_qryComplexes_perm147pairedChains.lst

    use strict;
    use File::Basename;

    my @QRYpairs = @{ shift @_ };
    my $jobDIR   = shift @_;
    my $outputFL = shift @_;        #"$jobDIR/statistics_$name.txt";

    my $num_pair = scalar @QRYpairs;
    my $seq_Q;                      #sequence of taskID
    my $seq_QP;
    my $int_Q;                      #interface of taskID, array ref
    my $int_QP;
    my $numInt_Q  = 0;              #number of interface residues of query
    my $numInt_QP = 0;
    my $length_Q;     #the length of the query sequence
    my $length_QP;    #the length of the binding partner of the query sequence
    my $flag = 0;

    my $num_residue1_Q;    #$queryLength
    my $num_residue2_Q;    #number of aa that have int info in pdb file

    print LOG "\nCollecting query-homolog sequence similarity information...\n";

    unlink $outputFL if ( -e $outputFL );

    my $intFL_Qry = "$jobDIR/input/seq_int_qryComplexes.lst";
    &writeHeader_statFL( $outputFL, $intFL_Qry )
      ;                    #write the header of the output file

    my $num = 0;           #number of query proteins
    foreach my $pair (@QRYpairs) {

        #'A:B'
        $num++;

        my ( $pairDIR, $chain1, $chain2 ) = &getPairDIR2( $pair, $jobDIR );

        if ( !-d $pairDIR ) {
            print LOG "$pairDIR does not exist. Next.\n";
            next;
        }

        print LOG "query pair: ----------------------------- $pair\n";

        my $localAlignFL1 =
          "$pairDIR/$chain1/seq_interface_homologs.localAlign";
        my $localAlignFL2 =
          "$pairDIR/$chain2/seq_interface_homologs.localAlign";

        &collectStat4onePair( $pair, $localAlignFL1,
            $localAlignFL2, $intFL_Qry, $outputFL );
    }

    print LOG
"Collecting the statistics of the PPI conservation between query sequences and theirs homologs is finished! \n";
    print LOG "There are totally $num_pair pairs tested.\n";
    print LOG "collectStat{}: $outputFL is generated.\n";
}

sub collectStat4onePair {

    # qry: A:B
    #
    # 1. read alignment quality metrics from localAlign files
    # 2. read homolog interface from localAlign files
    # 3. read query interface from intFL_Qry
    #
    #
    our $Rscript;

    my $qryPairID = shift @_;    #A:B
    my $localAlignFL1 =
      shift @_;    #"$pairDIR/$chain1/seq_interface_homologs.localAlign";
    my $localAlignFL2 =
      shift @_;    #"$pairDIR/$chain2/seq_interface_homologs.localAlign";
    my $intFL_Qry = shift @_;    #"$dataDIR/seq_int_qryComplexes.lst";
    my $statFL    = shift @_;    #output file

    #---- read alignments between A and homologs of A
    my ( $sbjIDs1, $len_A, $localAlignResults1 ) =
      &readLocalAlignFL($localAlignFL1)
      ;    #$alignmentLength->{"$homolog,$groupID"} =23; $homolog = '1o96Q|Q:Z'


    my $LALs1         = $localAlignResults1->{'LAL'};              #-- with gaps
    my $sbjctLength1  = $localAlignResults1->{'sbjctLen'};
    my $start_Hs1     = $localAlignResults1->{'start_S'};
    my $end_Hs1       = $localAlignResults1->{'end_S'};
    my $start_Qs1     = $localAlignResults1->{'start_Q'};
    my $end_Qs1       = $localAlignResults1->{'end_Q'};
    my $BitScore1     = $localAlignResults1->{'BitScore'};
    my $EVal1         = $localAlignResults1->{'EVal'};
    my $SID1          = $localAlignResults1->{'Identities'};
    my $Positives1    = $localAlignResults1->{'Positives'};
    my $Similarities1 = $localAlignResults1->{'Similarity'};
    my $hhpredProb1   = $localAlignResults1->{'hhpredProb'};
    my $hhpred_pVal1  = $localAlignResults1->{'hhpred_Pval'};
    my $hhpred_SS1    = $localAlignResults1->{'hhpred_SS'};
    my $alignedSeq_A  = $localAlignResults1->{'QrySeqAligned'};
    my $int_H1        = $localAlignResults1->{'HomologIntAligned'};

    #---- read alignments between B and homologs of B
    my ( $sbjIDs2, $len_B, $localAlignResults2 ) =
      &readLocalAlignFL($localAlignFL2)
      ;    #$alignmentLength->{"$homolog,$groupID"} =23; $homolog = '1o96Q|Q:Z'

    my $LALs2         = $localAlignResults2->{'LAL'};              #-- with gaps
    my $sbjctLength2  = $localAlignResults2->{'sbjctLen'};
    my $start_Hs2     = $localAlignResults2->{'start_S'};
    my $end_Hs2       = $localAlignResults2->{'end_S'};
    my $start_Qs2     = $localAlignResults2->{'start_Q'};
    my $end_Qs2       = $localAlignResults2->{'end_Q'};
    my $BitScore2     = $localAlignResults2->{'BitScore'};
    my $EVal2         = $localAlignResults2->{'EVal'};
    my $SID2          = $localAlignResults2->{'Identities'};
    my $Positives2    = $localAlignResults2->{'Positives'};
    my $Similarities2 = $localAlignResults2->{'Similarity'};
    my $hhpredProb2   = $localAlignResults2->{'hhpredProb'};
    my $hhpred_pVal2  = $localAlignResults2->{'hhpred_Pval'};
    my $hhpred_SS2    = $localAlignResults2->{'hhpred_SS'};
    my $alignedSeq_B  = $localAlignResults2->{'QrySeqAligned'};
    my $int_H2        = $localAlignResults2->{'HomologIntAligned'};

    #---- read qry interfaces
    my ( $A, $B ) = split( /:/, $qryPairID );
    my ( $comment, $seqQry, $intQry ) = &readSeqIntFL($intFL_Qry);
    my ( $num_residue1, $num_residue2, $num_int ) =
      &getIntNum($intFL_Qry)    #-- xue-fix
      ;

   # num_residue1: whole seq; num_residue2: structural seq
    if (   $len_A ne $num_residue1->{"$A|$A:$B"}
        || $len_B ne $num_residue1->{"$B|$B:$A"} )
    {
        my $resNum_A = $num_residue1->{"$A|$A:$B"};
        my $resNum_B = $num_residue1->{"$B|$B:$A"};

        die(
"qry ($A or $B) length ($len_A, $len_B) in local align file is different from $intFL_Qry ($resNum_A, $resNum_B) :$!"
        );
    }

    #---- For qry A, compare int_H and int_Q
    my ( $TP_1, $TN_1, $FP_1, $FN_1, $CC_1, $specificity_1, $sensitivity_1,
        $accuracy_1 )
      = &compare_qryInt_homologInt( $len_A, $alignedSeq_A,
        $intQry->{"$A|$A:$B"}, $sbjIDs1, $int_H1, $start_Qs1, $end_Qs1 );

    #- $TP_1->{"$homolog,$groupID"} = 12

    #---- For qry B, compare int_H and int_Q
    my ( $TP_2, $TN_2, $FP_2, $FN_2, $CC_2, $specificity_2, $sensitivity_2,
        $accuracy_2 )
      = &compare_qryInt_homologInt( $len_B, $alignedSeq_B,
        $intQry->{"$B|$B:$A"}, $sbjIDs2, $int_H2, $start_Qs2, $end_Qs2 );

    #- $TP_2->{"$homolog,$groupID"} = 14

    #------------ write the statistics file
    open( OUTPUT, ">>$statFL" ) or die("Cannot open $statFL:$!");
    foreach my $homolog_A ( @{$sbjIDs1} ) {
        foreach my $homolog_B ( @{$sbjIDs2} ) {

            # $homolog_A = '4dgaA|A:C,3'
            # $homolog_B = '1ak4D|D:E,80'

            #-------------
            # Check whether $homolog_A and $homolog_B is an interacting pair
            my $pdbID = substr( $homolog_A, 0, 4 );
            my ( $homolog1_tmp, $chnID, $chnID_partner, $groupID_tmp ) =
              split( /[\|,\:]/, $homolog_A );
            if ( $homolog_B !~ /$pdbID$chnID_partner\|$chnID_partner:$chnID,/ )
            {
                next;
            }

            #-------------

            my $len_H1 = $sbjctLength1->{$homolog_A};    # homolog for query A
            my $len_H2 = $sbjctLength2->{$homolog_B};    # homolog for query B

            my $numInt_Homolog_A = &countInt( $int_H1->{$homolog_A} )
              ;    #int of the aligned part of homolog A
            my $numInt_Homolog_B = &countInt( $int_H2->{$homolog_B} )
              ;    #int of the aligned part of homolog B

            #-- LAL without gaps for A and A'
            my $LAL1_Q =
              $end_Qs1->{$homolog_A} -
              $start_Qs1->{$homolog_A} + 1;    #-- without gap
            my $LAL1_S =
              $end_Hs1->{$homolog_A} -
              $start_Hs1->{$homolog_A} + 1;    #-- without gap

            #-- LAL without gaps for B and B'
            my $LAL2_Q =
              $end_Qs2->{$homolog_B} -
              $start_Qs2->{$homolog_B} + 1;    #-- without gap
            my $LAL2_S =
              $end_Hs2->{$homolog_B} -
              $start_Hs2->{$homolog_B} + 1;    #-- without gap

            my $frac_LAL1 =
              sprintf( "%.2f", ( $LAL1_Q / $len_A ) * ( $LAL1_S / $len_H1 ) );
            my $frac_LAL2 =
              sprintf( "%.2f", ( $LAL2_Q / $len_B ) * ( $LAL2_S / $len_H2 ) );

            #-- change 4dgaA|A:C,3:4dgaC|C:A,39 => 4dgaA,3:4dgaC,39
            my $homologPairID = "$homolog_A:$homolog_B";
            $homologPairID =~ s/\|\w+:\w+,/,/g;

            #-- temperarily put a column of pred_CC in the statFL
            my $pred_CC = 'nan';

            #--
            my @line_tmp = (
                $qryPairID,                   $homologPairID,
                $pred_CC,                     $len_A,
                $len_B,                       $len_H1,
                $len_H2,                      $num_int->{"$A|$A:$B"},
                $num_int->{"$B|$B:$A"},       $numInt_Homolog_A,
                $numInt_Homolog_B,            $EVal1->{$homolog_A},
                $EVal2->{$homolog_B},         $SID1->{$homolog_A},
                $SID2->{$homolog_B},          $Positives1->{$homolog_A},
                $Positives2->{$homolog_B},    $Similarities1->{$homolog_A},
                $Similarities2->{$homolog_B}, $hhpredProb1->{$homolog_A},
                $hhpredProb2->{$homolog_B},   $hhpred_pVal1->{$homolog_A},
                $hhpred_pVal2->{$homolog_B},  $hhpred_SS1->{$homolog_A},
                $hhpred_SS2->{$homolog_B},    $LAL1_Q,
                $LAL1_S,                      $LAL2_Q,
                $LAL2_S,                      $frac_LAL1,
                $frac_LAL2,                   $TP_1->{$homolog_A},
                $TN_1->{$homolog_A},          $FP_1->{$homolog_A},
                $FN_1->{$homolog_A},          $CC_1->{$homolog_A},
                $specificity_1->{$homolog_A}, $sensitivity_1->{$homolog_A},
                $TP_2->{$homolog_B},          $TN_2->{$homolog_B},
                $FP_2->{$homolog_B},          $FN_2->{$homolog_B},
                $CC_2->{$homolog_B},          $specificity_2->{$homolog_B},
                $sensitivity_2->{$homolog_B}
            );

            if ( !@line_tmp ) {
                die(
"Stat line not defined for qry $qryPairID and homolog $homologPairID:$!"
                );
            }

#            #---------------------------------------------------------
#            #--- for double checking
#            my $numInt_QryA = $num_int->{"$A|$A:$B"};
#            my $numInt_QryB = $num_int->{"$B|$B:$A"};
#
#            print "\n$qryPairID,\t$homolog_A:$homolog_B,
#                len_A: $len_A (whole seq),\tlen_B: $len_B (whole seq),
#                len_H1: $len_H1 (read from local_align file),\tlen_H2: $len_H2 (read from local_align file),
#                int_Qry1: $numInt_QryA,\tint_Qry2:$numInt_QryB,
#                int_H1: $numInt_Homolog_A,\tint_H2: $numInt_Homolog_B,
#                Eval1: $EVal1->{$homolog_A},\tEval2: $EVal2->{$homolog_B},
#                SID1: $SID1->{$homolog_A},\tSID2: $SID2->{$homolog_B},
#                postiveS: $Positives1->{$homolog_A},\tpostiveS: $Positives2->{$homolog_B},
#                hhpredProb: $hhpredProb1->{$homolog_A},\thhpredProb: $hhpredProb2->{$homolog_B},
#                LAL1: $LAL1,\tLAL2: $LAL2,
#                frac_LAL1: $frac_LAL1 ($LAL1/$len_A * $LAL1/$len_H1),\tfrac_LAL2: $frac_LAL2($LAL2/$len_B * $LAL2/$len_H2)\n";
#
#            #--- checking done.
#            #---------------------------------------------------------

            my $line = join( "\t", @line_tmp );
            print OUTPUT "$line\n";

        }
    }
    close OUTPUT;

}

#------------------------
sub rmSameProtPair_deleteFL {

    #Li Xue
    #Jul 14, 2012

#
#remove the homo-interolog that maybe the same interacting domain pair as the query pair
#rm the homolog that from intput file and write the remaining into the output file.

   #rmSameProtPair_deleteFL.pl :
   #	Given qry domain pair A:B and its homologous domain pair A':B',
   #	if A':B' are part of the homo-interologs in the deleteFL, A':B' is removed.

    # Input: statistics.txt and deleteFL
    # Output:statistics.txt without homo-interologs in deleteFL

# INPUT 1 (statistics.txt):
#
#   qryPairID  homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2   EVal1   EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL1    LAL2    frac_LAL1   frac_LAL2   TP_1    TN_1    FP_1    FN_1    CC_1    specificity_1sensitivity_1    TP_2    TN_2    FP_2    FN_2    CC_2    specificity_2   sensitivity_2
#   Q:P 2x2dC,42:2x2dD,69   164 137 165 147 Nan Nan 14  8   1.1E-45 3E-41   99.3    98.599.3    100.0   4.7 4.4 100.0   100.0 164 136 0.99    0.92
#   Q:P 2x2dC,42:2x2dD,38   164 137 165 147 Nan Nan 14  8   1.1E-45 3.7E-68 99.3    98.599.3    100.0   4.7 4.4 100.0   100.0 164 136 0.99    0.92

    # INPUT 2: (DeleteFL: Delete file contains Highly similar homo-interologs.)
    #	   1a2kA:1a2kC =>
    #	   1a14H:1a14N => 1NMB_H:1NMB_N,1A14_H:1A14_N,1NMA_H:1NMA_N
    #	   A:B => 3fcl_B:*,*:1lqm_A

#perl rmSameProtPair_deleteFL.pl ../data/test/statistics_test_protPair.txt  ../data/test/delete.lst

    use strict;
    use File::Basename;

    my $statFL_ori = shift @_;    #'statistics.txt'
    my $deleteFL   = shift @_;

    my $flag_Homolog = 1;

#$flag_Homolog=0: at least one of the protein pairs have homo-interologs in $FullStatFL
#$flag_Homolog=1: none of the query pairs have homo-interologs in $FullStatFL

    print LOG
"\n\n\n\n-------------------------------------------------------------\n\n";
    print LOG
      "\n\tNow remove highly similar homolog pairs from $statFL_ori....\n";
    print LOG "\tHighly similar homologs are listed in $deleteFL\n\n";
    print LOG
      "\n-------------------------------------------------------------\n\n\n";

    my $dirname = dirname($statFL_ori);
    my $basename = basename( $statFL_ori, ( '.txt', '.lst' ) );
    my $outputFL =
      "$dirname/$basename\_wo_sameProt.txt";    #'statistics_wo_sameProt.txt'
    my $num_sameProt = 0;

    #-------------------
    my %homoInterologsToBeDelMap =
      %{ &readDeleteFL($deleteFL) }
      ;    # $homoInterologsToBeDel->{QryA:QryB} = '1NMBH:1NMBN,1A14H:*,*:1NMAN'

    #--delete highly similar homolog pairs and write output file
    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    my $this_function = ( caller(0) )[3];
    print OUTPUT "#Generated by $this_function in PSHomPPI_resiPairs.pm\n";
    print OUTPUT "#Homologs listed in $deleteFL are removed from $statFL_ori\n";

    my $flag_inputFormat = 0;
    open( INPUT, "<$statFL_ori" ) or die("Cannot open $statFL_ori:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^qryPairID/) {
            print OUTPUT "$_\n";
            next;
        }
        if (/^(\w+):(\w+)\s+(\w{5},\d+):(\w{5},\d+)\s+/) {

#   Q:P 2x2dC,42:2x2dD,69   164 137 165 147 Nan Nan 14  8   1.1E-45 3E-41   99.3    98.599.3    100.0   4.7 4.4 100.0   100.0 164 136 0.99    0.92

            $flag_inputFormat = 1;

            my $qryA = $1;
            my $qryB = $2;
            my ( $protAprim, $groupID1 ) = split( /,/, $3 );    # $3 = 2x2dC,34
            my ( $protBprim, $groupID2 ) = split( /,/, $4 );    # $4 = 2x2dD,12

            #-- get the delete list for a specific qry pair
            my @homologsToBeDel;
            if ( defined $homoInterologsToBeDelMap{"$qryA:$qryB"} ) {
                @homologsToBeDel =
                  @{ $homoInterologsToBeDelMap{"$qryA:$qryB"} };
            }
            elsif ( defined $homoInterologsToBeDelMap{"$qryB:$qryA"} ) {
                @homologsToBeDel =
                  @{ $homoInterologsToBeDelMap{"$qryB:$qryA"} };
            }
            else {
                $flag_Homolog = 0;
                print LOG
"This query pair ($qryA:$qryB) has no homologs to be deleted (i.e., no homologs for this qry pair defined in $deleteFL).\n";
                print OUTPUT "$_\n";
                next;
            }

#            print "\nhomologs to be del for qry $qryB:$qryA: @homologsToBeDel\n";
#            print "A':B': $protAprim:$protBprim\n";
#            print "Check whether A':B' should be deleted\n";

            my $ans = &belong_new( "$protAprim:$protBprim", \@homologsToBeDel );

            if ( $ans == 0 ) {
                print OUTPUT "$_\n";
                $flag_Homolog = 0;
            }
            elsif ( $ans == 1 ) {
                $num_sameProt++;
                print LOG
"$protAprim,$groupID1:$protBprim,$groupID2 is removed from $statFL_ori\n";
            }
        }

    }
    close INPUT;
    close OUTPUT;

    if ( $flag_inputFormat == 0 ) {
        die("Input stat file ($statFL_ori) format is wrong:$!");
    }

    print LOG
"$outputFL is generated (flag_Homolog = $flag_Homolog. flag_Homolog= 0 means at least one of the protein pairs have homo-interologs in $outputFL). $num_sameProt homologs that maybe the same protein of the qry are not writen in this file.\n\n";

    return $flag_Homolog;

#$flag_Homolog=0: at least one of the protein pairs have homo-interologs in $FullStatFL
#$flag_Homolog=1: none of the query pairs have homo-interologs in $FullStatFL
}

sub readPairFL {
    my $protPairFL = shift @_;
    my @pairs;

    open( INPUT, "<$protPairFL" ) || die("Cannot open $protPairFL:$!");

    while (<INPUT>) {

        s/[\n\r]//mg;
        if (/^[\-\w\.]+:[\-\w\.]+/) {
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
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^\s*$/) {
            next;
        }
        if (/^>(.+)/) {
        #	>Actin

            if ( defined $ID && defined $seq ) {
                $seqs{$ID} = $seq;
                #                print LOG "$ID:$seq\n";
            }

            $seq = '';
            $ID  = $1;
            $ID  =~ s/\s+//g;
            next;
        }
        if (/^[a-zA-Z]+/) {
            s/\s+//;
            $seq = $seq . $_;
        }
    }
    close INPUT;

    #read the last protein
    if ( defined $ID && defined $seq ) {
        $seqs{$ID} = $seq;
        #        print LOG "$ID:$seq\n";
    }


    if ( !%seqs ) {
        die("No seq read from $protSeqFL. Check $protSeqFL:$!");
    }

    return \%seqs;

}

sub uniquePair {

    #A:B and B:A are considered the same pair

    my @pairs = @{ shift @_ };    #(A:B, B:A, B:C, E:F)

    #	print "\nThe original pairs are: @pairs\n\n";

    my @uniquePairs;
    my %seen;

    foreach my $pair (@pairs) {
        my ( $A, $B ) = split( /:/, $pair );

        if ( !defined $seen{"$A:$B"} && !defined $seen{"$B:$A"} ) {
            $seen{$pair} = 1;
        }

    }

    @uniquePairs = keys %seen;

    #	print "\nThe unique pairs are: @uniquePairs\n\n";

    return \@uniquePairs;

}

sub getPairDIR {
    my $pair    = shift @_;    #'A:B'
    my $dataDIR = shift @_;

    if ( !defined $dataDIR ) {
        die("dataDIR not defined:$!");
    }
    my $pairDIR;

    my ( $chain1, $chain2 ) = split( /:/, $pair );
    if ( !defined $chain1 || !defined $chain2 ) {
        die("Check $pair:$!");
    }

    my $pairDIR1 = "$dataDIR/$chain1:$chain2";    #possible pairDIR
    my $pairDIR2 = "$dataDIR/$chain2:$chain1";    #another possible pairDIR

    if ( $chain1 lt $chain2 ) {
        $pairDIR = $pairDIR1;
    }
    else {
        $pairDIR = $pairDIR2;
    }

    return ( $pairDIR, $chain1, $chain2 );

}

#-----------

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

sub write_hhpredConfigFL {

# For settings.config, use variables in /home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/runquery_abs.py
# See a config FL template: '/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/configExample.txt'

    my $configFL    = shift @_;
    my $hhpredParam = shift @_;    # hash ref

    #--- write the configFL based on user-specified parameters

    unlink $configFL if ( -e $configFL );
    open( CONFIG, ">>$configFL" ) or die("Cannot open $configFL:$!");
    foreach my $key ( keys %{$hhpredParam} ) {
        print CONFIG "$key = $hhpredParam->{$key}\n";
    }
    close CONFIG;
    print LOG "hhpred config file generated: $configFL .\n";
}

sub readHHpredOUTPUTFL_userProtDB {

#read the last round of iteration of PSI-BLAST results
#this subroutine is for *.blast file with format when Blast searches against pdb database.

 #search *.blast for the local alignment of homolog in PPIDB with query sequence
 #return $alignmentLength, $SbjctLength, $BitScore_Evalue, $Identities, etc.

# 1. Input file example,
#
#	 # PSIBLAST 2.2.25+
#	# Iteration: 1
#	# Query: 1avxA
#	# Database: C:\Users\lixue/ncbi-blast-2.2.25+/data/scop_v175_prot/scop_v175_prot
#	# Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, query acc.
#	# 485 hits found
#	pdb|1tx6|D	223	223	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIDVLEGNEQFINAAKIITHPNFNGNTLDNDIMLIKLSSPATLNSRVATVSLPRSCAAAGTECLISGWGNTKSSGSSYPSLLQCLKAPVLSDSSCKSSYPGQITGNMICVGFLEGGKDSCQGDSGGPVVCNGQLQGIVSWGYGCAQKNKPGVYTKVCNYVNWIQQTIAAN	1	223	1	223	 413	2e-116	223	100.00	100.00	1avxA
#
#	# PSIBLAST 2.2.27+
#	# Iteration: 1
#	# Query: PFK1
#	# Database: /home/lixue/programs/ncbi-blast-2.2.27+/data/nr_pdbaa_s2c/nr_pdbaa_s2c
#	# Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, query acc.
#	# 7 hits found
#	gi|321159870|pdb|3O8O|A;gi|321159872|pdb|3O8O|C;gi|321159874|pdb|3O8O|E;gi|321159876|pdb|3O8O|G 27      787     EYNKIGDILSGRLKLRAEVAALAAENK     EYNKIGDILSGRLKLRAEVAALAAENK     1       27      761     787     55.8    2e-10   27      100.00  100.00  PFK1
#	gi|169791680|pdb|2P4X|A;gi|169791681|pdb|2P4X|B 27      239     LKLRAEVAAL      LKLRAEVASL      13      22      152     161     27.8    0.31    10      90.00   90.00   PFK1
#	# BLAST processed 1 queries
#
# 2. When chain ID is in lowercase, Blast will output double uppercase chainID
# For example: pdb|3KIX|LL  Chain l,  >pdb|2QNH|MM Chain m,  >pdb|3KIS|LL Chain ...   259    3e-70
# So we need to extract chain ID from "Chain l" part not from "pdb|3KIX|LL".

    my $hhpredFL = shift @_;
    my $blastResults;
    my $qryID;
    my $qryLenth;
    my $flag = 0;    #0: at least one hit found. 1: no hit found.
    my $num  = 0;    # the number of alignment pairs in the hhpredFL

    open( INPUT, "<$hhpredFL" ) || die("Cannot open $hhpredFL:$!");
    while (<INPUT>) {

        if (/^# 0 hits found/i) {
            $flag = 1;
            last;
        }

        if (/^[(gi)]{0,}.+pdb\|/) {

#pdb|1si8|C      NGDKEQVAKWVNLAQKELVIKNFAKLSQSLETLDSQL   HGFGSHTFKWVNAAGEVFFVKYHFKTNQGIKNLESQL   68      104     189     225     26.1    8.9     37      35.14   54.05   2hrkA
#gi|321159870|pdb|3O8O|A;gi|321159872|pdb|3O8O|C;gi|321159874|pdb|3O8O|E;gi|321159876|pdb|3O8O|G 27      787     EYNKIGDILSGRLKLRAEVAALAAENK     EYNKIGDILSGRLKLRAEVAALAAENK     1       27      761     787     55.8    2e-10   27      100.00  100.00  PFK1

            #-- check format of the input file
            my @tmp = split( /\s+/, $_ );

            if ( scalar @tmp ne 19 ) {
                my $num = scalar @tmp;
                die(
"ERROR: The current line: $_\nERROR: The column number of the input file ($hhpredFL) is not 17. It has $num columns: $!"
                );

            }

            $num++;

            #--

            my (
                $sAllSeqID, $qryLen,      $slen,       $qseq,
                $sseq,      $qstart,      $qend,       $sstart,
                $send,      $bitscore,    $evalue,     $LAL,
                $pident,    $positiveS,   $similarity, $hhpred_prob,
                $qryID1,    $hhpred_pVal, $hhpredSS
            ) = split( /\s+/, $_ );

            #
            my @homologIDs = split( /;/, $sAllSeqID );

            foreach my $homologID (@homologIDs) {

                #				$homologID='pdb|1si8|C'
                ($homologID) = $homologID =~ /pdb\|(\w{4}\|\w+)/;
                my ( $pdbID, $chnID ) = split( /\|/, $homologID );

                if ( length($chnID) == 2 ) {

# When chain ID is in lowercase, Blast will output double uppercase chainID
# For example: pdb|3KIX|LL  Chain l,  >pdb|2QNH|MM Chain m,  >pdb|3KIS|LL Chain ...   259    3e-70

                    $chnID = lc( substr( $chnID, 0, 1 ) );

                }

                $pdbID     = lc($pdbID);
                $homologID = "$pdbID$chnID";

                #$positive_num: the total number of positive matches
                #$LAL: local alignment length
                $blastResults->{$homologID}->{$num}->{'slen'}     = $slen;
                $blastResults->{$homologID}->{$num}->{'qseq'}     = $qseq;
                $blastResults->{$homologID}->{$num}->{'sseq'}     = $sseq;
                $blastResults->{$homologID}->{$num}->{'qstart'}   = $qstart;
                $blastResults->{$homologID}->{$num}->{'qend'}     = $qend;
                $blastResults->{$homologID}->{$num}->{'sstart'}   = $sstart;
                $blastResults->{$homologID}->{$num}->{'send'}     = $send;
                $blastResults->{$homologID}->{$num}->{'bitscore'} = $bitscore;
                $blastResults->{$homologID}->{$num}->{'evalue'}   = $evalue;
                $blastResults->{$homologID}->{$num}->{'LAL'}      = $LAL;
                $blastResults->{$homologID}->{$num}->{'identityS'} =
                  $pident;    #[0-100]
                $blastResults->{$homologID}->{$num}->{'positiveS'} =
                  $positiveS;    #        [0-100]
                $blastResults->{$homologID}->{$num}->{'similarity'} =
                  $similarity;    #        [0-100]
                $blastResults->{$homologID}->{$num}->{'hhpredProb'} =
                  $hhpred_prob;    #[0-100]

                $blastResults->{$homologID}->{$num}->{'hhpred_Pval'} =
                  $hhpred_pVal;
                $blastResults->{$homologID}->{$num}->{'hhpred_SS'} = $hhpredSS;
            }

            $qryID    = $qryID1;
            $qryLenth = $qryLen;
        }
    }
    close INPUT;

    if ( $flag == 1 ) {

        #No hit in this blast file.

        return ( $qryID, $qryLenth, $blastResults );

    }

    if ( !defined $qryID ) {

        die("qryID undefined. Check $hhpredFL:$!");
    }

    if ( !defined $qryLenth ) {

        die("qryLenth undefined. Check $hhpredFL:$!");

    }

    if ( !defined $blastResults ) {

        die("Variable blastResults undefined. Check $hhpredFL:$!");
    }

    return ( $qryID, $qryLenth, $blastResults );

}

sub writeHomologLstfile {
    my ( $outfile, $homologRef, $hhpredProbThr ) = @_;
    unlink $outfile if ( -e $outfile );

    open( OUTFILE, ">>$outfile" ) || die("Cannot open $outfile:$!");

    my @homologs = @{$homologRef};
    my $num_seq  = scalar @homologs;

    #write header
    print OUTFILE
"#The number of homolog chains in the last round with >= hhpred prob. $hhpredProbThr is $num_seq\n";
    print OUTFILE "#Generated by parse_psiblast_userDB.pl\n";

    #write homolog IDs
    foreach my $homolog (@homologs) {

        print OUTFILE "$homolog\n";
    }

    close(OUTFILE);

    print LOG "$outfile generated.\n";

}

#------------
sub unique {
    my @a = @_;
    my %seen;
    @seen{@a} = 1 x scalar(@a);
    @a = keys(%seen);
    return @a;

}

sub readPDBHomologLstFL {
    use warnings;
    use warnings::register;
    my $homologLstFL = shift @_;
    my $list;
    open( INPUT, "<$homologLstFL" ) || die("Cannot open $homologLstFL:$!");
    foreach (<INPUT>) {

        #1auiB
        s/[\r\n]//mg;
        if (/^(\w+)/) {

            my $homolog = $1;
            my $pdbID   = substr( $homolog, 0, 4 );
            my $chainID = substr( $homolog, 4, 1 );

            push @{ $list->{$pdbID} }, $chainID;
        }
    }
    close(INPUT);

    if ( !defined $list ) {
        print LOG "\n\nWarning: no homologs read from $homologLstFL:$!";
#        warnings::warn("\n\nWarning: no homologs read from $homologLstFL:$!");
    }

    return $list;
}

sub headerHomcomplexFL {

    my $homComplexesFL = shift @_;

    unlink $homComplexesFL if ( -e $homComplexesFL );

    open( OUTPUT, ">>$homComplexesFL" )
      || die("Cannot open $homComplexesFL:$!");
    print OUTPUT "#Generated by compareHomlogLsts() in PSHomPPI_new.pm.\n";
    print OUTPUT
"#The chain IDs of the hom-complexes are ordered: first column of chains are A', and 2nd column of chains are B'.\n";
    close OUTPUT;

}

sub reformHomComplexPairs {

    #2greB:2greO -> 2gre B:O

    my @homComplexes = @{ shift @_ };
    my @reformed_homComplexes;

    foreach my $homComplex (@homComplexes) {
        my ( $protA, $protB ) = split( /:/, $homComplex );
        my $pdbID = substr( $protA, 0, 4 );
        my $chnA  = substr( $protA, 4, 1 );
        my $chnB  = substr( $protB, 4, 1 );

        push @reformed_homComplexes, "$pdbID\t$chnA:$chnB";
    }

    return \@reformed_homComplexes;

}

sub writeArrayIntoFL {

    #write the element of input array into a file.
    #one element one row.
    #
    my @inputArray = @{ shift @_ };   #$Qry1homologs_in_homComplexes = shift @_;
    my $outputFL = shift @_;    #	my $homologs_in_homComplexesFL   = shift @_;

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    foreach (@inputArray) {
        print OUTPUT "$_\n";
    }
    close(OUTPUT);
    print LOG "$outputFL generated.\n";

}

sub appendFL {
    use strict;

    #append file1 to the end of file2
    my $file1 = shift @_;
    my $file2 = shift @_;

    open( FL2, ">>$file2" ) || die("Cannot open $file2:$!");
    open( FL1, "<$file1" )  || die("Cannot open $file1:$!");

    foreach (<FL1>) {
        print FL2 $_;
    }
    close FL1;
    close FL2;

}

sub readHomologLstFL {

    # INPUT file (homologs_in_homComplexes.lst):
    #    1m9fA|A:D
    #    1m9yF|F:G

    my $homologsOfquerySeqLST = shift @_;
    my @homologs;
    open( homologsOfquerySeqLST, "<$homologsOfquerySeqLST" )
      || die("Cannot open $homologsOfquerySeqLST:$!");
    while (<homologsOfquerySeqLST>) {
        s/[\n\r]//mg;
        s/\s+//g;

        if (/^\w{5,6}\|\w{0,}:\w{0,}$/) {

            #		1t9gS|S:A

            push @homologs, "$_";    #1t9gS|S:A
        }
    }
    close(homologsOfquerySeqLST);

    if ( !@homologs ) {
        die("Nothing read from $homologsOfquerySeqLST:$!");
    }
    return @homologs;
}

sub writeSeqIntHomologs {

    # write the file of the seq and int for homologsof query

    our $intNum_cutoff;

    my $outputFL       = shift @_;
    my @homologs       = @{ shift @_ };    #1o94F|F:E
    my $seqs_complexes = shift @_;
    my $int_complexes  = shift @_;

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" );

    foreach my $homolog (@homologs) {
        if ( !$seqs_complexes->{$homolog} ) {

            print LOG
"localAlignmentGen.pl: $homolog or its binding partner has <= $intNum_cutoff int aa, and it is excluded in the local alignment file.\n";
        }
        else {

            print OUTPUT ">$homolog\n";
            print OUTPUT "$seqs_complexes->{$homolog}\n";
            print OUTPUT "$int_complexes->{$homolog}\n";
        }
    }
    close(OUTPUT);
}

#------------------------------
sub localAlign_sub_PS {

# author: Xue, Li
# date: Aug26, 2012
# This script is the pipeline of partner-specific PPI
# This script will be called by localAlignmentGen.pl
#
# Input:
# 1. "seq_interface_homologs.lst", which contains seq and interface info
# 2. BLAST output file: "$taskID.blast", which contains alignment info
#
# Ouput:
# "seq_interface_homologs.localAlign", a fasta file of homologs of query Sequence, in which only the locally aligned part of the homologs are included.
# The output file is put under the same directory with the input files.
#
# Usage:
# e.g.:  perl ./localAlign_sub.pl ../data/test/d2ez8b2-d2ez8b3/d2ez8b2/seq_interface_homologs.txt ../data/test/d2ez8b2-d2ez8b3/d2ez8b2/d2ez8b2.blast
# perl localAlign_sub.pl ../data/trans212/d1cs4b_-d1cs4c2/d1cs4c2/seq_interface_homologs.txt ../data/trans212/d1cs4b_-d1cs4c2/d1cs4c2/d1cs4c2.blast
# perl localAlign_sub.pl ../data/trans212/d1k5da_-d1k5db_/d1k5db_/seq_interface_homologs.txt ../data/trans212/d1k5da_-d1k5db_/d1k5db_/d1k5db_.blast

    use strict;
    use File::Basename;

    #Files
    my $SeqIntFL =
      shift @_;    #seq_interface file for the homologs (whole seq of homologs)
    my $BlastFL  = shift @_;
    my $QrySeqFL = shift @_;

    #Directories
    my $ProteinToBePredictedDIR = dirname($BlastFL);    #xue

    my $outputFL = "$ProteinToBePredictedDIR/seq_interface_homologs.localAlign";
    unlink $outputFL if ( -e $outputFL );

    #Program begin

    &Write_header_LocalAlignFL($outputFL);

    my ( $qryID, $totalLengthQry, $blastResults ) =
      &readHHpredOUTPUTFL_userProtDB($BlastFL);

    my ( $comment, $seqs_H, $ints_H ) =
      &readSeqIntFL($SeqIntFL);    #- $seqs_H->{'1ACBA|A:C'}=seq;

    open( OUTPUT, ">>$outputFL" ) or die("Cannot write to $outputFL:$!");
    print OUTPUT "Query: $totalLengthQry letters\n";

    foreach my $homologPartnerSpc ( keys %{$seqs_H} ) {

        # $homologPartnerSpc = '1ahjA|A:C'
        # $homologPartnerSpc = d1nb5.1|d1nb5.1:d1nb5i_
        my ($homolog) = $homologPartnerSpc =~ m/^([^\|]+)\|[^:]+:[^:]+/;

        print LOG
          "Searching the info of homolog $homologPartnerSpc in $BlastFL ...\n";
        foreach my $num ( keys %{ $blastResults->{$homolog} } ) {

            my $SbjctLength = $blastResults->{$homolog}->{$num}->{'slen'};

            my (
                $alignmentLength, $BitScore,     $Evalue,
                $Identities,      $Positives,    $similarity,
                $start_Q,         $end_Q,        $start_S,
                $end_S,           $localAlign_Q, $localAlign_S,
                $hhpred_prob,     $hhpred_pVal,  $hhpredSS
              )
              = (

                $blastResults->{$homolog}->{$num}->{'LAL'},
                $blastResults->{$homolog}->{$num}->{'bitscore'},
                $blastResults->{$homolog}->{$num}->{'evalue'},
                $blastResults->{$homolog}->{$num}->{'identityS'},
                $blastResults->{$homolog}->{$num}->{'positiveS'},
                $blastResults->{$homolog}->{$num}->{'similarity'},
                $blastResults->{$homolog}->{$num}->{'qstart'},
                $blastResults->{$homolog}->{$num}->{'qend'},
                $blastResults->{$homolog}->{$num}->{'sstart'},
                $blastResults->{$homolog}->{$num}->{'send'},
                $blastResults->{$homolog}->{$num}->{'qseq'},
                $blastResults->{$homolog}->{$num}->{'sseq'},
                $blastResults->{$homolog}->{$num}->{'hhpredProb'},
                $blastResults->{$homolog}->{$num}->{'hhpred_Pval'},
                $blastResults->{$homolog}->{$num}->{'hhpred_SS'}
              );

            my $groupID = $num;
            print OUTPUT
">$homologPartnerSpc,GroupID $groupID,Alignment_Length $alignmentLength,Sbjct_Length $SbjctLength,Sbjct_RESI $start_S - $end_S,Query_RESI $start_Q - $end_Q,Score = $BitScore bits, Expect = $Evalue, Identities = $Identities%, Positives = $Positives%, Similarity = $similarity, hhpred_prob = $hhpred_prob, hhpred_P-Val = $hhpred_pVal, hhpred_SS = $hhpredSS\n";
            print OUTPUT "$localAlign_Q\n";

#get the interface of the part of homologs that is aligned with the query sequence by BLASTp

            if ( !$start_S || !$end_S ) {
                die(
"Debug: start_S or end_S is not initiated. Check localAlign_sub.pl:$!"
                );
            }

            my $interface =
              &subSbjctInt( $SeqIntFL, $homologPartnerSpc, $start_S, $end_S );

            #insert spaces into $interface according to the BLASTp alignment
            $interface = &insertSpc( $localAlign_S, $interface );
            print OUTPUT "$interface\n";

            #		print "AlignedPart_query:          $localAlign_Q\n";
            #		print "AlignedPart_homolog:        $localAlign_S\n";
            #		print "int_homolog_aligned:$interface\n\n";

            if ( length($localAlign_Q) ne length($interface) ) {
                my $len_localAlignQ = length($localAlign_Q);
                my $len_interface   = length($interface);
                die(
"The length of localAlign_Q ($len_localAlignQ) is different from interface ($len_interface). Check localAlign_sub.pl:$!"
                );
            }

        }

    }

    close(OUTPUT);
    print LOG
"Generation of localAlignment Fasta file is done! The output file is $outputFL.\n\n";
}

sub writeHeader_statFL {

    my $outputFL  = shift @_;
    my $intFL_qry = shift @_;

    my $ans = &existIntInfo($intFL_qry);

    my @variables_tmp;

    if ( $ans == 1 ) {

        @variables_tmp = (
            "qryPairID",       "homologPair",
            "pred_CC",         "len_Qry1",
            "len_Qry2",        "len_H1",
            "len_H2",          "numInt_Q1",
            "numInt_Q2",       "numInt_Homolog1",
            "numInt_Homolog2", "EVal1",
            "EVal2",           "SID1",
            "SID2",            "Positives1",
            "Positives2",      "Similarity1",
            "Similarity2",     "hhpredProb1",
            "hhpredProb2",     "hhpred_pVal1",
            "hhpred_pVal2",    "hhpred_SS1",
            "hhpred_SS2",      "LAL_Q1",
            "LAL_H1",          "LAL_Q2",
            "LAL_H2",          "frac_LAL1",
            "frac_LAL2",       "TP_1",
            "TN_1",            "FP_1",
            "FN_1",            "CC_1",
            "specificity_1",   "sensitivity_1",
            "TP_2",            "TN_2",
            "FP_2",            "FN_2",
            "CC_2",            "specificity_2",
            "sensitivity_2"
        );
    }
    elsif ( $ans == 0 ) {
        @variables_tmp = (
            "qryPairID",       "homologPair", "pred_CC",
            "len_Qry1",        "len_Qry2",
            "len_H1",          "len_H2",
            "numInt_Q1",       "numInt_Q2",
            "numInt_Homolog1", "numInt_Homolog2",
            "EVal1",           "EVal2",
            "SID1",            "SID2",
            "Positives1",      "Positives2",
            "Similarity1",     "Similarity2",
            "hhpredProb1",     "hhpredProb2",
            "hhpred_pVal1",    "hhpred_pVal2",
            "hhpred_SS1",      "hhpred_SS2",
            "LAL_Q1",          "LAL_H1",
            "LAL_Q2",          "LAL_H2",
            "frac_LAL1",       "frac_LAL2",
        );
        my @tmp = ('') x 14;
        push (@variables_tmp, @tmp);

    }
    my $variables_line = join( "\t", @variables_tmp );

    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");

    print OUTPUT "#numInt_Homolog*: # of interface residues of the aligned part of the template.\n";
    print OUTPUT
"#Frac_LAL1: (LAL_Q1/len_Qry1) * (LAL_H1/len_Homolog1), where LAL is local alignment length without gaps\n";
    print OUTPUT
"#Frac_LAL2: (LAL_Q2/len_Qry2) * (LAL_H2/len_Homolog2), where LAL is without gaps\n";
    print OUTPUT "#Format:\n";
    print OUTPUT "$variables_line\n";
    close(OUTPUT);

}

sub readLocalAlignFL {

# INPUT file:
#
# >4lqwD|D:A,GroupID 85,Alignment_Length 136,Sbjct_Length 146,Sbjct_RESI 11 - 146,Query_RESI 1 - 136,Score = 269.4 bits, Expect = 1.5E-39, Identities = 98.5%, Positives = 100.0%, Similarity = 4.4, hhpred = 100.0
# VHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMY
# ?00000000000000000000000000000000000000000000000000000000000000000000000000011111111000000000000000000000000000000000000000000000000000
#
# ---------------------------------------------------------------------
#
# if the int of the cuttent homolog are all '?'s
# the information of this homolog will not be returned,
# so that it will not be used to make predictions.
#
# ---------------------------------------------------------------------

    my $alignFL = shift @_;    #$localAlignFL_P

    my @homologs;
    my $homolog;
    my $queryLength;
    my $num  = 0;    #the number of alignments in the local alignment file
    my $flag = 0;

    #hash refs
    my $localAlignResults;
    my $groupID;
    my $alignmentLength;
    my $sbjctLength;
    my $start_S;
    my $end_S;
    my $start_Q;
    my $end_Q;
    my $BitScore;
    my $EVal;
    my $Identities;
    my $Positives;
    my $hhpredProb;
    my $seq;
    my $int;

    open( ALIGNFL, "<$alignFL" ) || die("Cannot open $alignFL:$!");

    while (<ALIGNFL>) {

        s/[\n\r]//mg;

        if (/Query:\s+(\d+)\s+letters/) {
            $queryLength = $1;
            next;
        }

        if (/^>[\w\.]+|[\w\.]+:[\w\.]+,.+/) {

#>4lqwD|D:A,GroupID 85,Alignment_Length 136,Sbjct_Length 146,Sbjct_RESI 11 - 146,Query_RESI 1 - 136,Score = 269.4 bits, Expect = 1.5E-39, Identities = 98.5%, Positives = 100.0%, Similarity = 4.4,  hhpred_prob = 100.0, hhpred_P-Val = 4.4E-45, hhpred_SS = 26.3
            $_ =~ s/>//;

            #----------------------------------------------
            #-- check format
            my @tmp = split( ',', $_ );

            if ( scalar @tmp ne 14 ) {
                my $num_field = scalar @tmp;
                die(
"Format of $alignFL is wrong (the header has $num_field fields). Current line: $_:$!"
                );
            }

            #-- checking format done.
            #----------------------------------------------

            $num++;
            $flag = 1;
            my $pattern = "a-zA-Z0-9_\.\-";

            $homolog = shift @tmp;    #$homolog = '1o96Q|Q:Z'

            ($groupID) = $_ =~ m/GroupID\s+(\d+)/i
              ; #sometimes one protein aligns with another protein at two different places. So use $num to differentiate them

            ( $localAlignResults->{'LAL'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Alignment_Length\s+(\d+)/i;
            ( $localAlignResults->{'sbjctLen'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Sbjct_Length\s+(\d+)/i;
            ( $localAlignResults->{'start_S'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Sbjct_RESI\s+(\d+)\s*-\s*\d+/i;
            ( $localAlignResults->{'end_S'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Sbjct_RESI\s+\d+\s*-\s*(\d+)/i;
            ( $localAlignResults->{'start_Q'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Query_RESI\s+(\d+)\s*\-\s*\d+/i;
            ( $localAlignResults->{'end_Q'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Query_RESI\s+\d+\s*\-\s*(\d+)/i;
            ( $localAlignResults->{'BitScore'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Score\s*=\s*([$pattern]+)\s*bits/i;
            ( $localAlignResults->{'EVal'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Expect\s*=\s*([$pattern]+)/i;
            ( $localAlignResults->{'Identities'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Identities\s*=\s*([$pattern]+)\s*%/i;
            ( $localAlignResults->{'Positives'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Positives\s*=\s*([$pattern]+)/i;
            ( $localAlignResults->{'Similarity'}->{"$homolog,$groupID"} ) =
              $_ =~ m/Similarity\s*=\s*([$pattern]+)/i;
            ( $localAlignResults->{'hhpredProb'}->{"$homolog,$groupID"} ) =
              $_ =~ m/hhpred_prob\s*=\s*([$pattern]+)/i;
            ( $localAlignResults->{'hhpred_Pval'}->{"$homolog,$groupID"} ) =
              $_ =~ m/hhpred_P-Val\s*=\s*([$pattern]+)/i;
            ( $localAlignResults->{'hhpred_SS'}->{"$homolog,$groupID"} ) =
              $_ =~ m/hhpred_SS\s*=\s*([$pattern]+)/i;
            next;

        }

        if ( $flag == 1 && /^([A-Za-z\-]+)$/ ) {
            $localAlignResults->{'QrySeqAligned'}->{"$homolog,$groupID"} = $1;
            next;
        }

        if ( $flag == 1 && /^([01\?\-]+)$/ ) {

            my $int_tmp = $1;

            if ( $int_tmp =~ /^[\?]+$/ ) {

                #if the int of the cuttent homolog are all '?'s
                #the information of this homolog will not be returned,
                #so that it will not be used to make predictions.
                #
                foreach my $key ( keys %{$localAlignResults} ) {
                    delete $localAlignResults->{$key}->{"$homolog,$groupID"};
                }
            }

            else {

                $localAlignResults->{'HomologIntAligned'}->{"$homolog,$groupID"}
                  = $int_tmp;
                $flag = 0;
            }
            next;

        }

    }
    close ALIGNFL;

    if ( $num == 0 ) {
        die("Nothing read from $alignFL:$!");
    }

    #sort @homologs based on their hhpred prob. score
    my @homologs_sorted =
      sort {
        $localAlignResults->{'hhpredProb'}->{$b}
          <=> $localAlignResults->{'hhpredProb'}->{$a}
      } keys %{ $localAlignResults->{'hhpredProb'} };

    if (!defined $queryLength ){
        die("query length not read from $alignFL:$!");
    }

    return ( \@homologs_sorted, $queryLength, $localAlignResults );
}

sub readSeqIntFL {
    my $SeqIntFL = shift @_;
    my $comments ='';
    my $header;
    my %seqs;
    my %int;

    open( seq_int_FL, "<$SeqIntFL" )
      || die("Cannot open $SeqIntFL:$!");
    while (<seq_int_FL>) {
        s/[\n\r]//mg;

        if (/^#/){
            $comments = $comments."$_\n";
        }
        #        if (/^>(\w{4})(:|)(\w{1})(\|\w{0,}:\w{0,}|)$/) {
        if (/^>(.+)$/) {

            $header = $1;
            next;
        }

        if (/^[A-Za-z]+$/) {
            my $seq = $_;
            $seqs{$header} = $seq;
            next;
        }
        if (/^[01\?\-]+$/) {
            $int{$header} = $_;
            $header = '';
            next;
        }

    }
    close(seq_int_FL);

    if ( !%seqs || !%int ) {
        die("Nothing read from $SeqIntFL:$!");
    }

    return ( $comments, \%seqs, \%int );

}


sub readSeqIntFL_old {
    my $SeqIntFL = shift @_;
    my $header;
    my %seqs;
    my %int;

    open( seq_int_FL, "<$SeqIntFL" )
      || die("Cannot open $SeqIntFL:$!");
    while (<seq_int_FL>) {
        s/[\n\r]//mg;

        #        if (/^>(\w{4})(:|)(\w{1})(\|\w{0,}:\w{0,}|)$/) {
        if (/^>(.+)$/) {

            $header = $1;
            next;
        }

        if (/^[A-Za-z]+$/) {
            my $seq = $_;
            $seqs{$header} = $seq;
            next;
        }
        if (/^[01\?\-]+$/) {
            $int{$header} = $_;
            $header = '';
            next;
        }

    }
    close(seq_int_FL);

    if ( !%seqs || !%int ) {
        die("Nothing read from $SeqIntFL:$!");
    }

    return ( \%seqs, \%int );

}

sub getIntNum {

    #    >4dgaA|A:C
    #    MVNPTVFFDIAVDGEPLGRVSFELFADKVPKTAENFRALST
    #    ?00000000000000000????0000000000???000000
    #    >4dgaC|C:A
    #    MPIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPE
    #    ?000???????0000000000000000000000000
    #
    my $intFL = shift @_;
    my $num_residue1;    # the whole seq
    my $num_residue2;    # the structural seq
    my $num_int;
    my $header;
    my $format_checker = 0;

    open( INTFL, "<$intFL" ) || die("Cannot open $intFL:$!");
    while (<INTFL>) {
        s/[\n\r]//mg;

        if (/^>(.+)/) {
            $header                  = $1;
            $num_int->{$header}      = 0;
            $num_residue1->{$header} = 0;
            $num_residue2->{$header} = 0;

            if ( $format_checker != 0 && $format_checker != 2 ) {
                die(
"Format is wrong : $intFL (format_check = $format_checker :$!"
                );
            }

            $format_checker = 0;    #reset format_checker
            next;
        }
        if (/^[A-Za-z]+$/) {
            $num_residue1->{$header} = length($_);
            $format_checker++;
            next;
        }
        if (/^[01?\-]+$/) {
            $format_checker++;

            if (/^[\?]+$/) {

                #-- the interface line is all ?s
                $num_int->{$header}      = 'Nan';
                $num_residue2->{$header} = length($_);
                next;
            }

            my @int = split( //, $_ );
            foreach (@int) {
                if (/1/) {
                    $num_int->{$header}++;
                    $num_residue2->{$header}++;
                }
                elsif (/[0\?\-]/) {
                    $num_residue2->{$header}++;
                }
                else {
                    die(
"IntFL format error. Interface can only be 0, 1, ?, -. Check $intFL:$!"
                    );
                }
            }

        }
    }
    close(INTFL);

    if ( !defined $num_int || !defined $num_residue1 ) {
        die("num_int or num_residues not read from $intFL:$!");
    }

    return ( $num_residue1, $num_residue2, $num_int );

}

sub compare_qryInt_homologInt {

    #-- one homolog may have different alignment with the same query

    my $queryLength    = shift @_;
    my $aligned_QrySeq = shift @_
      ; # there may be dashes. $aligned_QrySeq->{$homolog,$groupID}='SEFEQEFQE----QEFEFE'
    my $int_Q      = shift @_;    #query interface on the whole qry seq
    my $homologIDs = shift @_;    # $homologIDs = ( "$homolog1,$groupID", ...)
    my $int_homologs =
      shift @_;    # $int_homologs->{"homolog,groupID"} = '0101010000'
    my $start_Q = shift @_;    # $start_Q->{"$homolog,$groupID"} = 34
    my $end_Q   = shift @_;

    #--

    my ( $TP, $TN, $FP, $FN, $CC, $specificity, $sensitivity, $accuracy );

    foreach my $homologID ( @{$homologIDs} ) {

        #-- one homolog may have different alignment with the same query
        # $homologID = '2o96Q|Q:Z,45'

#--- get the interface_H that is aligned with the query (3 steps)
#1) process query sequence line. Find the indices of '-' in the query sequence that is aligned with a homolog
        my @where = &findDashInSeq( $aligned_QrySeq->{$homologID} );

#2).process homolog interface line. Remove the remove the interface sign that corresponds to -, and add ? to the two ends

        my $int_H = join(
            '',
            @{
                &modifyIntLine(
                    $int_homologs->{$homologID}, $queryLength,
                    $start_Q->{$homologID},      $end_Q->{$homologID},
                    @where
                )
            }
        );    #$int_H = '00011000000'

        if ( $int_Q !~ /^\?+$/ ) {

#------ compare interfaces of query and homolog, and calculate TP, TN, FP, FN, CC, etc-------

            #            print "int_Q: $int_Q\n";
            #            print "int_H: $int_H\n";

            (
                $TP->{$homologID},          $TN->{$homologID},
                $FP->{$homologID},          $FN->{$homologID},
                $CC->{$homologID},          $specificity->{$homologID},
                $sensitivity->{$homologID}, $accuracy->{$homologID}
            ) = &TP_TN_FN_FP_CC_specif_sensi_acc( $int_Q, $int_H );

            #- both $int_Q and $int_H have to be strings: '000110'
        }
        else {

            #- if the interfaces are all questions marks
            (
                $TP->{$homologID},          $TN->{$homologID},
                $FP->{$homologID},          $FN->{$homologID},
                $CC->{$homologID},          $specificity->{$homologID},
                $sensitivity->{$homologID}, $accuracy->{$homologID}
            ) = ('') x 8;
        }

    }
    return ( $TP, $TN, $FP, $FN, $CC, $specificity, $sensitivity, $accuracy );
}

sub countInt {
    my $int = shift @_;
    my $num = 0;
    if ( !defined $int ) {
        die("interface is not defined:$!");
    }

    $int =~ s/[\n\r\s]//g;
    if ( $int =~ /^\?+$/ ) {

        #-- the whole int is ?
        $num = 'nan';
        return $num;
    }

    my @tmp = split( //, $int );
    foreach (@tmp) {
        if (/1/) {
            $num++;
        }
        elsif (/[0\?\-]+/) {
            next;
        }
        else {
            die(
"interface can only be 0, 1 , - or ?. The current interface: $int :$!"
            );
        }
    }
    return $num;

}

sub readDeleteFL {

    # DeleteFL Format:
    #	1a2kA:1a2kC =>
    #	1a14H:1a14N => 1NMBH:1NMBN,1A14H:1A14N
    #   A:B => 1NMBH:1NMBN,1NMAH:1NMAN
    #	A:B => 3fclB:*,1lqmA:*
    #	A:B => 3fcl*
    #	A:B => 3fcl*:3fcl*
    #	A:B => 3fclA:3fcl*
    #

    my $deleteFL = shift @_;
    my %homoInterologsToBeDelMap;
    my @homologs_tmp2;

    open( INPUT, "<$deleteFL" ) || die("Cannot open $deleteFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        s/\s+//g;

        if (/^\w+:\w+\s*=>/) {

            #	1a14H:1a14N => 1NMBH:1NMBN,1A14H:1A14N

            my ( $qryPair, $homologs ) = split( /=>/, $_ );
            $qryPair =~ s/\s+//g;     #1a14H:1a14N
            $homologs =~ s/\s+//g;    #1cgjE:1cgjI,1cgiE:1cgiI

            #---------
            #change pdb IDs to lower cases
            my @homologPairs = split( /,/, $homologs );
            my @homologpairs_new;
            foreach my $homologPair (@homologPairs) {

                #	1cgjE:1cgjI
                #   1lqmA:*
                #   1lqm*
                #   1lqmA:1lqm*

                my ( $prot1, $prot2 ) = split( /:/, $homologPair );

                #-- convert pdb ID to lower case
                $prot1 = &format_homologsToBeDel($prot1);
                if ( defined $prot2 ) {
                    $prot2 = &format_homologsToBeDel($prot2);
                }

                #-- double check
                if ( defined $prot2
                    && substr( $prot1, 0, 4 ) ne substr( $prot2, 0, 4 ) )
                {
                    print LOG
"WARNING: the homolog-pair-to-be-del $homologPair does not have the same pdb ID. It is discarded.\n";
                    next;
                }

                #--
                if ( defined $prot2 ) {
                    push @homologpairs_new, "$prot1:$prot2";
                }
                else {
                    # $homologPair = '1ahj*'
                    push @homologpairs_new, "$prot1";
                }
            }

            #---------
            $homoInterologsToBeDelMap{$qryPair} = \@homologpairs_new;
        }
    }
    close INPUT;

    if ( !%homoInterologsToBeDelMap ) {
        print(
"\n\n\nWARNING: No homologs to be deleted defined in $deleteFL. Check $deleteFL.\n\n\n"
        );
    }

    return \%homoInterologsToBeDelMap;

}

sub format_homologsToBeDel {

    #-- this function is called by &readDeleteFL{}
    #
    #-- For a homolog to be deleted, convert pdb ID of to lower case
    #
    #    1AHJA to 1ahjA
    #    * is kept not changed

    my $homolog = shift @_;    #1AHJA or *
    my $homolog_new;
    if ( length($homolog) == 5 ) {

        my $pdbID = lc( substr( $homolog, 0, 4 ) );
        my $chnID = substr( $homolog, 4, 1 );
        $homolog_new = "$pdbID$chnID";
    }
    elsif ( $homolog eq '*' ) {
        $homolog_new = $homolog;
    }
    else {
        die(
"Delete file format wrong (homolog-to-be-del = $homolog). The homolog has to be 5-letter long or is a wild card *:$!"
        );
    }

    #    print "homolog_ori: $homolog, homolog_new: $homolog_new\n\n";
    return $homolog_new;

}

#-----
sub Write_header_LocalAlignFL {
    our $serverName;    # 'PSHOMPPIv2.0';

    my $outputFL = shift @_;
    unlink($outputFL) if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    print OUTPUT
"# =========================================================================================================\n";
    print OUTPUT "# Auto generated file by $serverName\n";
    print OUTPUT "# Format description:\n";
    print OUTPUT
      "# >homolog's PDBIDchainID|ChainID:interacting ChainID|GroupID|....\n";
    print OUTPUT
"# Aligned part of QUERY Seq (Amino acid sequence extracted from new PPIDB).\n";
    print OUTPUT
"# Aligned part of HOMOLOG interface. Interface of the protein chain with the interacting chain: 1s denote interface residues, 0s denote non-interface residue\n";
    print OUTPUT
"# =========================================================================================================\n";
    close(OUTPUT);
}

sub subSbjctInt {

#to extract substring of interfaces of homologs, given the starting pos and ending position

    my ( $SeqIntFL, $homologPartnerSpc, $start, $end ) = @_;

    my $interface;
    my $flag = 0;
    open( SEQINTFL, "<$SeqIntFL" )
      || die("Can not open $SeqIntFL!\n");
    my ( $homolog, $chainID, $interactingChainID ) =
      split( /\||:/, $homologPartnerSpc );

    while (<SEQINTFL>) {
        s/[\n\r]//mg;
        if (/^>$homolog\|$chainID:$interactingChainID$/) {

            $flag = 1;
            next;
        }

        if ( $flag == 1 && /^[01\?\-]+$/ ) {
            $interface = substr( $_, $start - 1, $end - $start + 1 );
            last;
        }
    }
    close(SEQINTFL);

    if ( !defined $interface ) {
        print
"WARNING: Interface of homolog ($homologPartnerSpc) does not exist in $SeqIntFL (start $start - end $end)\n";
    }
    return $interface;
}

sub insertSpc {
    my ( $localAlign, $interface ) = @_;
    my $i = 0;                          #index of spaces in $localAlign
    my @temp = split '', $localAlign;
    foreach (@temp) {
        if (/\-/) {
            substr( $interface, $i, 0 ) = '-';
        }
        $i += 1;
    }
    return $interface;
}

sub rmSpc {

    my ( $localAlign, $interface ) = @_;
    my @interface_tmp = split( //, $interface );

    my @temp = split '', $localAlign;

    my $i = 0;
    foreach (@temp) {
        if (/\-/) {
            splice @interface_tmp, $i, 1;
        }
        $i = $i - 1;
    }

    $interface = join( '', @interface_tmp );

    return $interface;

}

sub findDashInSeq {
    my $seq = shift @_;

    my @where;

    #find the index for -

    my $index = index( $seq, '-' );

    while ( $index >= 0 ) {

        if ( $index >= 0 ) {
            push @where, $index;
        }
        $index = index( $seq, '-', $index + 1 );
    }
    return @where;
}

sub modifyIntLine {
    my $int         = shift @_;
    my $queryLength = shift @_;
    my $start_Q     = shift @_;
    my $end_Q       = shift @_;
    my @where       = @_;

    my @new_int;

    if ( !defined $int ) {
        die("Input interface not defined:$!");
    }

    #1. remove the interface sign that corresponds to -
    my $i = 0;
    if (@where) {
        foreach (<@where>) {
            substr( $int, $_ - $i, 1 ) = '';
            $i++;
        }
    }

    #2. add ? to the two ends
    my $multiX1 = "?" x ( $start_Q - 1 );            #$num_multiX question marks
    my $multiX2 = "?" x ( $queryLength - $end_Q );
    $int = "$multiX1$int$multiX2";

    if ( !defined $int ) {
        die("Output interface not defined:$!");
    }

    @new_int = split( //, $int );

    return \@new_int;
}

sub TP_TN_FN_FP_CC_specif_sensi_acc {

#----------compare interfaces of query and homolog, and calculate TP, TN, FP, FN, CC, etc-----------#

    my @int_Q = split( //, shift @_ );    #query interface as array ref
    my @int_H = split( //, shift @_ );    #homolog interface as array ref

    if ( scalar @int_Q ne scalar @int_H ) {
        die(
"Error: int_Q and int_H have different length.\nint_Q:@int_Q.\nint_H:@int_H:$!"
        );
    }

    my ( $TP, $TN, $FP, $FN, $CC, $specificity, $sensitivity, $accuracy ) =
      ( 0, 0, 0, 0, 0, 0, 0, 0 );

    for my $i ( 0 .. ( scalar @int_Q - 1 ) ) {
        if ( $int_Q[$i] =~ /1/ ) {
            if ( $int_H[$i] =~ /1/ ) {
                $TP++;
            }
            elsif ( $int_H[$i] =~ /0/ ) {
                $FN++;
            }
        }
        elsif ( $int_Q[$i] =~ /0/ ) {
            if ( $int_H[$i] =~ /1/ ) {
                $FP++;
            }
            elsif ( $int_H[$i] =~ /0/ ) {
                $TN++;
            }
        }

    }

    #------------calculate CC etc.------------#
    if ( ( $TP + $FN ) * ( $TP + $FP ) * ( $TN + $FP ) * ( $TN + $FN ) != 0 ) {
        $CC =
          ( $TP * $TN - $FP * $FN ) /
          sqrt( ( $TP + $FN ) * ( $TP + $FP ) * ( $TN + $FP ) * ( $TN + $FN ) );
    }
    if ( $TP + $FP != 0 ) {
        $specificity = $TP / ( $TP + $FP );
    }
    if ( $TP + $FN != 0 ) {
        $sensitivity = $TP / ( $TP + $FN );
    }
    if ( $TP + $TN + $FN + $FP != 0 ) {
        $accuracy = ( $TP + $TN ) / ( $TP + $TN + $FN + $FP );
    }

    #------
    $CC          = sprintf( "%.2f", $CC );
    $specificity = sprintf( "%.2f", $specificity );
    $sensitivity = sprintf( "%.2f", $sensitivity );
    $accuracy    = sprintf( "%.2f", $accuracy );

    return ( $TP, $TN, $FP, $FN, $CC, $specificity, $sensitivity, $accuracy );
}

sub writePredictionFL_noHomologFound {

#output file:
#                >A|No homolog found
#                GDKPIWEQIGSSFIQHYYQLFDNDRTQLGAIYIDASCLTWEGQQFQGKAAIVEKLSSLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKADEDPIMGFHQMFLLKNINDAWVCTNDMFRLALHNF
#Pred_score:     ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
#Prediction:     ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    use strict;

    my $jobDIR = shift @_;

    opendir( DIR, $jobDIR );
    my @taskDIRs = grep { -d "$jobDIR/$_" && /^\w+:\w+$/ } readdir(DIR);
    closedir(DIR);

    if ( scalar @taskDIRs == 0 ) {
        die("No taskDIRs found under $jobDIR:$!");
    }

    foreach my $pairDIR (@taskDIRs) {

        #$pairDIR = 'A:B'

        my @taskIDs = split( /:/, $pairDIR );

        for ( my $i = 0 ; $i < scalar @taskIDs ; $i++ ) {

            #file variables
            my $taskID  = $taskIDs[$i];
            my $taskDIR = "$jobDIR/$pairDIR/$taskID";
            my $fastaFL = "$taskDIR/$taskID.fasta.txt";
            my $intFL   = "$taskDIR/$taskID.int";         #interface file
            my $localAlignFL =
              "$taskDIR/seq_interface_homologs.localAlign"
              ;                                           #input file of KNN.pl

            my $predictionResultsDIR = "$taskDIR/predictionResults";
            my $outputFL = "$predictionResultsDIR/$taskID.KNN.prediction";

            mkdir $predictionResultsDIR unless ( -d $predictionResultsDIR );
            unlink $outputFL if ( -e $outputFL );

            #--
            my $partner;
            if ( $i == 0 ) {
                $partner = $taskIDs[1];

            }
            else {
                $partner = $taskIDs[0];
            }
            my $queryID = "$taskID\\|$taskID\\:$partner";    #'B|B:A'

            #--

            #            my $seq        = &readSeqFL($fastaFL);
            my $seq        = &readFastaFL($fastaFL);
            my $predScore  = join( ',', ('?') x length($seq) );
            my $prediction = '?' x length($seq);

            my $actualInt;

            if ( -e $intFL ) {

                ( my $seq_Q, my $actInt_ref, my $length_Q, my $numInt_Q ) =
                  &readINTFL($intFL);
                $actualInt = join( '', @$actInt_ref );
            }
            else {
                $actualInt = '';
            }

            my $Mode = "No homolog found";

            &writePredictionFL( $Mode, $taskID, $seq, $predScore, $prediction,
                $actualInt, $outputFL );

        }

    }

}

sub safeTwiDark {
    use strict;

    my $jobDIR         = shift @_;
    my $jobID          = shift @_;
    my $fullQryPairLst = shift @_;
    my $FullStatFL     = shift @_;

    #variables
    my $qryPredModeFL = "$jobDIR/predictionMode.txt";
    unlink $qryPredModeFL if ( -e $qryPredModeFL );

#-------------------------
#-- read stat file. The pred_CC column will be used to determine which homolog pairs belong to which zone.

    my ( $templateStats, $headers ) = &readStatisticsFL_all($FullStatFL);

    #--  $templateStats->{$qryPairID}->{$templateID}->{$header} = $stats[$i];

    #-------------------------
    # ----- Safe Zone

    print LOG "\n\n-------- start predicting...\n\n";

    my $Mode = 'SafeMode';

    #predict qrys under $Mode

    my $QRYpairs = &readTestChainPairs($fullQryPairLst);    #'A:B'

    my ( $qryPairWithSafeZoneHomolog, $qryPairCannotbePredicted_inSafe ) =
      &mode( $jobDIR, $QRYpairs, $templateStats, $Mode, $qryPredModeFL );

    my $num_qryWithSafeZoneHomolog = scalar @$qryPairWithSafeZoneHomolog;
    print LOG
"There are $num_qryWithSafeZoneHomolog qry pairs with Safe Zone Homologs. \n";

    #-------------------------
    # ----- Twilight Zone

    $Mode = 'TwilightMode1';

    my $num = scalar(@$qryPairCannotbePredicted_inSafe);
    my ( $qryPairWithTwiHomolog1, $qryPairCannotbePredicted_inTwi1 );

    if ( $num > 0 ) {

        ( $qryPairWithTwiHomolog1, $qryPairCannotbePredicted_inTwi1 ) =
          &mode( $jobDIR, $qryPairCannotbePredicted_inSafe,
            $templateStats, $Mode, $qryPredModeFL );

        my $num_qryWithTwiHom1 = scalar @$qryPairWithTwiHomolog1;
        print LOG
"There are $num_qryWithTwiHom1 qry pairs with only Twilight Zone 1 (or worse) Homologs.\n";

        #----twi2
        $Mode = 'TwilightMode2';
        my $num = scalar(@$qryPairCannotbePredicted_inTwi1);
        my ( $qryPairWithTwiHomolog2, $qryPairCannotbePredicted_inTwi2 );

        if ( $num > 0 ) {

            ( $qryPairWithTwiHomolog2, $qryPairCannotbePredicted_inTwi2 ) =
              &mode( $jobDIR, $qryPairCannotbePredicted_inTwi1,
                $templateStats, $Mode, $qryPredModeFL );

            my $num_qryWithTwiHom2 = scalar @$qryPairWithTwiHomolog2;
            print LOG
"There are $num_qryWithTwiHom2 qry pairs with only Twilight Zone 2 (or worse) Homologs.\n";

            #-------------------------
            # ----- Dark Zone

            $Mode = 'DarkMode';

            my $num = scalar(@$qryPairCannotbePredicted_inTwi2);
            my ( $qryPairWithDkZoneHomolog, $qryPairCannotbePredicted_inDk );

            if ( $num > 0 ) {

                ( $qryPairWithDkZoneHomolog, $qryPairCannotbePredicted_inDk ) =
                  &mode( $jobDIR, $qryPairCannotbePredicted_inTwi2,
                    $templateStats, $Mode, $qryPredModeFL );

                my $num_qryWithDkZoneHomolog =
                  scalar @$qryPairWithDkZoneHomolog;
                print LOG
"There are $num_qryWithDkZoneHomolog qry pairs with only Dark Zone Homologs.\n";
            }
        }
    }

}

sub cluster_PS_resiPair {

    #-- cluster templates and predict Ca-Ca based on clustered templates
    #
    our $logFL;
    our $atomDistThr;    # interface definition: CA-CA distance cutoff

    my $CaDistThr = $atomDistThr;

    #---
    my ( $jobDIR, $pairlstFL ) = @_;
    my $flag_prediction = 1;    #0: predictions made. 1: No predictions made.
    my $num             = 0;    #number of qry pairs that have Ca-Ca predicted

    my @QRYpairs = @{ &readTestChainPairs($pairlstFL) };    #'A:B' #xue-fix

    foreach my $queryPairID (@QRYpairs) {

        # $queryPairID = 'A:B'

        print LOG
"\n\n*** Predicting interacting residue pairs for $queryPairID ***\n\n";

        #        my ( $A, $B ) = split( /:/, $queryPairID );
        my ( $pairDIR, $A, $B ) = &getPairDIR2( $queryPairID, $jobDIR );

        if ( !-d $pairDIR ) {
            die("pairDIR:$pairDIR does not exist:$!");
        }

        #--
        my $command =
          "perl mainFun_superimpose_caca.pl $CaDistThr $pairDIR >> $logFL";
        print LOG "COMMAND: $command\n\n";
        system($command) == 0 or die("FAILED: $command:$!");

        #-- check whether there are Ca-Ca files generated for $A:$B
        my $Ca_CaDIR = "$pairDIR/superimposed_models/Ca_Ca_distances";
        my @ca_caFLs = `ls $Ca_CaDIR/cluster*_Ca_Ca_distance.*txt`;

        if ( !@ca_caFLs ) {
            print LOG
"No Ca-Ca files generated under : $Ca_CaDIR. No predictions made.\n";
        }
        else {
            #Ca-Ca predictions made for this qry pair
            $num++;
        }

    }

    print LOG "$num clusters of templates are generated for $jobDIR\n";

    if ( $num != 0 ) {
        $flag_prediction = 0;    #0: predictions made. 1: No predictions made.
    }

    return $flag_prediction;

}

sub addHeader2CaCaFLs {

    # add header to Ca-Ca files
    # save the Ca-Ca files to $outputDIR

    use strict;
    use File::Basename;

    our $logFL;

    my $jobDIR   = shift @_;
    my @QRYpairs = @{ shift @_ };    #'A:B' , '1fevA:1fevB'
    my ( $startTime, $endTime, $timeUsed ) = @_;
    my $outputDIR = "$jobDIR/Ca_Ca_distances";

    mkdir $outputDIR if ( !-d $outputDIR );

    my $taskID;
    my $taskID_partner;
    our $safe_filename_characters;

    print LOG "\n\n\n\n Add header to Ca-Ca files...\n";

    #process each task ID

    foreach my $QRYpair (@QRYpairs) {

        #'A:B'

        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $jobDIR );

        #1. write A|A:B
        $taskID         = $c1;
        $taskID_partner = $c2;
        my $Ca_CaDIR = "$pairDIR/superimposed_models/Ca_Ca_distances";
        my @ca_caFLs = `ls $Ca_CaDIR/cluster*_Ca_Ca_distance.*txt`;

        if ( !@ca_caFLs ) {
            print LOG
"No Ca-Ca files found under : $Ca_CaDIR. So no predictions made for $QRYpair\n";
            next;
        }

        #--get prediction zone for $QRYpair
        my $predictFL =
          "$pairDIR/$taskID/predictionResults/$taskID.KNN.prediction";
        my $zone = &getZone($predictFL);

        foreach my $ca_caFL_ori (@ca_caFLs) {
            $ca_caFL_ori =~ s/[\n\r]//mg;
            my $basename       = basename($ca_caFL_ori);
            my $outputDIR_pair = "$outputDIR/$c1:$c2";
            mkdir($outputDIR_pair) if ( !-d $outputDIR_pair );
            my $outputFL = "$outputDIR_pair/$basename";
            &finalPredictionFLheader( $outputFL, $startTime, $endTime,
                $timeUsed );
            &writeFinalCaCaPredictionFL( $outputFL, $ca_caFL_ori, $zone );


            #-- copy tbl and pml to $outputDIR
            my @tbl_pml_files = glob ("$Ca_CaDIR/*tbl $Ca_CaDIR/*pml");
            foreach (@tbl_pml_files){
               copy($_, $outputDIR_pair) or die "Copy failed: $!";;
            }


            print LOG
"All predictions for proteins in $jobDIR are put into $outputFL(usr-friendly format).\n ";
        }

    }

}

sub getZone {

    #    INPUT:
    #     >A|SafeMode
    #VNPTVFFDIAVDGEPLGRVSF

    my $predictFL = shift @_;
    my ( $seq, $predictionScore, $pInt, $zone ) = &ReadPredictionFL($predictFL);

    if ( !defined $zone ) {
        die("No zone info read from $predictFL:$!");
    }
    return $zone;

}

sub writeFinalCaCaPredictionFL {

    #INPUT:
    #    A   B   A_template1 B_template1 mean    min max
    #    103 78  104 88  6.955   5.555   7.760
    #    103 77  104 87  8.074   6.741   9.018

    my $outputFL    = shift @_;
    my $ca_caFL_ori = shift @_;
    my $zone        = shift @_;

    my $command = "echo \" Mode = $zone \" >> $outputFL";
    system($command ) == 0 or die("FAILED: $command:$!");
    $command =
"egrep \"^[^#!]\" $ca_caFL_ori | awk \'{print \$1 \" \" \$2 \" \" \$3 \" \" \$4 \" \" \$9 \" \" \$10 \" \" \$11}\' >>$outputFL";
    system($command ) == 0 or die("FAILED: $command:$!");
    print LOG "$outputFL generated\n";

}

sub readStatisticsFL_all {

#Input $statisticsFL example:
#
#qryPairID homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2 EVal1  EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL1    LAL2    frac_LAL1 frac_LAL2   TP_1    TN_1    FP_1    FN_1    CC_1    specificity_1sensitivity_1  TP_2    TN_2    FP_2    FN_2    CC_2 specificity_2   sensitivity_2
#Q:P 2x2dC,42:2x2dD,69   164 137 165 147 Nan Nan 14  8   1.1E-45 3E-41   99.3    98.599.3    100.0   4.7 4.4 100.0   100.0   164    136 0.99    0.92

    use strict;

    my $statFL = shift @_;
    my $templateStats;
    my @headers;
    my @stats;

    open( INPUT, "<$statFL" ) or die("Cannot open $statFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^qryPairID.+$/) {

#qryPairID homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2 EVal1  EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL1    LAL2    frac_LAL1 frac_LAL2   TP_1    TN_1    FP_1    FN_1    CC_1    specificity_1sensitivity_1  TP_2    TN_2    FP_2    FN_2    CC_2 specificity_2   sensitivity_2
            @headers = split( /\s+/, $_ );
            shift @headers;
            shift @headers;
            next;
        }
        if (/^\w+:\w+\s+/) {

#Q:P 2x2dC,42:2x2dD,69   164 137 165 147 Nan Nan 14  8   1.1E-45 3E-41   99.3    98.599.3    100.0   4.7 4.4 100.0   100.0   164    136 0.99    0.92
            @stats = split( /\s+/, $_ );
            my $qryPairID  = shift @stats;
            my $templateID = shift @stats;

            if ( scalar @stats ne scalar @headers ) {
                my $num_header = scalar @headers;
                my $num_stat   = scalar @stats;
                die(
"ERROR! Current line: $_\nstatFL's header field number ($num_header) is different from the value fields ($num_stat). Check statFL $statFL:$!"
                );
            }

            for ( my $i = 0 ; $i < scalar @headers ; $i++ ) {
                my $header = $headers[$i];
                $templateStats->{$qryPairID}->{$templateID}->{$header} =
                  $stats[$i];
            }

        }
    }
    close INPUT;

    if ( !defined $templateStats ) {
        die("Nothing read from $statFL:$!");
    }

    return ( $templateStats, \@headers );

}

sub mode {

    use strict;

    my ( $jobDIR, $QRYpairs, $templateStats, $Mode, $qryPredModeFL ) = @_;


    print LOG "\n\n------------------------------\n";
    print LOG "\t***   $Mode   ***\n";
    print LOG "------------------------------\n\n";


    my ( $qryPairCanbePredicted, $qryPairCannotbePredicted ) =
      &batch( $jobDIR, $QRYpairs, $templateStats, $Mode );    #-- xue-fix

    #write prediction mode file
    &writePredictionModeFL( $Mode, $qryPairCanbePredicted, $qryPredModeFL );

    print LOG "\n$Mode finished.\n\n";

    return ( $qryPairCanbePredicted, $qryPairCannotbePredicted );

}

sub readTestChainPairs {

    #input format 'protA:protB'
    #return Qrypairs =  (A:B, ...) or ("$pdbID$chainID1:$pdbID$chainID2",...)
    #

    my $inputFL = shift @_;
    my @QRYpairs;

    open( INPUT, "<$inputFL" )
      || die("Cannot open $inputFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^(\w+):(\w+)/) {

            #1a14N:1a14H

            my $prot1 = $1;
            my $prot2 = $2;

            my $pair = "$prot2:$prot1";

            #remove redundant pairs
            #'A:B' and 'B:A' are redundant pairs

            my $lastpair = $QRYpairs[ scalar @QRYpairs - 1 ];
            if ( $QRYpairs[ scalar @QRYpairs - 1 ] ) {

                if ( $pair ne $lastpair ) {
                    push @QRYpairs, $pair;
                }
            }
            else {
                push @QRYpairs, "$prot1:$prot2";
            }

        }
    }
    close(INPUT);

    if ( !@QRYpairs ) {
        die("No qry pairs read from $inputFL:$!");
    }

    return ( \@QRYpairs );

}

sub finalPredictionFLheader {

    use strict;
    use globalVariables;
    our $serverName;

    my $finalPredictionOutputFL = shift @_;
    my ( $startTime, $endTime, $timeUsed ) = @_;

    unlink $finalPredictionOutputFL if ( -e $finalPredictionOutputFL );

    open( OUTPUT, ">>$finalPredictionOutputFL" );
    print OUTPUT
"Prediction results of Partner-specific interface residues by PS-HomPPI (http://ailab1.ist.psu.edu/$serverName).\n\n";
    print OUTPUT "\t*****************\n\n";
    print OUTPUT "\tInterface prediction starts at $startTime.\n";
    print OUTPUT "\tInterface prediction ends at $endTime.\n";
    print OUTPUT "\tTotal time used: $timeUsed seconds.\n\n";
    print OUTPUT "\t*****************\n\n";

#	print OUTPUT "\nNote:\n\t1. Chain IDs in the docked models are renamed, because one docked model may have the same chain IDs.\n";
#	print OUTPUT "\tPlease see the map file for the mapping of the Chain IDs in this file and Chain IDs in the original docked models.\n";
#	print OUTPUT "\t2. DockRank only uses the predicted interfaces between recptor and ligands to rank the docked models.\n\n\n";
    print OUTPUT "Notations:\n";
    print OUTPUT
"\t1. A|A:B: the interface residues of protein A that interact with protein B.\n";
    print OUTPUT "\t2. MODE: \n";
    print OUTPUT
"\t\t Mode = SafeMode: the query protein can find homologous interacting pairs in Safe Zone.\n";
    print OUTPUT
"\t\t Mode = TwilightMode1: the query protein can find homologous interacting pairs in Twilight Zone 1.\n";
    print OUTPUT
"\t\t Mode = TwilightMode2: the query protein can find homologous interacting pairs in Twilight Zone 2.\n";
    print OUTPUT
"\t For more details about the Safe/Twilight/Dark Zone, please refer to the paper for PS-HomPPI:\n\tXue, L. C., Dobbs, D., & Honavar, V. (2011). HomPPI: A Class of sequence homology based protein-protein interface prediction methods. BMC Bioinformatics, 12, 244.\n";

#	print OUTPUT "\t3. resSeq: residual sequence number. It is extracted from docked models.\n";
    print OUTPUT
"\t3. pINT: predicted interface residues. 1: interface. 0: non-interface. ?: no prediction can be made.\n";
    print OUTPUT
"\t4. SCORE: prediction score from PS-HomPPI. The higher the score the higher prediction confidence.\n";
    print OUTPUT "\n\n-----------------------------------------\n\n\n";
    close OUTPUT;
}

sub ReadPredictionFL {

    #read the prediction FL of PS-HomPPI
    #
    #>1m9xE|SafeMode
    #MVNPTVFFDIAVDGEPLGRVSFELF
    #Pred_score:    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    #Prediction:    00000000000000000000000000000000
    #Actual:        00110000000001110000000000000000

    my $predictFL = shift @_;
    my ( @seq, @predictionScore, @pInt, $mode );

    open( PREDIC, "<$predictFL" )
      || die("Cannot open $predictFL:$!");
    while (<PREDIC>) {
        s/[\n\r]//mg;

        if (/^(\s*)>[\w\:\s\t\|\.]+\|([\w\s]+)/) {
            $mode = $2;
            next;
        }
        if (/^\s*([a-zA-Z]+)[\s\t]{0,}$/) {

            @seq = split( //, $1 );
            next;
        }
        if (/^Pred_score:\t+([,\?,\.\de]+)/) {
            @predictionScore = split( /,/, $1 );
            next;
        }
        if (/^Prediction:\t+([01\?]+)/) {
            @pInt = split( //, $1 );

        }

    }
    close(PREDIC);

    if ( !@seq ) {
        die("seq unread.  Check $predictFL:$!");
    }

    if ( !@pInt ) {
        die("pInt unread. Check $predictFL:$!");
    }

    if ( !@predictionScore ) {
        die("predictionScore unread. Check $predictFL:$!");

    }

    return ( \@seq, \@predictionScore, \@pInt, $mode );

}

sub writePredictionModeFL {
    my ( $mode, $qryPairs, $qryPredModeFL ) = @_;
    open( OUTPUT, ">>$qryPredModeFL" ) || die("Cannot open $qryPredModeFL:$!");
    foreach my $pair (@$qryPairs) {
        print OUTPUT "$pair $mode\n";
    }
    close OUTPUT;
}

sub batch {

    #    $templateStats->{$qryPairID}->{$templateID}->{$header}
    #    my ($templateStats , $headers) = &readStatisticsFL_all($FullStatFL);
    #

    my $jobDIR        = shift @_;         # .../uploadData/jobID
    my @qryPairs      = @{ shift @_ };    # (A:B, P:Q)
    my $templateStats = shift @_;
    my $Mode          = shift @_;

    my $scoreThr = 0.5;                   #shift @_;    #0.5

    #    my @taskIDs = &readTaskIDs($qryLst);

    my @qryCanbePredicted;
    my @qryCannotbePredicted;

#&getModeParameters : returns $EvalueThr, $logEvalThr, $positiveScoreThr, $logAlignmentLengthThr, $prop_alignLenQryHomThr for KNN2.pl
#&getModeParameters: returns $CC_cutoff for KNN2_weighted.pl

    my $CC_cutoff = &getModeParameters($Mode);

    print LOG "\nmode: $Mode; CC_cutoff: $CC_cutoff\n\n";

    foreach my $qryPair (@qryPairs) {

        print LOG "\n\n ** Int prediction for $qryPair  **\n\n";

        my ( $qry1, $qry2 ) = split( /:/, $qryPair );

#        my $ProteinToBePredictedDIR =
#          "$jobDIR/$qryPair";    #each query is given a specific directory

        my ( $ProteinToBePredictedDIR, $c1, $c2 ) = &getPairDIR2( $qryPair, $jobDIR );

        my $predictionResultsDIR = "$ProteinToBePredictedDIR/predictionResults";

        my $statPredFL   = "$predictionResultsDIR/prediction.stat";
        my $predictionFL = "$predictionResultsDIR/KNN.prediction";

        #-------
        my $templateStats_oneQryPair = $templateStats->{$qryPair};
        if ( !defined $templateStats_oneQryPair ) {
            $templateStats_oneQryPair = $templateStats->{"$qry2:$qry1"};
        }

        if ( !defined $templateStats_oneQryPair ) {
            die("qry pair $qryPair does not have template info in statFL:$!");
        }

        #-------
        my $flag = &extractWeightedKNNpredict(
            $Mode,      $ProteinToBePredictedDIR,
            $CC_cutoff, $templateStats_oneQryPair
        );

        # collect qrys can be predicted and cannot be predicted
        if ( $flag == 0 ) {

            #$taskID can find homologs in $Mode
            push @qryCanbePredicted, $qryPair;
            print LOG "$qryPair can find homologs in $Mode mode.\n";
        }
        elsif ( $flag == 1 ) {
            push @qryCannotbePredicted, $qryPair;
            print LOG "$qryPair can NOT find homologs in $Mode mode.\n";
        }

    }

    return ( \@qryCanbePredicted, \@qryCannotbePredicted );

}

sub extractWeightedKNNpredict {

    #!/usr/bin/perl -w

#Author: Xue, Li
#Date: May, 2008
#This script is part of HomoPPI pipeline
#Extract prediction for each residue of the query sequence, and write into a file, e.g. "1luk.prediction" under the directory "predictionResults"
#
#INPUT  :  localAlign file
#OUTPUT 1. ..data/taskID/predictionResults/taskID.KNN.prediction
#OUTPUT 2. homologs_used_in_prediction.txt
#
#If this script is used as evalution purpose, that is, there exists $taskID.int file that stores the actual interface info, then combine the predicted and actual interface residue together into the output file.

    use strict;
    use File::Basename;

    our $scoreThr;

    #variables
    my $Mode                     = shift @_;
    my $qrypairDIR               = shift @_;    # .../uploadData/jobDIR/A:B
    my $CC_cutoff                = shift @_;
    my $templateStats_oneQryPair = shift @_;

    my $flag =
      0;   #0: this taskID can be predicted. 1: this taskID cannot be predicted.

    #file variables
    my $qryPair = basename($qrypairDIR);    # A:B
    my ( $qry1, $qry2 ) = split( /:/, $qryPair );
    my $localAlignFL1 =
      "$qrypairDIR/$qry1/seq_interface_homologs.localAlign"
      ;                                     #input file of KNN.pl
    my $localAlignFL2 =
      "$qrypairDIR/$qry2/seq_interface_homologs.localAlign"
      ;                                     #input file of KNN.pl

    #--- weighted KNN
    my ( $predictionScores_A, $predictionScores_B,
        $similarityScores_finalcomplex )
      = &KNN2_weighted( $CC_cutoff, $templateStats_oneQryPair, $localAlignFL1,
        $localAlignFL2 )
      ; #- $predictionScores_A = '0.3,0.4,0...', $templates_used_in_prediction = ('1ahjA,2:1ahjB,34', ...)
    my $predictionScores->{$qry1} = join( ',', @$predictionScores_A );
    $predictionScores->{$qry2} = join( ',', @$predictionScores_B );

    #--- write the file of homologs_used_in_prediction.txt
    my $finalHomologFL = "$qrypairDIR/templates_used_in_prediction.lst";
    my @templates_used_in_prediction_sorted = sort {
        $similarityScores_finalcomplex->{$b}
          cmp $similarityScores_finalcomplex->{$a}
    } keys %{$similarityScores_finalcomplex};
    my @similarityScores = map { $similarityScores_finalcomplex->{$_} }
      @templates_used_in_prediction_sorted;

    &write_finalHomologFL( \@templates_used_in_prediction_sorted,
        \@similarityScores, $finalHomologFL );

    my $num_finalHomologs = scalar(@templates_used_in_prediction_sorted);
    print LOG "Num of final homologs used in prediction: $num_finalHomologs \n";
    print LOG
"Final homologs and their similarity scores are saved in $finalHomologFL.\n\n";

    my $prediction_binary->{$qry1} =
      &getBinaryPrediction( $predictionScores_A, $scoreThr );

    $prediction_binary->{$qry2} =
      &getBinaryPrediction( $predictionScores_B, $scoreThr );

    #--- write prediction file

    foreach my $qry ( ( $qry1, $qry2 ) ) {

        # $qry = '1ahjA'
        my $fastaFL = "$qrypairDIR/$qry/$qry.fasta.txt";
        my $intFL_Q = "$qrypairDIR/$qry/$qry.int";         #interface file
        my $seq_qry = &readFastaFL($fastaFL);
        my $actInt_tmp;
        if ($intFL_Q) {
            ( my $seq_Q, $actInt_tmp, my $length_Q, my $numInt_Q ) =
              &readINTFL($intFL_Q);
        }
        my $actInt = join( '', @{$actInt_tmp} );

        #-- prepare result folder
        my $predictionResultsDIR = "$qrypairDIR/$qry/predictionResults";
        my $predictionFL         = "$predictionResultsDIR/$qry.KNN.prediction";
        mkdir $predictionResultsDIR unless ( -d $predictionResultsDIR );

        #-- write prediction files
        &writePredictionFL(
            $Mode, $qry, $seq_qry,
            $predictionScores->{$qry},
            $prediction_binary->{$qry},
            $actInt, $predictionFL
        );

        #--calculate TP, TN, FP, FN, CC
        my $statPredFL = "$predictionResultsDIR/prediction.stat";
        &statisticsOfPrediction( $predictionFL, $scoreThr, $statPredFL );
    }

    #---
    if ( $prediction_binary->{$qry1} =~ /^\?+$/ ) {

        print LOG
"$qryPair cannot find homologs in $Mode zone. No prediction is made.\n";
        $flag = 1;
    }

    print LOG "\n\n\nWEIGHTED KNN PREDICTION DONE.\n\n\n\n";

    return $flag;

}

sub getModeParameters {
    use Switch;
    our $safeMode_thr;    # 0.7
    our $twiMode_thr1;    # 0.5
    our $twiMode_thr2;    # 0.4
    our $twiMode_thr3;    # 0.2
    our $darkMode_thr;    # 0

    my $Mode = shift @_;
    my $CC_cutoff;

    switch ($Mode) {

        #	SafeMode, TwilightMode, DarkMode
        case 'SafeMode' {

            $CC_cutoff = $safeMode_thr;    #0.7
        }
        case 'TwilightMode1' {
            $CC_cutoff = $twiMode_thr1;    #0.5
        }
        case 'TwilightMode2' {
            $CC_cutoff = $twiMode_thr2;    #0.4
        }

        case 'TwilightMode3' {
            $CC_cutoff = $twiMode_thr3;    #0.2
        }
        case 'DarkMode' {
            $CC_cutoff = $darkMode_thr;    #0.0
        }
        else {
            die(
"Please set a valid mode: SafeMode, TwilightMode1, TwilightMode2, TwilightMode3, DarkMode:$!"
            );
        }

    }

    return $CC_cutoff;

}

sub KNN2_weighted {

#Xue, Li
#May, 2008
#
#this script is part of HomoPPI pipeline.
# weighted version of KNN2.pl
#
#use the mean of K-nearest neighbour of the query sequence as prediction
#Only homologs with Similarity score > = $CC_cutoff will be used in predition.
#
#If no homologs has local alignment length larger than alignmentLengthThreshold, no prediction will be made.
#
#read statistics_$basename.txt for PositiveS, Eval,etc to determine which homolog to be used as template.
#
#Usage: perl ./KNN2_weighted.pl similarity_cutOff statisticsFL seq_interface_homologs.localAlign k log_localAlignmentLen_Thr $logEvalThr $PositiveSThr
#$PositiveSThr : the threshold for mean of PositiveS_H and PositiveS_HP. range: [0-100]
#e.g.,
#perl ./KNN2_weighted.pl 0.5 ../data/trans212/statistics_trans212pairedChains_onlyTransInterfaces.txt ../data/trans212/1jiw/1jiw_I_P/1jiwP_P_I/seq_interface_homologs.localAlign

    use File::Basename;
    use List::Util qw(sum min max);

    #---
    my $CC_cutoff = shift @_;
    my $templateStats_oneQryPair =
      shift @_;    #   $templateStats->{$templateID}->{pred_CC} = $pred_CC;
    my $alignFL1 = shift
      @_;    #-- localAlignFL for qry1, used for extracting homolog's interface
    my $alignFL2 = shift @_;    #-- localAlignFL for qry2

    my %similarityScore;
    my %homologs_Seq; #to store sequences of all the homologs in the input files
    my %homologs_int1;    #to store interfaces of  all the homologs (A')
    my %homologs_int2;    #to store interfaces of  all the homologs (B')

    #if qry is peptide
    #if($Length->{$queryID} <= 150){
    #	$logAlignmentLengthThr=0.5*$logAlignmentLengthThr;
    #
    #}

    #get homolog IDs that meet the thresholds
    my @finalHomologIDs;    #the homologs that pass the thresholds
    my $similarityScores_finalcomplex;
    my ( $qryLen1, $qryLen2 );

    foreach my $homologPair ( keys %{$templateStats_oneQryPair} ) {

        #- $homologPair = '4lqwA,1:4lqwC,26'

        #--- get query lengths
        $qryLen1 = $templateStats_oneQryPair->{$homologPair}->{'len_Qry1'};
        $qryLen2 = $templateStats_oneQryPair->{$homologPair}->{'len_Qry2'};

#--- if the aligned part does not contain interfacial residues, exclude this homolog pair
        if ( $templateStats_oneQryPair->{$homologPair}->{'numInt_Homolog1'} <= 3
            && $templateStats_oneQryPair->{$homologPair}->{'numInt_Homolog2'}
            <= 3 )
        {
            print LOG
"$homologPair aligned part has <= 3 interfacial residues. It is not used in prediction.\n";
            next;
        }

        my $predCC = $templateStats_oneQryPair->{$homologPair}->{'pred_CC'};

        if ( !defined $predCC ) {
            die(
"$homologPair does not have predicted CC in the statFL \(statistics_wo_sameProt.txt or statistics.txt\) :$!"
            );
        }

        if ( $predCC >= $CC_cutoff ) {
            push @finalHomologIDs, $homologPair;
            $similarityScores_finalcomplex->{$homologPair} = $predCC;
        }

    }

    my ( @predictionScores_A, @predictionScores_B );
    if ( !@finalHomologIDs ) {
        @predictionScores_A = ('?') x $qryLen1;
        @predictionScores_B = ('?') x $qryLen2;

        return ( \@predictionScores_A, \@predictionScores_B,
            $similarityScores_finalcomplex );
    }

    #-- read alignFL to get %homologs_int
    my %similarityScores;
    foreach my $homologPair (@finalHomologIDs) {

        # $homologPair = '1ahjA,3:1ahjB,20'
        my ( $homolog1_tmp, $homolog2_tmp ) = split( /:/, $homologPair );
        my $Aprime_chn = substr( $homolog1_tmp, 4, 1 );
        my $Bprime_chn = substr( $homolog2_tmp, 4, 1 );
        my ( $homolog1, $groupID1 ) = split( /,/, $homolog1_tmp );
        my ( $homolog2, $groupID2 ) = split( /,/, $homolog2_tmp );

        my $Aprime =
          "$homolog1|$Aprime_chn:$Bprime_chn,$groupID1";    # 1ahjA|A:B,3
        my $Bprime =
          "$homolog2|$Bprime_chn:$Aprime_chn,$groupID2";    # 1ahjB|B:A,20

        $homologs_int1{$Aprime} = &getHomologInt( $Aprime, $alignFL1 );
        $homologs_int2{$Bprime} = &getHomologInt( $Bprime, $alignFL2 );

        $similarityScores{$Aprime} =
          $templateStats_oneQryPair->{$homologPair}->{'pred_CC'};
        $similarityScores{$Bprime} =
          $templateStats_oneQryPair->{$homologPair}->{'pred_CC'};

    }

    #--- calculate weighted avarge of homolog interfaces
    @predictionScores_A =
      @{ &weightedAverageHomologInt( \%homologs_int1, \%similarityScores ) };
    @predictionScores_B =
      @{ &weightedAverageHomologInt( \%homologs_int2, \%similarityScores ) };

    return ( \@predictionScores_A, \@predictionScores_B,
        $similarityScores_finalcomplex );
}

sub weightedAverageHomologInt {

# Given multiple homologs' interface for one query, return the average score for each column.
# Note: all homologs' interface have equal length of the query.
    use List::Util qw (sum);

    my %homologInts = %{ shift @_ };

    # $homologInts{1ahjA|A:B,2} = '0101000000111000000'

    my %similarityScores = %{ shift @_ };

    # $similarityScores{1ahjB|B:A,34} = 0.8

    #-- predict using KNN

#- 1. store interfaces as a matrix, and the last column of the matrix is the similarity score
    my @homologs = keys %homologInts;
    my @matrix
      ;    #each row is an interface, the last column is the similarity score

    my $num_row = scalar( keys %homologInts );
    my $num_col;
    for my $i ( 0 .. scalar @homologs - 1 ) {

        my @tmp = split( //, $homologInts{ $homologs[$i] } );
        my $similarScore = $similarityScores{ $homologs[$i] };
        push @tmp, $similarScore;
        $num_col = scalar @tmp;
        push @matrix, [@tmp];
    }

    #- 2. multiply each element in @matrix by weight (Similarity score)
    my $qry_length = $num_col - 1;  #-- the last column of matrix is the pred_CC
    my @predictionScores = (0) x $qry_length;
    for my $y ( 0 .. $qry_length - 1 ) {
        for my $x ( 0 .. $num_row - 1 ) {
            my $weight = $matrix[$x][-1];

            if ( $matrix[$x][$y] eq '?' || $matrix[$x][$y] eq '-' ) {
                $predictionScores[$y] = $predictionScores[$y] + 0;
            }
            else {
                $predictionScores[$y] =
                  $predictionScores[$y] + $matrix[$x][$y] * $weight;
            }
        }
    }

    #- 3. normalize @weighted_matrix by the sum of weights
    my $sum_weights = sum( values %similarityScores );
    @predictionScores = map {
        $_ = $_ / $sum_weights;
        if ( $_ ne '0' ) { sprintf( "%.2f", $_ ) }
        else             { $_; }
    } @predictionScores;

    #--
    return \@predictionScores;

}

sub write_finalHomologFL {

    my @finalHomologIDs  = @{ shift @_ };
    my @similarityScores = @{ shift @_ };
    my $finalHomologFL   = shift @_;

    my $num_finalHomologs = scalar @finalHomologIDs;

    unlink $finalHomologFL if ( -e $finalHomologFL );
    open( OUT, ">> $finalHomologFL" )
      or die("Cannot open $finalHomologFL: $!");
    print OUT "#homologous complexes used in final prediction.\n";
    print OUT
      "#Num of final homologs used in prediction: $num_finalHomologs \n";
    print OUT "#\n#Template_ID\tpredicted_CC\n";
    for ( my $i = 0 ; $i < scalar @finalHomologIDs ; $i++ ) {
        my $similarityScore = sprintf( "%.3f", $similarityScores[$i] );
        print OUT "$finalHomologIDs[$i]\t$similarityScore\n";
    }
    close(OUT);

    print LOG "$finalHomologFL generated.\n";
}

sub getBinaryPrediction {
    my @scores      = @{ shift @_ };    #(0.25,0.25,0.25,0)
    my $scoreCutoff = shift @_;
    my $prediction  = '';

    foreach my $score (@scores) {

        if ( $score eq '?' ) {

            $prediction = $prediction . '?';
        }
        elsif ( $score >= $scoreCutoff ) {
            $prediction = $prediction . '1';
        }
        else {
            $prediction = $prediction . '0';
        }
    }

    #double check
    if ( length($prediction) ne scalar @scores ) {
        die(
            "Binary prediction and prediction scores have different lengths:$!"
        );
    }

    return $prediction;

}

sub readINTFL {

    #read in $taskID's interface and sequences

    my $interfaceFL = shift @_;
    my $seq;
    my @int;    # array ref
    my $num_int = 0;    #number of interface residues
    my $length;         #the length of the sequence

    open( INT_Q_FL, "<$interfaceFL" ) || die("Cannot open $interfaceFL:$!");
    while (<INT_Q_FL>) {
        s/[\n\r]//mg;
        if (/^([A-Za-z]+)/) {
            $seq = $1;
            next;
        }
        if (/^([01\?\-]+)/) {
            @int = split //, $1;
            last;
        }
    }

    close INT_Q_FL;

    $length = length($seq);

    foreach (@int) {
        if (/1/) {
            $num_int++;
        }
    }

    my $int_string = join( '', @int );

    if ( $int_string =~ /^\?+$/ ) {

        #int is all question marks
        $num_int = 'Nan';
    }

    if ( scalar @int ne $length ) {

        die(
"int and seq have different length.\nint:@int\nseq:$seq\nCheck $interfaceFL:$!"
        );
    }

    return ( $seq, \@int, $length, $num_int );

}

sub statisticsOfPrediction {

#Author: Xue, Li
#Date: Jan 10, 2012
#
#This script is part of HomDDI
#To get the prediction statistics, such as TP, TN, FP, FN, CC (coorelation coefficient)
#
#Question marks in predictions are replaced with zeros.
#
#Input file format:
#
#				    >d2ez8b3
#	                KQEGPLQAYQVLRAVNKIAEPDAIYSIDVGDINLNANRHLKLTPSNRHITSNLFATMGVGIPGAIAAKLNYPERQVFNLAGDGGASMTMQDLATQVQYHLPVINVVFTNCQYGWIKDEQEDTNQNDFIGVEFNDIDFSKIADGVHMQAFRVNKIEQLPDVFEQAKAIAQHEPVLIDAVITGDRPLPAEKLRLDSAMSSAADIEAFKQRYEAQDLQPLSTYLKQFGLDD
#	Pred_score:     ?,?,?,?,?,?,?,?,?,?,0,0,0,0,0,0,0,0,0,0,0,0,0,0.35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0.5,0.8,1,1,0.5,1,0.5,0,0.5,0,0,0.5,0.5,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0.5,0.5,0,0,0.5,0.5,0,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0.5,0,0,0.5,0.3,0,0,0,0.5,0,0,0.5,0.45,0.5,0.5,0.5,0,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,?,?,0,0,0,0,0,0,0,0,0,0,0,0,0,0,?,?,0,0,0,0,0,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
#	Prediction:     ??????????000000000000000000000000000000000000001011111110100110111011100000000000100110011011000000000000000001001001000010010111011000000000000000000000000000000??00000000000000??00000??????????????????????????????????????????
#	Actual:         000000000000000000000001000000000000000000000010111111011110011001101110000000000000000101101100110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

#
#Output: TP, TN, FP, FN, CC
#
#Usage: ./statPrediction.pl ../data/$taskID/predictionResults/$taskID.prediction scoreThr
#perl ./statPrediction.pl  ../data/test/1a2kA_1a2kC/1a2kC/predictionResults/1a2kC.KNN.prediction 0.5

    use strict;

    if ( scalar @_ ne 3 ) {
        die "statPrediction.pl: Please input prediction file and scoreThr:$!";
    }

    my $inputFL  = shift @_;
    my $scoreThr = shift @_;
    my $outputFL = shift @_;
    my $TN       = 0;
    my $FP       = 0;
    my $FN       = 0;
    my $TP       = 0;
    my $F1;
    my $CC;
    my $specificity;
    my $sensitivity;
    my $accuracy;
    my @predictionScore;
    my @prediction;
    my @actualInterface;

    #my $flag=0;#flag "=== Error on test data ==="

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");

    open( INPUT, "<$inputFL" ) || die("Cannot open $inputFL:$!");
    foreach (<INPUT>) {

        if (/^Pred_score:[\s\t]+([01\d\.\?,\-e\w!]+)/) {

            @predictionScore = split( /[,\s]+/, $1 );
            next;
        }
        if (/Actual:[\s\t]+([01\?\-]+)/) {
            @actualInterface = split( //, $1 );
            last;
        }
    }
    close(INPUT);

    if ( !@predictionScore || !@actualInterface ) {
        die("Check $inputFL:$!");
    }

    #change prediction score to 0 or 1

    foreach (@predictionScore) {
        if ( $_ eq '?' ) {

#replace QNs with zeros, i.e., regarding these residues are predicted as non-interface sites.
            push @prediction, 0;

            #		push @prediction, '?';
            next;
        }
        if ( $_ >= $scoreThr ) {
            push @prediction, 1;
        }
        else {
            push @prediction, 0;
        }
    }

    if ( scalar @prediction ne scalar @actualInterface ) {
        die("the number of prediction is different from that of actual int:$!");
    }

    #count TP,TN, FP, FN

    if (@prediction) {

        my $i = 0;
        foreach (@actualInterface) {
            if ( $_ ne '?' && $prediction[$i] ne '?' ) {

                if (/1/) {

                    if ( $prediction[$i] == 1 ) {
                        $TP++;
                    }
                    else {
                        $FN++;
                    }
                }
                else {
                    if ( $prediction[$i] == 0 ) {

#xue,check here:					Quantifier follows nothing in regex; marked by <-- HERE in m/? <-- HERE / at statPrediction.pl line 48.

                        $TN++;
                    }
                    else {
                        $FP++;
                    }
                }
            }
            $i++;
        }

        if (
            ( $TP + $FN ) * ( $TP + $FP ) * ( $TN + $FP ) * ( $TN + $FN ) != 0 )
        {
            $CC =
              ( $TP * $TN - $FP * $FN ) /
              sqrt(
                ( $TP + $FN ) * ( $TP + $FP ) * ( $TN + $FP ) * ( $TN + $FN ) );
        }

        elsif (
            ( $TP + $FN ) * ( $TP + $FP ) * ( $TN + $FP ) * ( $TN + $FN ) == 0 )
        {

#when CC is not defined, use the approximation of CC (Burset and Guigo 1996, Genomics; Baldi 2000)
#
            if ( $TP == 0 && $TN == 0 && $FP == 0 && $FN == 0 ) {

                #no prediction is made. all '?'.
                $CC = '';
            }
            else {
                print OUTPUT
"#The approximation of CC is used (Burset and Guigo 1996, Genomics; Baldi 2000)\n";
                $CC = &AC( $TP, $TN, $FP, $FN );
            }

        }

  #	elsif ( $TP + $FP == 0 ) {
  #
  #		#When all the residues are predicted as non-interface
  #		$CC = 0;
  #	}
  #	elsif ( $TN + $FN == 0 && $FP != 0 ) {
  #
  #		#When all the residues are predicted as interfaces
  #		# and there are some acutral non-interfaces
  #		$CC = 0;
  #	}
  #	elsif ( $TN == 0 && $FP == 0 && $FN == 0 ) {
  #
  #	  #When all the residues are actually interfaces and the prediction is right
  #		$CC = 0;
  #	}
  #	else {
  #		$CC = '';
  #	}

        if ( $TP + $FP != 0 ) {
            $specificity = $TP / ( $TP + $FP );
        }
        elsif ( $TP + $FP == 0 ) {
            $specificity = 0;
        }
        else {
            $specificity = 1;
        }

        if ( $TP + $FN != 0 ) {

            $sensitivity = $TP / ( $TP + $FN );
        }
        else {
            $sensitivity = 1;
        }

        if ( $TP + $TN + $FN + $FP == 0 ) {
            $TN = 'Nan';
            $TP = 'Nan';
            $FN = 'Nan';
            $FP = 'Nan';

            $CC          = 'Nan';
            $specificity = 'Nan';
            $sensitivity = 'Nan';
            $accuracy    = 'Nan';

        }
        else {

            $accuracy = ( $TP + $TN ) / ( $TP + $TN + $FN + $FP );
        }

        if ( $sensitivity eq 'Nan' || $specificity eq 'Nan' ) {
            $F1 = 'Nan';
        }
        elsif ( $sensitivity != 0 && $specificity != 0 ) {

            $F1 =
              2 * $specificity * $sensitivity / ( $specificity + $sensitivity );
        }
        else {
            $F1 = 'Nan';
        }

    }

    else {
        $TN = 'Nan';
        $TP = 'Nan';
        $FN = 'Nan';
        $FP = 'Nan';

        $F1          = 'Nan';
        $CC          = 'Nan';
        $specificity = 'Nan';
        $sensitivity = 'Nan';
        $accuracy    = 'Nan';

        print OUTPUT "No prediction is made.\n";
    }

    #-----------------
    #
    my $predictionScores = join( '', @predictionScore );
    if ( $predictionScores =~ /^\?+$/ ) {
        $TN = 'Nan';
        $TP = 'Nan';
        $FN = 'Nan';
        $FP = 'Nan';

        $F1          = 'Nan';
        $CC          = 'Nan';
        $specificity = 'Nan';
        $sensitivity = 'Nan';
        $accuracy    = 'Nan';

        print OUTPUT "No prediction is made.\n";

    }

    #-----------------
    #
    ( $F1, $CC, $specificity, $sensitivity, $accuracy ) = map {

        if (/Nan/i) {
            $_ = $_;
        }
        $_ = sprintf( "%.2f", $_ );

    } ( $F1, $CC, $specificity, $sensitivity, $accuracy );

    print OUTPUT
"#Predictions marked with question marks are replaced with zeros, i.e., these residues are considered predicted non-interface.\n";
    print OUTPUT "TN=$TN; TP=$TP; FN=$FN; FP=$FP\n";
    print OUTPUT
"F1= $F1; CC= $CC; specificity= $specificity; sensitivity = $sensitivity; accuracy= $accuracy\n";

}

#--------------------------------
sub getHomologInt {

#   get the int for $homolog from $alignFL
#
# $alignFL:
# >1m9xA|A:G,GroupID 50,Alignment_Length 164,Sbjct_Length 165,Sbjct_RESI 2 - 165,Query_RESI 2 - 165,Score = 280.0 bits, Expect = 1.1E-45, Identities = 100.0%, Positives = 100.0%, Similarity = 4.7, hhpred = 100.0
# VNPTVFFDIAVDGEPLGRVSFELFADKVPKTAENFRALSTGEKGFGYKGSCFHRIIPGFMCQGGDFTRHNGTGGKSIYGEKFEDENFILKHTGPGILSMANAGPNTNGSQFFICTAKTEWLDGKHVVFGKVKEGMNIVEAMERFGSRNGKTSKKITIADCGQLE
# 10000000000000000000000011100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
#

    use strict;

    my $homolog      = shift @_;    #1im9A|A:D,3
    my $localAlignFL = shift @_;

    #    print LOG "Read localAlignFL for homolog's int (homolog: $homolog)\n";

    #---
    #read local align file
    my ( $sbjIDs, $queryLength, $localAlignResults ) =
      &readLocalAlignFL($localAlignFL)
      ;    #$alignmentLength->{"$homolog,$groupID"} =23; $homolog = '1o96Q|Q:Z'

    my $seq_alignedPart = $localAlignResults->{'QrySeqAligned'}->{$homolog};
    my $int_alignedPart = $localAlignResults->{'HomologIntAligned'}->{$homolog};

    if ( !defined $seq_alignedPart ) {
        die(
"Aligned part of seq of homolog $homolog is not read from $localAlignFL:$!"
        );
    }

    if ( !defined $int_alignedPart ) {
        die(
"Aligned part of int of homolog $homolog is not read from $localAlignFL:$!"
        );
    }

    #find the index for - in the $seq_alignedPart
    my @where = ();
    my $index = index( $seq_alignedPart, '-' );

    while ( $index >= 0 ) {

        if ( $index >= 0 ) {
            push @where, $index;
        }
        $index = index( $seq_alignedPart, '-', $index + 1 );
    }

    #remove the interface sign that corresponds to - in $seq_alignedPart

    my $start_Q = $localAlignResults->{'start_Q'}->{$homolog};
    my $end_Q   = $localAlignResults->{'end_Q'}->{$homolog};

    my @tempINT = split( //, $int_alignedPart );

    if (@where) {
        foreach (<@where>) {
            $tempINT[$_] = '';
        }

    }
    $int_alignedPart = join( '', @tempINT );

    #add ? to the two ends
    my $multiX1 = "?" x ( $start_Q - 1 );            #$num_multiX question marks
    my $multiX2 = "?" x ( $queryLength - $end_Q );
    my $homologInt = $multiX1 . $int_alignedPart . $multiX2;

    #-----

    if ( !defined $homologInt ) {

        die("Int of homolog $homolog is not read from $localAlignFL: $!");
    }

    return $homologInt;
}

#---------------

sub AC {

    #the approximation of CC (Burset and Guigo 1996, Genomics; Baldi 2000)
    my ( $TP, $TN, $FP, $FN ) = @_;

    #Average Conditional Probability (ACP)
    my $ACP = &ACP( $TP, $TN, $FP, $FN );
    my $AC = ( $ACP - 0.5 ) * 2;
    return $AC;

}

sub ACP {

#Average Conditional Probability (ACP) (Burset and Guigo 1996, Genomics; Baldi 2000)

    use strict;
    use List::Util qw(sum);
    my ( $TP, $TN, $FP, $FN ) = @_;
    my $ACP;

    my @p;
    if ( $TP + $FN != 0 ) {
        push @p, $TP / ( $TP + $FN );
    }
    if ( $TP + $FP != 0 ) {
        push @p, $TP / ( $TP + $FP );
    }
    if ( $TN + $FP != 0 ) {
        push @p, $TN / ( $TN + $FP );
    }
    if ( $TN + $FN != 0 ) {
        push @p, $TN / ( $TN + $FN );
    }

    $ACP = sum(@p) / scalar(@p);

}
sub rmChainsWithLessThan5intAA_PS {


    #remove chains with interface residues less than 5 aa.
    #Li Xue
    #Apr 12th, 2010
    #
    #remove chains with <= 5 interfacial AAs
    #If a chain is removed, it interacting parner is also removed.
    #
    #OUTPUT 1: the chains with less than $num_intAA_thr int are removed from the input file
    #OUTPUT 2: seq_int_homologWithLessThan5int.lst (contains seq and int info of the removed chains)
    #
    #perl ./rmChainsWithLessThan5intAA.pl 3 ../data/seq_int_qryComplexes_test.lst

    use strict;
    use File::Copy;
    use File::Basename;

    my $num_intAA_thresh =
      shift @_
      ; #seq with less than 5 interface residues are removed from the input file.
    my $inputFL = shift @_; # seq_int_templates.lst

    #-- output 1
    my $outputFL1 = "$inputFL";

    #-- output 2
    my $dirname = dirname($inputFL);
    my $seq_int_lessThan5intFL =
      "$dirname/seq_int_templatesWithLessThan5int.lst"
      ;    #(contains seq and int info of the removed chains)


    #--- read input file: seq_int_templates.lst
    my ( $comment, $seqs_complexes, $int_complexes ) =
      &readSeqIntFL($inputFL);

    #--- remove homologs with <= 5 int
    my @keep;
    my @remove;
    foreach my $homologID (keys %{$int_complexes}){
        # $homologID = '3heiI|I:J'
            my $int     = $int_complexes->{$homologID};
            my $numOnes = &numOnes($int);

            if ( $numOnes < $num_intAA_thresh ) {
                my $partner = &partnerID($homologID); # '3heiJ|J:I'
                push @remove, $homologID ;
                push @remove, $partner ;
            }
    }

    #-- if a homolog exists in @remove, remove it from @keep
    @keep = keys %{$int_complexes};
    if (@remove){
        @keep = &furtherRm(\@keep, \@remove);
    }

    #---
    &write_final_seq_intFL(\@keep, $seqs_complexes, $int_complexes, $comment, $outputFL1);
    &write_final_seq_intFL(\@remove, $seqs_complexes, $int_complexes, $comment , $seq_int_lessThan5intFL);
}

sub write_final_seq_intFL{
    my @homologIDs = @{shift @_};
    my %seqs = %{shift @_};
    my %ints = %{shift @_};
    my $comment = shift @_;
    my $outputFL = shift @_;

    unlink $outputFL if (-e $outputFL);

    open (OUTPUT, ">>$outputFL") or die ("Cannot open $outputFL:$!");
    print OUTPUT "$comment\n";
    foreach my $homologID (@homologIDs){
        my $seq = $seqs{$homologID};
        my $int = $ints{$homologID};
        print OUTPUT ">$homologID\n";
        print OUTPUT "$seq\n";
        print OUTPUT "$int\n";
    }
    close OUTPUT;

    print LOG "$outputFL generated\n";
}

sub furtherRm{
    #-- remove any elements in @remove from @keep
    my @keep = @{shift @_};
    my @remove = @{shift @_};

    my %seen;
    @seen{@remove}= (1) x scalar @remove;

    my @final = grep {! defined $seen{$_}} @keep;

    if (! @final){
        die("ERROR: No homolog is kept.\nRemove: @remove\nKeep_ori: @keep:$!");
    }

    return @final;
}


#sub rmChainsWithLessThan5intAA_PS_old {
#
##remove chains with interface residues less than 5 aa.
##Li Xue
##Apr 12th, 2010
##
##remove chains with <= 5 interfacial AAs
##If a chain is removed, it interacting parner is also removed.
##
##OUTPUT 1: the chains with less than $num_intAA_thr int are removed from the input file
##OUTPUT 2: seq_int_homologWithLessThan5int.lst (contains seq and int info of the removed chains)
##
##perl ./rmChainsWithLessThan5intAA.pl 3 ../data/seq_int_qryComplexes_test.lst
#
#    use strict;
#    use File::Copy;
#    use File::Basename;
#
#    my $num_intAA_thresh =
#      shift @_
#      ; #seq with less than 5 interface residues are removed from the input file.
#    my $inputFL = shift @_;    # seq_int_templates.lst
#
#    #-- output 1
#    my $outputFL1 = $inputFL;
#
#    #-- output 2
#    my $dirname = dirname($inputFL);
#    my $seq_int_lessThan5intFL =
#      "$dirname/seq_int_templatesWithLessThan5int.lst"
#      ;    #(contains seq and int info of the removed chains)
#
#    #--- read input file: seq_int_templates.lst
#    my ( $comment, $seqs_complexes, $int_complexes ) = &readSeqIntFL($inputFL);
#
#    #--- remove homologs with <= 5 int
#    my @keep;
#    my @remove;
#    foreach my $homologID ( keys %{$int_complexes} ) {
#
#        # $homologID = '3heiI|I:J'
#        my $int     = $int_complexes->{$homologID};
#        my $numOnes = &numOnes($int);
#
#        if ( $numOnes < $num_intAA_thresh ) {
#            my $partner = &partnerID($homologID);    # '3heiJ|J:I'
#            push @remove, $homologID;
#            push @remove, $partner;
#        }
#    }
#
#    #-- if a homolog exists in @remove, remove it from @keep
#    @keep = keys %{$int_complexes};
#    if (@remove) {
#        @keep = &furtherRm( \@keep, \@remove );
#    }
#
#    #---
#    &write_final_seq_intFL( \@keep, $seqs_complexes, $int_complexes,
#        $outputFL1 );
#    &write_final_seq_intFL( \@remove, $seqs_complexes, $int_complexes,
#        $seq_int_lessThan5intFL );
#}
#
sub write_final_seq_intFL_old {
    my @homologIDs = @{ shift @_ };
    my %seqs       = %{ shift @_ };
    my %ints       = %{ shift @_ };
    my $outputFL   = shift @_;

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    foreach my $homologID (@homologIDs) {
        my $seq = $seqs{$homologID};
        my $int = $ints{$homologID};
        print OUTPUT ">$homologID\n";
        print OUTPUT "$seq\n";
        print OUTPUT "$int\n";
    }
    close OUTPUT;

    print LOG "$outputFL generated\n";
}

sub furtherRm_old {

    #-- remove any elements in @remove from @keep
    my @keep   = @{ shift @_ };
    my @remove = @{ shift @_ };

    my %seen;
    @seen{@remove} = (1) x scalar @remove;

    my @final = grep { defined $seen{$_} } @keep;

    if ( !@final ) {
        die("ERROR: No homolog is kept.\nRemove: @remove\nKeep_ori: @keep:$!");
    }

    return @final;
}

sub partnerID {
    my $homologID = shift @_;

    # $homologID = '3heiI|I:J'

    if ( $homologID !~ /^\w{5}\|\w{1}:\w{1}$/ ) {
        die("homolog ID ($homologID) format wrong. Correct format: 3heiI|I:J:$!"
        );
    }

    my ( $prot, $chn1, $chn2 ) = split( /[\|:]/, $homologID );
    my $pdbID = substr( $prot, 0, 4 );
    my $partnerID = "$pdbID$chn2|$chn2:$chn1";

    return $partnerID;
}

sub rmChainsWithLessThan5intAA_PS_old {

#remove chains with interface residues less than 5 aa.
#Li Xue
#Apr 12th, 2010
#
#remove chains with <= interfacing AAs
#If a chain is removed, it interacting parner is also removed.
#
#output1: the chains with less than $num_intAA_thr int are removed from the input file
#output2: seq_int_homologWithLessThan5int.lst (contains seq and int info of the removed chains)
#perl ./rmChainsWithLessThan5intAA.pl 3 ../data/seq_int_qryComplexes_test.lst

    use strict;
    use File::Copy;
    use File::Basename;

    my $num_intAA_thresh =
      shift @_
      ; #seq with less than 5 interface residues are removed from the input file.
    my $inputFL = shift @_;
    my $dirname = dirname($inputFL);

    #--- tmp file
    my $rand  = rand(2);
    my $tmpFL = "$dirname/$rand.tmp";

    #--
    my $seq_int_lessThan5intFL =
      "$dirname/seq_int_templatesWithLessThan5int.lst"
      ;    #(contains seq and int info of the removed chains)
    my $header;
    my $seq;
    my $int;
    my $numOnes;
    my $flag           = 0;    #flag=1: sequences and interfaces starts
    my $num_totalinput = 0;    #number of seqs in the input file
    my $num_final      = 0;    #number of sequences in the final output file
    my @partnersToberemoved
      ; #these proteins have > 5 int, but their partner have <= 5 int. So they are removed, too.

    print LOG
"\n\n\nRemoving any chains with interface residues less than $num_intAA_thresh amino acids...\n\n";

    unlink $seq_int_lessThan5intFL if ( -e $seq_int_lessThan5intFL );

    open( OUTPUT, ">>$seq_int_lessThan5intFL" );
    print OUTPUT "# Homologs with interface <= $num_intAA_thresh amino acids\n";
    close OUTPUT;

    open( INPUT, "<$inputFL" ) || die("Cannot open $inputFL:$!");
    unlink($tmpFL) if ( -e $tmpFL );

    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^>\w+/) {

            #		>1efvA|A:B
            $num_totalinput++;
            $header = $_;
            $flag   = 1;
            next;
        }
        if (/^[A-Za-z]+$/) {
            $seq = $_;
            next;
        }
        if (/^[01\?]+$/) {
            $int     = $_;
            $numOnes = &numOnes($int);

            if ( $numOnes >= $num_intAA_thresh ) {

                $num_final++;

                open( TMP, ">>$tmpFL" );
                print TMP "$header\n";
                print TMP "$seq\n";
                print TMP "$int\n";
                close(TMP);
            }
            else {
                open( OUTPUT, ">>$seq_int_lessThan5intFL" );
                print OUTPUT "$header\n";
                print OUTPUT "$seq\n";
                print OUTPUT "$int\n";
                close(OUTPUT);
                print LOG
"$header has less than $num_intAA_thresh interface residues and is removed from $inputFL.\n";

                #-
                if ( $header =~ />(\w{4})(\w{1})\|\w{0,}:(\w{0,})/ ) {
                    my $pdbID  = $1;
                    my $chain1 = $2;
                    my $chain2 = $3;
                    push @partnersToberemoved, "$pdbID$chain2|$chain2:$chain1"
                      ;    #The partner will also be removed from $inputFL
                }

            }

            $header = '';
            $seq    = '';
            $int    = '';
            next;
        }

        #print header lines of the input file
        if ( $flag == 0 ) {
            open( TMP, ">>$tmpFL" );
            print TMP "$_\n";
            close(TMP);
            next;
        }
    }

    close(INPUT);

    #remove the partners
    #remove @partnersToberemoved from $tmpFL to $seq_int_lessThan5intFL
    my $num_partnersRemoved = 0;
    if (@partnersToberemoved) {

        foreach (@partnersToberemoved) {
            print LOG
">The binding partner of $_  has less than $num_intAA_thresh interface residues and $_ is also removed from $inputFL.\n";

            my $num = &rmProteinFromFL( $tmpFL, $seq_int_lessThan5intFL, $_ );
            $num_partnersRemoved = $num_partnersRemoved + $num;
        }

    }
    $num_final = $num_final - $num_partnersRemoved;

    unlink($inputFL);
    move( $tmpFL, $inputFL ) or die "Move $tmpFL failed: $!";

    print LOG
"Chains with less than $num_intAA_thresh interface residues are removed from $inputFL to $seq_int_lessThan5intFL.\n";
    print LOG
"There are total $num_totalinput input sequences. And now $num_final seqs in $inputFL.\n\n";
}

sub numOnes {
    my $int     = shift @_;
    my $numOnes = 0;
    my @temp    = split //, $int;

    foreach (@temp) {
        if (/1/) {
            $numOnes = $numOnes + 1;
        }
    }
    return $numOnes;

}

sub rmProteinFromFL {

    #remove $proteinIDtoBeRemove from $file1
    #write $proteinIDtoBeRemove to $file2

    my $file1              = shift @_;
    my $file2              = shift @_;
    my $proteinToBeRemoved = shift @_;    #4monB|B:A

    my $dirname = dirname($file1);
    my $rand    = rand(1);
    my $tmpFL   = "$dirname/$rand.tmp";
    my $num     = 0;
    my ( $protein, $chain1, $chain2 ) =
      split( /[\|:]/, $proteinToBeRemoved );

    my $header;
    my $seq;
    my $int;
    my $flag = 0;

    open( INPUT,  "<$file1" )  || die("Cannot open $file1.\n");
    open( OUTPUT, ">>$file2" ) || die("Cannot open $file2:$!");
    unlink $tmpFL if ( -e $tmpFL );

    foreach (<INPUT>) {
        s/[\n\r]//mg;

        if (/^>$protein\|$chain1:$chain2/) {

            #		>1efvA|A:B
            $num    = 1;
            $header = $_;
            print OUTPUT "$header\n";
            $flag = 1;
            print LOG " $_ is removed from $file1.\n";
            next;
        }
        if ( $flag == 1 && /^[A-Za-z]+$/ ) {
            $seq = $_;
            print OUTPUT "$seq\n";
            next;
        }
        if ( $flag == 1 && /^[01\?]+$/ ) {
            $int = $_;
            print OUTPUT "$int\n";
            $flag = 0;
            next;
        }

        #print the rest of lines of the input file to $tmpFL
        if ( $flag == 0 ) {
            open( TMP, ">>$tmpFL" ) || die("Cannot open $tmpFL:$!");
            print TMP "$_\n";
            close(TMP);
            next;
        }
    }

    close(INPUT);

    move( $tmpFL, $file1 ) || die("Cannot move $tmpFL:$!");
    return $num;

    #	print "$proteinToBeRemoved is removed from $file1 and put in $file2.\n";

}

sub existIntInfo {

#    >Q|Q:P
#    VNPTVFFDIAVDGEPLGRVSFELFADKVPKTAENFRALSTGEKGFGYKGSCFHRIIPGFMCQGGDFTRHNGTGGKSIYGEKFEDENFILKHTGPGILSMANAGPNTNGSQFFICTAKTEWLDGKHVVFGKVKEGMNIVEAMERFGSRNGKTSKKITIADCGQLE
#    ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#    >P|P:Q
#    VHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSP
#    ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

    my $intFL = shift @_;
    my $ans =
      0;    # 1: interface information exists. 0: no interface info (all ?s).

    open( INTFL, "<$intFL" ) || die("Cannot open $intFL:$!");
    while (<INTFL>) {
        s/[\n\r]//mg;

        if (/^[01?\-]+$/) {

            if (/^[\?]+$/) {

                #-- the interface line is all ?s
                next;
            }
            else {
                #-- some queries have interface info
                $ans = 1;
                return $ans;
            }

        }
    }
    return $ans;

}

sub readFastaFL {
    my $queryFastaFile = shift @_;    #-- should only contain one sequence
    my $seq            = '';

    open( INPUT, "<$queryFastaFile" )
      || die("Cannot open $queryFastaFile:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^[A-Za-z]+$/) {
            $seq = $seq . $_;
        }
    }
    close INPUT;

    if ( !defined $seq ) {
        die("No seq read from $queryFastaFile:$!");
    }
    return $seq;

}

sub writePredictionFL {

    #write qryID.KNN.prediction file
    use strict;

    my ( $Mode, $qryID, $seq, $predScore, $prediction, $actualInt, $outputFL )
      = @_;

    unlink $outputFL if ( -e $outputFL );

    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    print OUTPUT "\t\t>$qryID|$Mode\n";
    print OUTPUT "\t\t$seq\n";
    print OUTPUT "Pred_score:\t$predScore\n";
    print OUTPUT "Prediction:\t$prediction\n";
    if ( $actualInt ne '' ) {
        print OUTPUT "Actual:\t\t$actualInt\n";
    }
    close OUTPUT;

    print LOG "$outputFL generated.\n";
}

sub getTemplate1_ID {

#-- extract template 1's ID from templates_used_in_prediction.lst
#
# INPUT (.../TemplatesUsed/templates_used_in_prediction.lst):
# #qryPairID homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2 EVal1  EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL1    LAL2    frac_LAL1 frac_LAL2   TP_1    TN_1    FP_1    FN_1    CC_1    specificity_1sensitivity_1  TP_2    TN_2    FP_2    FN_2    CC_2 specificity_2   sensitivity_2
# 1m9xE:1m9xH 4lqwA,1:4lqwC,26    165 146 177 146 14  9   9   11  1.5E-47 2.3E-72 66.297.9    88.9    98.6    3.4 4.3 100.0   100.0  163 146 0.91    1.00    2   142 7   12  0.120.22    0.14    0   115 11  9   -0.08   0.00    0.00
#
# OUTPUT:
#
#            $template1_ID->{$qry1}= $template1_rec;
#            $template1_ID->{$qry2}= $template1_lig;
#
    my $homcomplexes_used_in_predictionFL = shift @_;
    my $template1_ID;
    open( INPUT, "<$homcomplexes_used_in_predictionFL" )
      or die("Cannot open $homcomplexes_used_in_predictionFL:$!");
    while (<INPUT>) {
        s/[\n\r]//gm;
        if (/^\w+:\w+\s+/) {
            my @tmp = split( /\s+/, $_ );
            my ( $qry1,          $qry2 )          = split( /:/, $tmp[0] );
            my ( $template1_rec, $template1_lig ) = split( /:/, $tmp[1] );
            $template1_ID->{$qry1} = $template1_rec;
            $template1_ID->{$qry2} = $template1_lig;

            last;
        }
    }
    close(INPUT);

    if ( !defined $template1_ID ) {
        die(
"template 1's ID not read from $homcomplexes_used_in_predictionFL:$!"
        );
    }

    return $template1_ID;

}

sub genAtomResNumMapFL_new {

#--output:
# A A'
# 2 91
# 3 92
# 4 -
# 5 93
# - 94
# 6 95
#
# --
# Note: we do not know the atomResNum for the query, so atomResNum (the 1st column of the output file) is seqResNum starting from 1.

    use strict;

    our $PPIDBseqDIR;

    my ( $homolog, $blastAlignFL, $outputFL ) = @_;

    #-- $homolog = '1m9xA,1'

    my $pdbID = substr( $homolog, 0, 4 );
    my $chnID = substr( $homolog, 4, 1 );
    my @atomResNums_homolog = @{ &AtomResNums4wholeSeq( $pdbID, $chnID ) };
    my ( $start_Q, $end_Q, $start_H, $end_H, $alignedSeq_Q, $alignedSeq_H ) =
      &searchHHpredFL( $blastAlignFL, $homolog );

    #-- get atomResNums for the aligned part of the query
    my @atomResNums_Q = ( 1 .. 9999 );
    my @atomResNums_aligned_Q =
      @{ &getAlignedAtomResNum( $start_Q, $alignedSeq_Q, \@atomResNums_Q ) };

    #-- get atomResNums for the aligned part of the homolog
    my @atomResNums_aligned_H =
      @{ &getAlignedAtomResNum( $start_H, $alignedSeq_H, \@atomResNums_homolog )
      };

    #--
    if ( scalar @atomResNums_aligned_Q ne scalar @atomResNums_aligned_H ) {
        die(
"Lengths of atomResNums_aligned_H and atomResNums_aligned_Q are not equal:$!"
        );
    }

    #-- write the aligned atomResNums for Q and H into a file
    unlink $outputFL if ( -e $outputFL );

    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    print OUTPUT "#Generated by genAtomResNumMapFL_new()\n";
    print OUTPUT
"#Note: we do not know the atomResNum for the query, so atomResNum (the 1st column of the output file) is seqResNum starting from 1.\n\n";
    print OUTPUT "QRY\t$homolog\n";

    for ( my $i = 0 ; $i < scalar @atomResNums_aligned_Q ; $i++ ) {
        print OUTPUT "$atomResNums_aligned_Q[$i]\t$atomResNums_aligned_H[$i]\n";
    }
    close OUTPUT;

    print "$outputFL generated.\n";
}

sub searchHHpredFL {

    my $hhpredFL  = shift @_;
    my $homologID = shift @_;    # '1ahjA,2'

    my ( $qryID, $totalLengthQry, $hhpredResults ) =
      &readHHpredOUTPUTFL_userProtDB($hhpredFL);

    my ( $homolog, $num ) = split( /,/, $homologID );
    my (
        $alignmentLength, $BitScore,   $Evalue,       $Identities,
        $Positives,       $similarity, $start_Q,      $end_Q,
        $start_S,         $end_S,      $localAlign_Q, $localAlign_S,
        $hhpred_prob
      )
      = (

        $hhpredResults->{$homolog}->{$num}->{'LAL'},
        $hhpredResults->{$homolog}->{$num}->{'bitscore'},
        $hhpredResults->{$homolog}->{$num}->{'evalue'},
        $hhpredResults->{$homolog}->{$num}->{'identityS'},
        $hhpredResults->{$homolog}->{$num}->{'positiveS'},
        $hhpredResults->{$homolog}->{$num}->{'similarity'},
        $hhpredResults->{$homolog}->{$num}->{'qstart'},
        $hhpredResults->{$homolog}->{$num}->{'qend'},
        $hhpredResults->{$homolog}->{$num}->{'sstart'},
        $hhpredResults->{$homolog}->{$num}->{'send'},
        $hhpredResults->{$homolog}->{$num}->{'qseq'},
        $hhpredResults->{$homolog}->{$num}->{'sseq'},
        $hhpredResults->{$homolog}->{$num}->{'hhpredProb'}
      );

    return ( $start_Q, $end_Q, $start_S, $end_S, $localAlign_Q, $localAlign_S );

}

sub AtomResNums4wholeSeq {

# read s2c file to return atomResNum mapped to a whole seq (including residues that are not in the ATOM section)
#
# s2c file format:
#		HEADSC 1m9x
#		COMMNT S2C correlation file created: Mon Apr 15 22:17:08 EDT 2013
#		COMMNT
#		COMMNT If you use this database, please cite:
#		COMMNT
#		COMMNT Guoli Wang, Jonathan W. Arthur, and Roland L. Dunbrack, Jr.
#		COMMNT "S2C: A database correlating sequence and atomic
#		COMMNT   coordinate numbering in the Protein Data Bank"
#		COMMNT   dunbrack.fccc.edu/Guoli/s2c
#		COMMNT   Copyright (c) February 2000, April 2002.
#		COMMNT
#		COMMNT SEQCRD columns are as follows:
#		COMMNT
#		COMMNT Column Positions Item
#		COMMNT      1       1-6 Record identifier
#		COMMNT      2         8 Chain
#		COMMNT      3        10 One letter residue code
#		COMMNT      4     12-14 SEQRES three letter residue code
#		COMMNT      5     16-18 ATOM three letter residue code
#		COMMNT      6     20-24 SEQRES residue number
#		COMMNT      7     26-31 ATOM residue number
#		COMMNT      8        33 PDB secondary structure
#		COMMNT      9        35 STRIDE secondary structure
#		COMMNT     10     37-43 Error flags
#		COMMNT
#		COMMNT Secondary structrue annotation:
#		COMMNT     H: Helix      E: Strand     T: Turn
#		COMMNT     B: Bridge     G: 310Helix   C: Coil
#		COMMNT
#		SEQCRD E M MET MET     1      1 C T 5
#		SEQCRD E V VAL VAL     2      2 C T 5
#		SEQCRD E N ASN ASN     3      3 C T 5
#		SEQCRD E P PRO PRO     4      4 C T 5
#		SEQCRD E T THR THR     5      5 E E -
#		SEQCRD E V VAL VAL     6      6 E E -
#		SEQCRD E F PHE PHE     7      7 E E -
#		SEQCRD E F PHE PHE     8      8 E E -
#		SEQCRD E D ASP ASP     9      9 E E -
#		SEQCRD E I ILE ILE    10     10 E E -
#		SEQCRD E A ALA ALA    11     11 E E -

    our $s2cDIR;
    my $pdbID = lc( shift @_ );
    my $chnID = shift @_;

    my $s2cFL = "$s2cDIR/$pdbID.sc";

    print
"Read s2c file for atomResNums mapped to the whole seq: $pdbID$chnID...\n";

    my @atomResNums_mapped2Seq;

    open( INPUT, "<$s2cFL" ) or die("Cannot open $s2cFL:$!");
    while (<INPUT>) {
        s/[\n\r]//gm;

        if (/^SEQCRD\s+$chnID/) {

            my @tmp = split( /\s+/, $_ );
            my $atomResNum = $tmp[6];

            push @atomResNums_mapped2Seq, $atomResNum;

        }
    }
    close INPUT;

    if ( !@atomResNums_mapped2Seq ) {
        die("Nothing read from $s2cFL:$!");
    }

    return \@atomResNums_mapped2Seq;

}

sub getAlignedAtomResNum {
    use strict;
    my $start       = shift @_;         #-- from the blast output
    my $alignedSeq  = shift @_;         #-- from the blast output
    my @atomResNums = @{ shift @_ };    #for the whole protein

    my @alignedSeq_array = split( //, $alignedSeq );
    my $LAL = length($alignedSeq);
    my @atomResNum_aligned;
    my $index = $start - 1;
    for ( my $i = 0 ; $i < $LAL ; $i++ ) {
        if ( $alignedSeq_array[$i] ne '-' ) {

            if ( !defined $atomResNums[$index] ) {
                die(
"index $index does not have elements defined in @atomResNums:$!"
                );
            }

            push @atomResNum_aligned, $atomResNums[$index];
            $index = $index + 1;

        }
        else {

            push @atomResNum_aligned, '-';
        }
    }

    if ( !@atomResNum_aligned ) {
        die(
"\n**Error!Aligned atomResNum not defined!\n start=$start\naligned seq: $alignedSeq:$!"
        );
    }

    return \@atomResNum_aligned;

}

sub belong_new {

#The difference between belong() and belong_new() is that @homointerologsToBeDel in belong_new() may contain wild cards

    use strict;

    my $pair = shift @_;    # homolog pair: A':B'
    my @homointerologsToBeDel =
      @{ shift @_ };        #(1okbB:*,1lqjC:*,3fclB:*,1lqmA:*)
    my $ans = 0;

    #print LOG "Whether $pair belongs to @homointerologsToBeDel...\n";

    #---------- convert pdb ID of $pair to lower case
    my ( $a, $b ) = split( /:/, $pair );
    my $pdbID = lc( substr( $a, 0, 4 ) );
    my $chn = substr( $a, 4, 1 );
    $a = "$pdbID$chn";

    $pdbID = substr( $b, 0, 4 );
    $chn   = substr( $b, 4, 1 );
    $b     = "$pdbID$chn";

    $pair = "$a:$b";

    #---------------------------------------------------

    foreach my $_ (@homointerologsToBeDel) {

        #		case 0: $_ = '1okbB:1okbC'
        #		case 1: $_ = '1okbB:*', '1okbB:1okb*'
        #		case 2: $_ = '*:1lqmH', '1ahj*:1ahjB'
        #		case 3: $_ = '1ahj*', '1ahj*:1ahj*'

   #In perl,  the star (*) means to match the preceding item zero or more times.
   #print LOG "\npair: $pair\n";
   #print LOG "homointerologToBeDel: $_\n";
   #print LOG "ans = $ans (0:no. 1:yes)\n";

        s/\s+//g;

        if (/\w{5}:\w{5}/) {

            #		case 0: $_ = '1okbB:1okbC'
            $ans = &belong_case0( $pair, $_ );
        }
        elsif ( /\w{5}:\*$/ || /\w{5}:\w{4}\*$/ ) {

            #		case 1: $_ = '1okbB:*'
            $ans = &belong_case1( $pair, $_ );
        }
        elsif ( /^\*:\w{5}/ || /\w{4}\*:\w{5}/ ) {

            #		case 2: $_ = '*:1lqmH', '1ahj*:1ahjB'
            $ans = &belong_case2( $pair, $_ );
        }
        elsif ( /^\w{4}\*$/ || /^\w{4}\*:\w{4}\*$/  ) {

            #		case 3: $_ = '1ahj*','1ahj*:1ahj*'
            $ans = &belong_case3( $pair, $_ );
        }
        else {
            die("templates to be deleted ($_) format is wrong:$!");
        }

    #--- as long as $pair is equal to ONE homolog-pair-to-be-del, return $ans= 1
        if ( $ans == 1 ) {
            return $ans;
        }
    }

    return $ans;

}

sub belong_case0 {

    #		case 0: $homointerologToBeDel = '1okbB:1okbC'

    my $pair =
      shift @_;    # the pdb ID of the pair is already converted to lower case
    my $homointerologToBeDel = shift
      @_; #-- pdb ID of $homointerologToBeDel is already converted to lower case
    my $ans = 0;

#   print "case0: pair = $pair, homointerologToBeDel = $homointerologToBeDel\n";
    if ( $pair eq $homointerologToBeDel ) {

        $ans = 1;
    }

    return $ans;

}

sub belong_case1 {

    #		case 1: $homointerologToBeDel = '1okbB:*'
    my ( $pair, $homointerologToBeDel ) = @_;
    my $ans = 0;

#   print "case1: pair = $pair, homointerologToBeDel = $homointerologToBeDel\n";

    my ( $A_del, $B_del ) = split( /:/, $homointerologToBeDel );
    my ( $a,     $b )     = split( /:/, $pair );
    if ( $a =~ /$A_del/ ) {

        # $a is the first protein of $pair
        $ans = 1;
    }

    return $ans;
}

sub belong_case2 {

    #		case 2: $homointerologToBeDel = '*:1lqmH'

    my ( $pair, $homointerologToBeDel ) = @_;
    my $ans = 0;

#   print "case2: pair = $pair, homointerologToBeDel = $homointerologToBeDel\n";

    my ( $A_del, $B_del ) = split( /:/, $homointerologToBeDel );
    my ( $a,     $b )     = split( /:/, $pair );
    if ( $b =~ /$B_del/ ) {

        # $a is the first protein of $pair
        $ans = 1;
    }

    return $ans;
}

sub belong_case3 {

    #		case 3: $homointerologToBeDel = '1ahj*'

    my ( $pair, $homointerologToBeDel ) = @_;
    my $ans = 0;

#    print "case3: pair = $pair, homointerologToBeDel = $homointerologToBeDel\n";

    my $pdbID_pair = substr( $pair, 0, 4 );
    my $pdbID_del = substr( $homointerologToBeDel, 0, 4 );

    if ( $pdbID_pair eq $pdbID_del ) {
        $ans = 1;
    }

    return $ans;
}

sub collectPredictionFLsIntoOneFL1 {

    #put all $taskID.KNN.prediction files into one output file
    #format is usr friendly.

    use strict;
    use File::Basename;

    our $logFL;

    my $jobDIR   = shift @_;
    my @QRYpairs = @{ shift @_ };    #'A:B' , '1fevA:1fevB'
    my $outputFL = shift @_;         #final prediction to show usrs
    my ( $startTime, $endTime, $timeUsed ) = @_;

    my $predictFL;
    my $taskID;
    my $taskID_partner;
    our $safe_filename_characters;

    print LOG "\n\n\n\nCollecting PredictionFLs Into One FL...\n";

    &finalPredictionFLheader( $outputFL, $startTime, $endTime, $timeUsed );

    #process each task ID

    foreach my $QRYpair (@QRYpairs) {

        #'A:B'

        open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
        print OUTPUT "\n\n***\tChain Pair: $QRYpair\t***\n";
        close OUTPUT;

        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $jobDIR );

        #1. write A|A:B
        $taskID         = $c1;
        $taskID_partner = $c2;
        $predictFL =
          "$pairDIR/$taskID/predictionResults/$taskID.KNN.prediction";

        &writeFinalPredictionFL( $taskID, $taskID_partner, $predictFL,
            $outputFL );

        #2. write B|B:A

        $taskID         = $c2;
        $taskID_partner = $c1;
        $predictFL =
          "$pairDIR/$taskID/predictionResults/$taskID.KNN.prediction";

        &writeFinalPredictionFL( $taskID, $taskID_partner, $predictFL,
            $outputFL );

    }

    print LOG
"All predictions for proteins in $jobDIR are put into $outputFL(usr-friendly format).\n ";

}

sub writeFinalPredictionFL {

    #write the interface prediction file for users
    #Format:
    #residual, resSeq, pInt, pInt-score

    use strict;

    my ( $taskID, $taskID_partner, $predictFL, $outputFL ) = @_;

    if ( -e $predictFL ) {
        my ( $seq, $predictionScore, $pInt, $mode ) =
          &ReadPredictionFL($predictFL);

        open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
        print OUTPUT "\n\nQry = $taskID|$taskID:$taskID_partner\tMODE: $mode\n";
        print OUTPUT "SEQ\tpINT\tSCORE\n";

        for ( my $i = 0 ; $i < scalar @$seq ; $i++ ) {
            print OUTPUT "$seq->[$i]\t$pInt->[$i]\t$predictionScore->[$i]\n";
        }
        close OUTPUT;
    }
    else {

        open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
        print OUTPUT
"\n\nQry = $taskID|$taskID:$taskID_partner\tMODE: No Prediction Made.\n";

        close OUTPUT;

    }

}

sub collectPredictionFLsIntoOneFL2 {

    #put all $taskID.KNN.prediction files into one output file
    #format is machine friendly.

    use strict;
    use File::Basename;

    my $jobDIR   = shift @_;
    my @QRYpairs = @{ shift @_ };    #'A:B' , '1fevA:1fevB'
    my $outputFL = shift @_;         #final prediction to show usrs
    my ( $startTime, $endTime, $timeUsed ) = @_;

    my $predictFL;
    my $taskID;
    my $taskID_partner;

    print LOG "Collecting PredictionFLs Into One FL...\n";

    &finalPredictionFLheader( $outputFL, $startTime, $endTime, $timeUsed );

    #process each task ID
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");

    foreach my $QRYpair (@QRYpairs) {

        #'A:B'

        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $jobDIR );

        if ( !-d $pairDIR ) {
            print(
"pairDIR $pairDIR does not exist. No prediction made. Next pair.\n"
            );
            next;

        }

        #1. write A|A:B
        $taskID         = $c1;
        $taskID_partner = $c2;
        $predictFL =
          "$pairDIR/$taskID/predictionResults/$taskID.KNN.prediction";

        open( PRED, "<$predictFL" ) || die("Cannot open $predictFL:$!");
        foreach (<PRED>) {
            s/[\n\r]//mg;
            if (/^([\s\t]{0,})>[^\|]+\|([\w\s]+)/) {

                #header line
                #				 >B|TwilightMode1
                #>B|No homolog found

                my $mode = $2;
                print OUTPUT "$1>$taskID|$taskID:$taskID_partner|$mode\n";
            }
            else {
                print OUTPUT "$_\n";
            }
        }
        close PRED;

        #2. write B|B:A

        $taskID         = $c2;
        $taskID_partner = $c1;
        $predictFL =
          "$pairDIR/$taskID/predictionResults/$taskID.KNN.prediction";

        open( PRED, "<$predictFL" ) || die("Cannot open $predictFL:$!");
        foreach (<PRED>) {
            s/[\n\r]//mg;
            if (/^([\s\t]{0,})>.+\|([\w\s]+)/) {

                #header line
                my $mode = $2;
                print OUTPUT "$1>$taskID|$taskID:$taskID_partner|$mode\n";
            }
            else {
                print OUTPUT "$_\n";
            }
        }
        close PRED;

    }
    close OUTPUT;

    print LOG
"All predictions for proteins in $jobDIR are put into $outputFL (machine-friendly format).\n ";

}

sub collectPredictedResiPairFLsIntoOneFolder {
    my ( $jobDIR, $QRYpairs, $outputCACA_DIR, $startTime, $endTime, $timeUsed )
      = @_;
    # add header to CA-CA files and save them to $outputDIR
    &addHeader2CaCaFLs( $jobDIR, $QRYpairs, $startTime, $endTime, $timeUsed );
    system("tar -czf $outputCACA_DIR.tar.gz -C $outputCACA_DIR .") == 0
      or die("compress final ca-ca dir failed: $outputCACA_DIR:$!");
    unlink "$jobDIR.tar.gz" if ( -e "$jobDIR.tar.gz" );
    system("tar -czf $jobDIR.tar.gz -C $jobDIR . ") == 0
      or die("compress final ca-ca dir failed: $jobDIR:$!");

}

#sub collectPredictedResiPairFLsIntoOneFolder_old {
#
## move all the predicted residue-residue distance files for each query chain pairs into one folder
##
#    use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
#
#    my $jobDIR   = shift @_;
#    my @QRYpairs = @{ shift @_ };    #'A:B' , '1fevA:1fevB'
#
#    my $outputResiPairDIR = shift @_;
#    mkdir $outputResiPairDIR if ( !-d $outputResiPairDIR );
#
#    print LOG "\n\n\n\nCollecting Ca-Ca-distance files Into One Folder ...\n";
#
#    foreach my $QRYpair (@QRYpairs) {
#        my ( $pairDIR, $c1, $c2 ) = &getPairDIR2( $QRYpair, $jobDIR );
#        my $resiPairDIR = "$pairDIR/resiPair";
#        if ( !-d $resiPairDIR ) {
#            next;
#        }
#        my $predictedResiPairFL =
#          "$resiPairDIR/Ca_Ca_distance*.txt";    #"Ca_Ca_distance.A_B.txt"
#        system("cp  $predictedResiPairFL $outputResiPairDIR") == 0
#          or die("Cannot copy $predictedResiPairFL to $outputResiPairDIR:$!");
#
#        system("cp -r $resiPairDIR $outputResiPairDIR/resiPair.$c1\_$c2") == 0
#          or die("Cannot copy $resiPairDIR to $outputResiPairDIR:$!");
#    }
#
#    #--compress the folder
#    use Cwd;
#    my $cwd = getcwd;
#    chdir($outputResiPairDIR) or die("Cannot enter $outputResiPairDIR:$!");
#
#    #- check whether this folder is empty
#    #
#    my $num_files = `ls |wc -l`;
#    if ( $num_files =~ /^0/ ) {
#        print LOG "$outputResiPairDIR is empty.\n";
#        return;
#    }
#
##	bzip2 "<$outputResiPairDIR/*>"  => "predicted_resiPair_distance.bz2" or die ("bizp2 failed:$Bzip2Error\n");
#    system("tar -czf predicted_resiPair_distance.tar.gz *") == 0
#      or die("Cannot compress $outputResiPairDIR:$!");
#    chdir($cwd);
#    print LOG "tar -czf predicted_resiPair_distance.tar.gz .\n";
#
#    print LOG
#"All predicted residue pair Ca distances for proteins in $jobDIR are put into $outputResiPairDIR.\n\n\n";
#}
#
sub writeArray2FL {
    my @a        = @{ shift @_ };
    my $outputFL = shift @_;

    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    foreach (@a) {

        print OUTPUT "$_\n";
    }
    close OUTPUT;

    print LOG "$outputFL generated.\n";

}

sub addPredCC {

    #- predict CC for each template and add the pred_CC to $ori_statFL

    use globalVariables;
    our $Rscript;

    #-- predict CC: mick's model
    #        my $intercept = -0.53162;
    #        my $a         = 0.00051917;
    #        my $b         = 0.00517;
    #        my $c         = 0.60002;
    #        my $d         = 0.08854;
    #        my $similarityScore =
    #          $intercept +
    #          $a * $mean_logEval +
    #          $b * $mean_PositiveS +
    #          $c * $aligLen_Qry * $aligLen_Hom +
    #          $d * log($alignmentLen);
    #
    print LOG "\n---------------------------------\n";
    print LOG "CALL RF MODEL TO PREDICT CC ...";
    print LOG "\n---------------------------------\n";

    #-- random forest
    my $statFL            = shift @_;
    my $statFL_withPredCC = "$statFL.tmp";    #- output of predCC.R
    my $command = "$Rscript predCC.R $statFL $statFL_withPredCC";
    system($command) == 0 or die("FAILED: $command:$!");

    #-- reformat $statFL_withPredCC into the format of $statFL
    &reformat_Routput( $statFL, $statFL_withPredCC );

    #-- clean up
    unlink "$statFL.tmp";

    #--
    print LOG "\nPredicted CC added to $statFL\n\n";

}

sub reformat_Routput {

# original format:
#
#    "qryPairID"	"homologPair"	"pred_CC"	"len_Qry1"	"len_Qry2"	"len_H1"	"len_H2"	"numInt_Q1"	"numInt_Q2"	"numInt_Homolog1"	"numInt_Homolog2"	"EVal1"	"EVal2"	"SID1"	"SID2"	"Positives1"	"Positives2"	"Similarity1"	"Similarity2"	"hhpredProb1"	"hhpredProb2"	"hhpred_pVal1"	"hhpred_pVal2"	"hhpred_SS1"	"hhpred_SS2"	"LAL_Q1"	"LAL_H1"	"LAL_Q2"	"LAL_H2"	"frac_LAL1"	"frac_LAL2"	"TP_1"	"TN_1"	"FP_1"	"FN_1"	"CC_1"	"specificity_1"	"sensitivity_1"	"TP_2"	"TN_2"	"FP_2"	"FN_2"	"CC_2"	"specificity_2"	"sensitivity_2"
#    "1acbE:1acbI"	"4b2aC,452:4b2aD,8"	0.727130813492064	245	70	223	66	23	15	28	17	3.5e-41	3.3e-28	40.6	98.4	73.5	98.4	1.7	4.6	100	99.9	4.4e-45	5e-32	26.3	6.6	229	222	66	66	0.93	0.94	18	186	9	5	0.69	0.67	0.78	12	46	2	3	0.78	0.86	0.8
#
#    This function 1) removes quotation marks, and 2) add headers.

    my $statFL_reference       = shift @_;
    my $statFL_to_be_formatted = shift @_;

    #-- remove quotation marks
    my $command =
"perl -ple 's/\"//g;' $statFL_to_be_formatted > $statFL_to_be_formatted.tmp";
    system($command) == 0 or die("FAILED: $command:$!");

    #-- add header lines
    $command = "egrep '^#' $statFL_reference > $statFL_reference.header ";
    system($command) == 0 or die("FAILED: $command:$!");
    $command =
"cat $statFL_reference.header $statFL_to_be_formatted.tmp > $statFL_reference";
    system($command) == 0 or die("FAILED: $command:$!");

    #-- clean up
    unlink "$statFL_reference.header";
    unlink "$statFL_to_be_formatted.tmp";

}

#----------------------------------------------------
sub hasTemplate {
    my $qryPairDIR = shift @_;
    my $templates_used_in_predictionFL =
      "$qryPairDIR/templates_used_in_prediction.lst";
    my $ans = 1;

    my $predCCs =
      &readTemplates_used_in_predictionFL($templates_used_in_predictionFL);
    if ( !defined $predCCs ) {
        $ans = 0;
    }
    return $ans;

}

sub readTemplates_used_in_predictionFL {

    # INPUT: templates_used_in_prediction.lst
    #
    #    Template_ID  predicted_CC
    #    1m9xB,50:1m9xC,9    0.998
    #    4lqwA,1:4lqwD,77    0.986

    my $homcomplexes_used_in_predictionFL = shift @_;
    my $predCCs;

    open( INPUT, "<$homcomplexes_used_in_predictionFL" )
      or die("Cannot open $homcomplexes_used_in_predictionFL:$!");
    while (<INPUT>) {
        s/[\n\r]//gm;

        if (/^\w+/) {
            my ( $templateID, $pred_CC ) = split( /\s+/, $_ );
            $predCCs->{$templateID} = $pred_CC;
        }
    }
    close INPUT;

    #  if (!defined $predCCs){
    #      die("nothing read from $homcomplexes_used_in_predictionFL:$!");
    #  }

    return $predCCs;
}


sub unique_keepOrder{
    #-- keep the unique elements and keep the original order

    my @a = @{shift @_};
    my @b;
    my %seen;

    foreach(@a){

        if (! defined $seen{$_} ){
            push @b, $_;
            $seen{$_}=1;
        }

    }

    return \@b;

}
1;

