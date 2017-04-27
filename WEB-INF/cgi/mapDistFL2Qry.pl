#!/usr/bin/env perl
#===============================================================================
#
#         FILE: mapDistFL2Qry.pl
#
#        USAGE: ./mapDistFL2Qry.pl
#
#        INPUT 1 (Ca-Ca distance file on template 1. The structure of Template 1 was used to generate superimposed models. ChnID1 and chnID2 are always A and B, because chain IDs of template 1 are set to A and B):
#
#             ChnID1 aa1 ChnID2 aa2 mean min max
#             A 128 B 90 14.437 13.6659 15.1642
#             A 128 B 92 12.270 12.163 12.4211
#
#        INPUT 2 (two maps between qry seqResNum and Template 1 atomResNum):
#
#              QRY    3bn1A
#              2   7
#              3   8
#
#        INPUT 3 (ID map between qry IDs and template1 chain IDs):
#
#              QRY superimposed_models
#              C  A
#              D  B
#
#        OUTPUT:
#
#            chnID_Qry1  aa_Qry1 chnID_Qry2  aa_Qry2 A_template1 aa1_template1 B_template1 aa2_template1 mean    min max
#            P   89  Q   103 A 99 B 104 13.254 12.0263 14.5945
#            P   89  Q   102 A 99 B 103 13.036 12.2848 14.1158
#            P   89  Q   104 A 99 B 105 16.454 14.9065 17.8477
#
#      CREATED: 02/05/2015 12:01:19 PM
#===============================================================================

use strict;
use warnings;
use utf8;
use File::Basename;
use File::Copy;
use File::Spec::Unix;

my $PSHomPPI_result_onePair_DIR = shift @ARGV;
#  '../../uploadData/Thu_Feb__5_17_23_01_2015.72.test/C:D';

my ( $template_pdbID_rec, $template_pdbID_lig ) =
  &getTemplate1ChnIDs($PSHomPPI_result_onePair_DIR);    #2x83B

#------
my $qry_template1_chnIDmap =
  "$PSHomPPI_result_onePair_DIR/superimposed_models/IDmap_qry_supModels.txt";
my $chnMap = &getQryTemplateChnMap($qry_template1_chnIDmap)
  ;    #chnMap->{A}=qry_chnID1; chnMap->{B}=qry_chnID2
my $qryChn1 = $chnMap->{'A'};    #template chain IDs are always A and B
my $qryChn2 = $chnMap->{'B'};    #template chain IDs are always A and B
my $resiNumMapFL_A =
"$PSHomPPI_result_onePair_DIR/TemplatesUsed/qry_template1_resNumMap/$qryChn1\_$template_pdbID_rec.resiNumMap"
  ;                              # C_2x83B.resiNumMap
my $resiNumMapFL_B =
"$PSHomPPI_result_onePair_DIR/TemplatesUsed/qry_template1_resNumMap/$qryChn2\_$template_pdbID_lig.resiNumMap"
  ;                              # C_2x83B.resiNumMap

print "\nResidue map between qry and template 1: \n1. $resiNumMapFL_A \n2. $resiNumMapFL_B\n\n";

# $resiNumMapFL_A->{template}=qry

my $resiNumMap_A = &readResiNumMapFL($resiNumMapFL_A)
  ;    # $resiNumMap_A->{template_chnA_atomResNum}=Qry_chnRec_seqResNum
my $resiNumMap_B = &readResiNumMapFL($resiNumMapFL_B)
  ;    # $resiNumMap_B->{template_chnB_atomResNum}=Qry_chnLig_seqResNum

#-------
my @Ca_Ca_FLs_ori =
`ls $PSHomPPI_result_onePair_DIR/superimposed_models/Ca_Ca_distances/*_Ca_Ca_distance*txt`;
if ( !@Ca_Ca_FLs_ori ) {
    die(
"No Ca-Ca distance files in $PSHomPPI_result_onePair_DIR/superimposed_models/Ca_Ca_distances:$!"
    );
}
print
"\n\nAdd Qry atomResNums to Ca-Ca distance files under superimposed_models/Ca_Ca_distances\n\n";

foreach my $Ca_Ca_FL_ori (@Ca_Ca_FLs_ori) {

#my $Ca_Ca_FL_ori = "$PSHomPPI_result_onePair_DIR/superimposed_models/Ca_Ca_restraints/cluster1.Ca_Ca_distance.txt";

    $Ca_Ca_FL_ori =~ s/[\r\n]//mg;

    my $rel_path_CACA_FL_ori = File::Spec::Unix->abs2rel($Ca_Ca_FL_ori);

    print "\n\n**** original CA-CA file (w/o qry atomResNum) $rel_path_CACA_FL_ori\n";
    my $Ca_Ca_FL_tmp = "$Ca_Ca_FL_ori.tmp";
    &header($Ca_Ca_FL_tmp);

    my $flag               = 0;
    my $num_contacts_final = 0;
    open( OUTPUT, ">>$Ca_Ca_FL_tmp" ) or die("Cannot open $Ca_Ca_FL_tmp:$!");
    open( INPUT,  "<$Ca_Ca_FL_ori" )  or die("Cannot open $Ca_Ca_FL_ori:$!");
    while (<INPUT>) {
        s/[\n\r]//gm;
        if (/^(#|!)/ ) {
            print OUTPUT "$_\n";
            next;
        }
        if (/^ChnID1\s+aa1/) {
#             ChnID1 aa1 ChnID2 aa2 mean min max
            $flag++;
            print OUTPUT "chnID_Qry1\taa_Qry1\tchnID_Qry2\taa_Qry2\tA_template1\taa1_template1\tB_template1\taa2_template1\tmean\tmin\tmax\n";
            next;
        }
        if (/^\w+\s+[\w\-]+\s+/) {
            #A 99 B 104 13.254 12.0263 14.5945
#            print "input: $_\n";

            $flag ++;

            #-- check input format
            my @tmp=split(/\s+/,$_);
            if (scalar @tmp != 7){
                die ("Input format wrong. Check $Ca_Ca_FL_ori:$!");
            }

            #--
            my ($chnID1, $atomResNum_templateA, $chnID2, $atomResNum_templateB, $mean, $min,
                $max )
              = @tmp;
            my $seqResNum_qry1 = $resiNumMap_A->{$atomResNum_templateA};
            my $seqResNum_qry2 = $resiNumMap_B->{$atomResNum_templateB};

            if ( !defined $seqResNum_qry1 || !defined $seqResNum_qry2 ) {

               print "WARNING: seqResNum for queries not defined.cacaFL_ori: $Ca_Ca_FL_ori ; atomResNum_template A: $atomResNum_templateA ; atomResNum_template_B: $atomResNum_templateB !!\n";
                next;
            }
            if ( $seqResNum_qry1 eq '-' || $seqResNum_qry2 eq '-' ) {
                next;
            }

            print OUTPUT "$qryChn1\t$seqResNum_qry1\t$qryChn2\t$seqResNum_qry2\t$_\n";
            $num_contacts_final++;

        }
    }
    close INPUT;
    close OUTPUT;

    if ( $flag < 2 ) {

        print "\n\n::WARNING: Nothing read from the Ca-Ca file. The format might be wrorng or the file is empty:  $rel_path_CACA_FL_ori. Next ca-ca file.\n\n";

        #-- label the ca-ca filename as empty
        move($Ca_Ca_FL_ori,"$Ca_Ca_FL_ori"."_emtpy");
        unlink ($Ca_Ca_FL_tmp);

        #-- unlink the tbl file generated from the $Ca_Ca_FL_ori
        my $dirname = dirname ($Ca_Ca_FL_ori);
        my @tmp = split(/_/,basename($Ca_Ca_FL_ori)); #$Ca_Ca_FL_ori ='dir/cluster1_Ca_Ca_distance.txt'
        my $clusterID = shift @tmp;
        my $haddockRestFL = "$dirname/$clusterID\_restraints.tbl";
        unlink($haddockRestFL) if (-e $haddockRestFL);

        #-- go to next CA-CA file
        next;

    }

    my $minimal_contactNum = 5;
    if ( $num_contacts_final < $minimal_contactNum ) {
        print
"WARNING: The final Ca-Ca file  ($rel_path_CACA_FL_ori) with Qry atomResNums has < $minimal_contactNum contact pairs. This Ca-Ca file is not reliable. Lable it as notUsedAsFinal.\n";
        unlink $Ca_Ca_FL_ori;
        move($Ca_Ca_FL_tmp, "$Ca_Ca_FL_ori".".notUsedAsFinal") or die ("FAILed:$!");
    }
    else {

        move( $Ca_Ca_FL_tmp, $Ca_Ca_FL_ori ) or die("Failed:$!");
        print
"query atomResNums added to Ca-Ca distance file: $rel_path_CACA_FL_ori\n";
    }

}

#--------
sub header {
    my $ca_ca_fl_new = shift @_;
    unlink $ca_ca_fl_new if ( -e $ca_ca_fl_new );
    open( OUTPUT, ">>$ca_ca_fl_new" ) or die("cannot open $ca_ca_fl_new:$!");
    print OUTPUT "# generated by $0\n";
    close OUTPUT ;

}

sub getQryTemplateChnMap {

    #all superimposed models and template 1 uses a and b as chain ids.
    #
    #              qry superimposed_models
    #              c  a
    #              d  b
    my $qry_template1_chnidmap = shift @_;
    my $chnmap;
    open( INPUT, "<$qry_template1_chnidmap" )
      or die("cannot open $qry_template1_chnidmap:$!");
    while (<INPUT>) {
        s/[\r\n]//mg;

        if (/^qry/) {

            #              qry superimposed_models
            next;
        }
        if (/^\w+\s+\w+/) {
            my ( $qry_chnid, $superimposed_model_chnid ) = split( /\s+/, $_ );
            $chnmap->{$superimposed_model_chnid} = $qry_chnid;

        }

    }
    close(INPUT);

    if ( !defined $chnmap ) {
        die("nothing read from $qry_template1_chnidmap:$!");
    }
    return $chnmap;

}

sub readResiNumMapFL {

    #              qry    3bn1a
    #              2   7
    #              3   8
    my $resiNumMapFL = shift @_;
    my $map;

    open( INPUT, "<$resiNumMapFL" ) or die("cannot open $resiNumMapFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^qry/) {
            next;
        }
        if (/^[\w\-]+\s+[\w\-]+/) {
            my ( $atomResNum_qry, $atomresnum_template ) = split( /\s+/, $_ );
            $map->{$atomresnum_template} = $atomResNum_qry;

        }
    }
    close INPUT;

    if ( !defined $map ) {
        die("nothing read from $resiNumMapFL:$!");
    }
    return $map;
}

sub getTemplate1ChnIDs {
    #return the pdb ID and chain IDs for template 1

# The format of hom-complexes_used_in_prediction.lst:
#
# qryPairID  homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2  EVal1   EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2   LAL1    LAL2    frac_LAL1   frac_LAL2
# P:Q 1ak4D,80:1ak4A,16   137 164 145 165 Nan Nan 7   18  1.5E-39 9.6E-47 98.5    100.0   100.0   100.0   4.4 4.7 100.0    100.0   135 163 0.92    0.98


    my $PSHomPPI_result_onePair_DIR = shift @_;

    #--return values
    my $template_rec ;    #1ak4A
    my $template_lig ;

    #-- input file
    my $homComplexes_used_in_predictionFL =
"$PSHomPPI_result_onePair_DIR/TemplatesUsed/templates_used_in_prediction.stat";
    open( INPUT, "<$homComplexes_used_in_predictionFL" )
      or die("Cannot open $homComplexes_used_in_predictionFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^\w+/) {

        # P:Q 1ak4D,80:1ak4A,16   137 164 145 165

            #read template 1 info
            my @a = split( /\s+/, $_ );
            my $template1 = $a[1]; # -- 1ak4D,80:1ak4A,16
             ($template_rec, my $group1, $template_lig, my $group2) = split(/[:,]/, $template1);
            last;
        }

    }
    close(INPUT);

    if ( !defined $template_rec || !defined $template_lig ) {
        die(
"Template 1 info NOT read from $homComplexes_used_in_predictionFL:$!"
        );
    }

    return ( $template_rec, $template_lig );

}
sub getTemplate1ChnIDs_ori {

#    #pdbID  recA:ligB   local_IdentityS_AA' local_IdentityS_BB' global_IdentityS_AA'    global_IdentityS_BB'
#    1ak4    A:D 100 98.52   99.6969696969697    94.4036345330984
#    1ak4    B:C 100 98.52   99.6969696969697    94.4036345330984

    #return the pdb ID and chain IDs for template 1

    my $PSHomPPI_result_onePair_DIR = shift @_;
    my $template_pdbid;
    my $template_rec_lig;
    my $homComplexes_used_in_predictionFL =
"$PSHomPPI_result_onePair_DIR/TemplatesUsed/hom-complexes_used_in_prediction.lst";
    open( INPUT, "<$homComplexes_used_in_predictionFL" )
      or die("Cannot open $homComplexes_used_in_predictionFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^\w+/) {

            #    1ak4    A:D 100 98.52   99.6969696969697    94.4036345330984

            #read template 1 info
            my @a = split( /\s+/, $_ );
            $template_pdbid   = $a[0];    #1ak4
            $template_rec_lig = $a[1];    #A:D
            last;
        }

    }
    close(INPUT);

    if ( !$template_rec_lig || !defined $template_pdbid ) {
        die(
"Template 1 info NOT read from $homComplexes_used_in_predictionFL:$!"
        );
    }
    my ( $rec_ID, $lig_ID ) = split( /:/, $template_rec_lig );
    my $template_rec = "$template_pdbid$rec_ID";    #1ak4A
    my $template_lig = "$template_pdbid$lig_ID";

    return ( $template_rec, $template_lig );

}
