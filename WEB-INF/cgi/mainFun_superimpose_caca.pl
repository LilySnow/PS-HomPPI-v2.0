#!/usr/bin/perl -w
# Li Xue
# July 3rd, 2014
#
# 0) generate superimposed models
# 1) cluster the superimposed models
# 2) calculate Ca-Ca for each cluster
# 3) write haddock restraint files
#

use strict;
use File::Copy;
use File::Basename;
use PSHomPPI_resiPairs;

my $caDistThr = shift @ARGV
  ; #15 #-- only include contacts with CA-CA distance <=8.5 angstroms into the final output file
my $PSHomPPI_result_onePair_DIR = shift @ARGV
  ; #PSHomPPI_result_onePair_DIR=../../uploadData/Thu_Feb__5_17_23_01_2015.72.test/C:D
my $supUnbound2TemplatePDB_DIR =
  "$PSHomPPI_result_onePair_DIR/superimposed_models/supUnbound2TemplatePDB";


#------------------------------------
# 0. -- check whether this query pair has templates

if (! &hasTemplate($PSHomPPI_result_onePair_DIR)){
    my $qrypair = basename ($PSHomPPI_result_onePair_DIR);
    print "Qry pair ($qrypair) has no templates in Safe/Twi/Dark Zones. No prediction is made.\n";
    exit;
}



#------------------------------------
# 0. -- generate superimposed models

my $skip = 0;

if ( $skip ne 1 ) {

    printf "\n\n======================================================";
    printf "\nGenerating superimposed models ...........................\n";
    printf "======================================================\n\n";

    my $command =
"perl genSuperimposePDBFL.pl  $PSHomPPI_result_onePair_DIR > $PSHomPPI_result_onePair_DIR/log_superimposition.txt";

    print "COMMAND: $command\n";
    system("$command") == 0 or die("FAILED: $command :$!");

    print
"Finish generating superimposed models. LOG file: $PSHomPPI_result_onePair_DIR/log_superimposition.txt\n";

}

#------------------------------------
# 0. -- copy the final superimposed models from different templates into one folder called "final"

$skip = 0;
if ( $skip ne 1 ) {

    printf "\n\n======================================================\n";
    printf
"Copy the final superimposed models from different templates into one folder called final ......\n";
    printf "======================================================\n\n";

    if ( !-d $supUnbound2TemplatePDB_DIR ) {

        die("ERROR: $supUnbound2TemplatePDB_DIR does not exit:$!");
    }

    if ( !-d "$supUnbound2TemplatePDB_DIR/final" ) {
        mkdir "$supUnbound2TemplatePDB_DIR/final";
    }

    if ( !-d "$supUnbound2TemplatePDB_DIR/final/template_pdbs" ) {
        mkdir "$supUnbound2TemplatePDB_DIR/final/template_pdbs";
    }

    foreach my $j (`ls $supUnbound2TemplatePDB_DIR/template*/finalSup*`) {

        $j =~ s/[\n\r]//mg;
        my $filename   = basename ($j, ".pdb");
        my $templateID = dirname ($j);
        $templateID = basename ($templateID);


        #--
        copy( $j,
            "$supUnbound2TemplatePDB_DIR/final/$filename.$templateID.pdb" );

        #--copy alignment symbol files to 'final'
        copy(
"$supUnbound2TemplatePDB_DIR/$templateID/protein1.sup.aligned_resiNum",
"$supUnbound2TemplatePDB_DIR/final/protein1.$templateID.aligned_resiNum"
        );
        copy(
"$supUnbound2TemplatePDB_DIR/$templateID/protein2.sup.aligned_resiNum",
"$supUnbound2TemplatePDB_DIR/final/protein2.$templateID.aligned_resiNum"
        );

        #--copy template pdb files to 'final/template_pdb'
        my $command =
          "echo \'$filename\'| perl -ne '\@a=split(/\\./, \$_); END{print \$a[1]}'";
        my $template_pdbIDchnIDs = `$command`;
        $template_pdbIDchnIDs=~ s/[\n\r]//g;

        if (! defined $template_pdbIDchnIDs || $template_pdbIDchnIDs =~/^\s*$/){
            die("template_pdbIDchnIDs are not defined. COMMAND: $command:$!");
        }


        $command =
"cat $supUnbound2TemplatePDB_DIR/$templateID/$templateID.rec*pdb $supUnbound2TemplatePDB_DIR/$templateID/$templateID.lig*pdb > $supUnbound2TemplatePDB_DIR/final/template_pdbs/$templateID.rec_lig.$template_pdbIDchnIDs.pdb";


        system("$command") == 0 or die("FAILED: $command:$!");
    }

    printf
"Finish copy final superimposed models and template pdb files to the folder called final.\n";

}


#------------------------------------
# 1. -- cluster superimposed models

$skip = 0;

if ( $skip ne 1 ) {

    printf "\n\n======================================================";
    printf "\nCluster superimposed models ........................\n";
    printf "======================================================\n\n";

    my $command =
      "bash clusterStructures_oneCase.sh $supUnbound2TemplatePDB_DIR/final";
    print "COMMAND: $command\n";
    system($command) == 0 or die("FAiled: $command:$!");

    printf "\nFinish cluster superimposed models.\n\n";
}



#------------------------------------
# 2. -- calculate the distance restraints from each cluster of templates and combine them

$skip = 0;

if ( $skip ne 1 ) {

    printf "\n\n======================================================";
    printf
"\nCalculate CA-CA restraints for each cluster of templates ........................\n";
    printf "======================================================\n\n";

    my $atomDistThr = $caDistThr + 3
      ; #20 #-- atom distance threshold for calculating the initial contact file, which is later used to extract CA-CA distances

    foreach my $clusterDIR (`ls $supUnbound2TemplatePDB_DIR/final/cluster* -d`)
    {
        $clusterDIR =~ s/[\r\n]//mg;

        my $clusterNo = basename ($clusterDIR);

        printf "\n----------------------------------------------------------\n";
        printf "\n** Calculate CA-CA file for cluster $clusterNo **\n";
        printf "\n----------------------------------------------------------\n";

        my $command =
"bash distanceRestraints_oneCluster.sh $clusterDIR $atomDistThr $caDistThr";
        print "COMMAND: $command\n";
        system($command) == 0 or die("FAILED $command:$!");

    }

    printf
      "\nFinish calculating the contacts for each suerpoimposed models.\n\n";
}



#------------------------------------
# 3. -- move calculated predicted Ca-Ca distances for a specific folder
$skip = 0;

if ( $skip ne 1 ) {

    printf "\n\n\n======================================================";
    printf
"\nCopy predicted Ca-Ca distances to an output folder (where users can download) ..................\n";
    printf "======================================================\n\n";

    my $dirname   = dirname ($supUnbound2TemplatePDB_DIR);
    my $outputDIR = "$dirname/Ca_Ca_distances";

    printf "\n\nmove predicted Ca-Ca distances to $outputDIR\n\n";

    if ( !-d $outputDIR ) {

        `mkdir -p $outputDIR`;
    }

    my $num_good_cluster = 0;

    foreach my $clusterDIR (`ls $supUnbound2TemplatePDB_DIR/final/cluster* -d`)
    {
        $clusterDIR =~ s/[\n\r]//mg;

        my $clusterID = basename ($clusterDIR);
        print "clusterID: $clusterID\n";

        if ( -e "$clusterDIR/final_caca.txt" ) {
            copy( "$clusterDIR/final_caca.txt",
                "$outputDIR/$clusterID\_Ca_Ca_distance.txt" );
            copy(
                "$clusterDIR/restraints_4haddock.CaCa.txt",
                "$outputDIR/$clusterID\_restraints.tbl"
            );
            $num_good_cluster++;
        }
        else {
            printf "Cluster $clusterID does not have CA-CA file generated.\n";
        }
    }

    if ( $num_good_cluster > 0 ) {

        print
          "\nRestraint files from $num_good_cluster cluster(s) generated.\n";
        print "$outputDIR generated\n";
    }
    else {

        print "\nWARNING: no CA-CA file(s) generated for this query pair !!\n";
    }

}


#------------------------------------
# If user uploaded pdb files, then skip step 4 and 5.

my $skip_4;
my $skip_5;

my $jobDIR=dirname($PSHomPPI_result_onePair_DIR);
my $uploadedPDBdir="$jobDIR/pdb";

if (!-d $uploadedPDBdir){
    # user did not upload query pdb files
    # ps-homppi used tempalte1 as query pdb files in generating supimposed models
    # Do NOT skip step 4 and 5
    $skip_4 = 0;
    $skip_5 = 0;
}
else{
    # user uploaded query pdb files
    # Skip step 4 and 5
    $skip_4=1;
    $skip_5=1;
}



#------------------------------------
# 4. -- write the ID map file for query and superimposed models

my $qryPair = basename($PSHomPPI_result_onePair_DIR);    #'C:D'
my ($qry1, $qry2) = split(/:/, $qryPair);

if($skip_4 ==0){
    printf "\n\n\n======================================================\n";
    printf
    "Write ID maps for qry and superimposed models ...........................\n";
    printf "======================================================\n\n";

    my $IDmapFL =
    "$PSHomPPI_result_onePair_DIR/superimposed_models/IDmap_qry_supModels.txt";
    `echo "# generated by $0" > $IDmapFL`;
    `echo -e "QRY\tsuperimposed_models" >> $IDmapFL`;
    `echo -e "$qry1\tA" >> $IDmapFL`;
    `echo -e "$qry2\tB" >> $IDmapFL`;
    print "$IDmapFL generated\n";

}



#------------------------------------
# 5. -- add Qry chain IDs and atomResNums to predicted Ca-Ca distance file

if ($skip_5 ==0){
    printf "\n\n\n======================================================\n";
    printf "Add qry chain IDs and atomResNums to predicted Ca-Ca distance file .....................\n";
    printf "======================================================\n\n";

    #--step 1: generate residue number map between query and template 1.
    # output: qry1_4lqwA.resiNumMap, qry2_4lqwB.resiNumMap

    my $TemplatesUsed_DIR = "$PSHomPPI_result_onePair_DIR/TemplatesUsed";
    my $homcomplexes_used_in_predictionFL = "$TemplatesUsed_DIR/templates_used_in_prediction.stat";
    my $template1_ID = &getTemplate1_ID ($homcomplexes_used_in_predictionFL); # $template1_ID ->{qry1} = '1ahjA,23'

    my $hhpredFL_qry1 = "$PSHomPPI_result_onePair_DIR/$qry1/$qry1.hhpred";
    my $hhpredFL_qry2 = "$PSHomPPI_result_onePair_DIR/$qry2/$qry2.hhpred";

    my $outputDIR= "$TemplatesUsed_DIR/qry_template1_resNumMap";
    mkdir ($outputDIR) if ( ! -d $outputDIR );

    my $template1_rec = $template1_ID->{$qry1}; #--- 4lqwA,1
    my $template1_lig = $template1_ID->{$qry2}; #--- 4lqwC,26
    $template1_rec =~ s/,\d+//; #--- 4lqwA
    $template1_lig =~ s/,\d+//; #--- 4lqwC

    my $template1_qry1_resNumMapFL = "$outputDIR/$qry1\_$template1_rec.resiNumMap";
    my $template1_qry2_resNumMapFL = "$outputDIR/$qry2\_$template1_lig.resiNumMap";

    &genAtomResNumMapFL_new ($template1_ID->{$qry1}, $hhpredFL_qry1, $template1_qry1_resNumMapFL);
    &genAtomResNumMapFL_new ($template1_ID->{$qry2}, $hhpredFL_qry2, $template1_qry2_resNumMapFL);

    #--step 2: add qry residue number and qry chnID to the predicted Ca-Ca distance file
    my $command = " perl mapDistFL2Qry.pl $PSHomPPI_result_onePair_DIR";
    print "COMMAND: $command\n";
    system($command) ==0 or die ("FAILED: $command:$!");

    print "Qry chain IDs and atomResNums added to the CA-CA distance file.\n";

}

#------------------------------------
# 6. --  write haddock restraint files

my @cacaFLs = `ls $PSHomPPI_result_onePair_DIR/superimposed_models/Ca_Ca_distances/cluster*_Ca_Ca_distance.txt`;
foreach my $cacaFL (@cacaFLs){
    # $cacaFL = cluster1_Ca_Ca_distance.txt
    $cacaFL =~s/[\n\r]//gm;

    #-- set the name of output file
    my $filename=basename($cacaFL);
    my @tmp = split(/_/, $filename);
    my $clusterID = shift @tmp; # $clusterID = 'cluster1'
    my $dirname=dirname($cacaFL);
    my $haddockRestFL = "$dirname/$clusterID\_restraints.tbl";

    #-- write $haddockRestFL

    my $command = "perl genHaddockRestrFL.pl $cacaFL $haddockRestFL $caDistThr ";
    system($command) ==0 or die ("FAILED: $command:$!");
}


