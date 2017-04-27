#!/usr/bin/env perl -w
#
# Li Xue
# Mar 1st, 2013
#
# Superimpose unbound receptor and ligand to a template complex,
# write the final PDB.
#
#  ******
#  NOTE: Given $PSHomPPIoutputDIR='.../C_D', the chain IDs for rec (A) and lig (B) in the final superimposed models corresponds to C and D, respectively.
#  ******
#
# This script prepares files for clustering of templates.
#
# step1: prepare template PDB files using TemplateUsedByPSHomPPI()
# step2: generate superimposed pdb file of the first template structure over other templates
#
#
# TMalign is used for the superimposition.
# 	>TMalign prot1_pdbFL prot2_pdbFL
# 1. The coordinates of prot2 is kept unchanged.
# 2. If prot1 contains multiple chains, TMalign automatically combine them into one chain and use chain ID A in the output.
# Same for prot2, chain ID B is used in the output.
# 3. TMalign output file has nothing after the cordinates. It lost the information about ocupancy,  B-factor, etc.
# 4. TMalign output file will lose the HETATM section if the original pdb file has this section.
# 5. TMalign does not change the atomResNum.
# 6. If a non-amino acid is written as ATOM instead of HETATM, TMalign will write this chemical as X in the alignment file.
# 7. If an atom has alternative location (pos. 17), TMalign will only use the first location and discard other locations.

# perl genSuperimposePDBFL.pl  $PSHomPPIoutputDIR
# perl genSuperimposePDBFL.pl /home/lixue/DockRank_HADDOCK/callHADDOCK_CaCa/data/BM4_dimer_Ca_Ca/predicted_int_PSHomPPI
# perl genSuperimposePDBFL.pl /data/web_servers/lxue/PSHOMPPIv1.3/uploadData/Tue_Jan_27_09_37_36_2015.15.T80/A_B
#
use strict;
use File::Path;
use File::Basename;

use globalVariables;
our $perl5LibDIR;
use lib $perl5LibDIR;
use Cwd;


#--- global variables
our $unboundPDB_DIR ; # '/data/lixue/DBs/ppidev/ppiInfoNew/asymmetric/pdb';    #--to get pdb files for templates
our $pdb_chainPY    ; # '/home/lixue/tools/pdb-tools/pdb_chain.py'; #set chain ID for a pdb file
our $TMalign        ; # '/home/lixue/tools/TMtool/TMalign';
our $pdb2atomresnumPL; # '/home/lixue/tools/PDB2AtomResNum.pl';

#--
if ( ! -d $perl5LibDIR ){
    die("$perl5LibDIR does not exist:$!");
}

if ( !-d $unboundPDB_DIR ) {
    die("$unboundPDB_DIR does not exist:$!");
}

#--input variables
my $PSHomPPIoutputDIR =
  shift @ARGV;    #../../uploadData/Tue_Jan_27_09_37_36_2015.15.T80/C:D
#$PSHomPPIoutputDIR=File::Spec::Unix->abs2rel($PSHomPPIoutputDIR);

my $jobDIR = dirname($PSHomPPIoutputDIR); # ../../uploadData/Tue_Jan_27_09_37_36_2015.15.T80
my $IDpair = basename ($PSHomPPIoutputDIR); #C:D


#--
my $superimDIR = "$PSHomPPIoutputDIR/superimposed_models";
my $outputDIR = "$superimDIR/supUnbound2TemplatePDB";
rmtree($superimDIR) if (-d $superimDIR );
mkdir($superimDIR);
mkdir $outputDIR;

#--prepare template pdb files

&TemplateUsedByPSHomPPI_oneCase_new($PSHomPPIoutputDIR);

#--superimpose unbound pdb to template pdb
#
if ( !-d "$PSHomPPIoutputDIR/TemplatesUsed" ) {
    die(
"$PSHomPPIoutputDIR/TemplatesUsed does not exist. Need to run \&TemplateUsedByPSHomPPI:$!"
    );
}
print
"\nSuperimpose template 1 (as if it is the query unbound pdbs) to other templates  ...\n\n";

print "\n\n\t============================================================\n\n";
print "\t\tGenerate superimposed model for this case   ... ";
print "\n\n\t============================================================\n\n";

my $caseOutputDIR  = "$outputDIR";
my $templatePDBDIR = "$PSHomPPIoutputDIR/TemplatesUsed/pdb";

#--
rmtree($caseOutputDIR) if ( -d $caseOutputDIR );

#--
if ( !-d $templatePDBDIR ) {
    print
"\n**Warning: this case does not have template PDB DIR: $templatePDBDIR. Maybe this case is NOT a dimer OR do not templates from PS-HomPPI. We skip such cases so far!!!\n\n";
    exit 0;

}

#--

opendir( DIR, $templatePDBDIR )
  || die("Cannot open folder $templatePDBDIR:$!");
my @templateFLnames = grep { /^template.+.pdb/ } readdir(DIR);
closedir(DIR);

if ( scalar @templateFLnames == 0 ) {

    print "**Warning: No template files under $templatePDBDIR:$!";
    next;
}

#--
mkdir $caseOutputDIR if ( !-d $caseOutputDIR );

#--
my $templateFLname
  ;    # $templateFLname->{'template1'}->{'rec'}= 'template1.rec.2qyi.pdb';
foreach my $FLname (@templateFLnames) {
    my ( $templateID, $recOrLig, $pdbID ) =
      split( /\./, basename( $FLname, '.pdb' ) );

    $templateFLname->{$templateID}->{$recOrLig} = $FLname
      ;    # $templateFLname->{'template1'}->{'rec'}= 'template1.rec.2qyi.pdb';
}

#--------------------------------------------------
#-- choose the query PDB files: either use User-uploaded ones or use template 1

my ( $unbound_r_PDBFL, $unbound_l_PDBFL );
if (!-d "$jobDIR/pdb"){

    #-- user did NOT upload query PDB file
    #-- treat Template 1 as the structure of the query and superimpose it to other templates

    my $rec_pdbFL_template1 =
    "$templatePDBDIR/$templateFLname->{'template1'}->{'rec'}"
    ;    #check. rec is aligned only partly with the template
    my $lig_pdbFL_template1 =
    "$templatePDBDIR/$templateFLname->{'template1'}->{'lig'}";

     ( $unbound_r_PDBFL, $unbound_l_PDBFL ) =
    &preparePDBFLs4qry( $rec_pdbFL_template1, $lig_pdbFL_template1,
        $caseOutputDIR );

    #--$unbound_r_PDBFL and $unbound_l_PDBFL is the pdb file that will be used as query PDB, which are used later to generate superimposed models
}
else{
    #-- user provided query PDB files
    my ($rec,$lig) = split(/:/, $IDpair);
    my $rec_pdbFL_qry = "$jobDIR/pdb/$rec.pdb";
    my $lig_pdbFL_qry = "$jobDIR/pdb/$lig.pdb";

     ( $unbound_r_PDBFL, $unbound_l_PDBFL ) =
    &preparePDBFLs4qry( $rec_pdbFL_qry, $lig_pdbFL_qry,
        $caseOutputDIR );

    #--$unbound_r_PDBFL and $unbound_l_PDBFL is the pdb file that will be used as query PDB, which are used later to generate superimposed models
}


#-------------------------------------------------

foreach my $templateID ( sort keys %$templateFLname ) {

    #$templateID = 'template2'
    print
"\n\n\t--------------------------------------------------------------\n\n";
    print "\t\tGenerate superimposed model based on $templateID ... ";
    print
"\n\n\t--------------------------------------------------------------\n\n";

    my $caseTemplateOutputDIR = "$caseOutputDIR/$templateID";
    my $template_r_PDBFL =
      "$templatePDBDIR/$templateFLname->{$templateID}->{'rec'}"
      ;    #check. rec is aligned only partly with the template
    my $template_l_PDBFL =
      "$templatePDBDIR/$templateFLname->{$templateID}->{'lig'}";

    if ( !-e $template_r_PDBFL || !-e $template_l_PDBFL ) {
        die("$template_r_PDBFL and/or $template_l_PDBFL do not exist:$!");
    }
    print "Template rec PDB: $template_r_PDBFL\n";
    print "Template lig PDB: $template_l_PDBFL\n";

    &genSuperimposePDBFL_oneTemplate(
        $caseTemplateOutputDIR, $template_r_PDBFL, $template_l_PDBFL,
        $unbound_r_PDBFL,       $unbound_l_PDBFL
    );
}

print "\n\n--------\n$outputDIR generated.\n";

#------------------

sub genSuperimposePDBFL_oneTemplate {

#Li Xue
#Jan 12, 2013
#
#Superimpose unbound receptor and ligand to a template complex,
#write the final PDB, for which chain A is rec, and B is lig.
#
#TMalign is used for the superimposition.
#	>TMalign prot1_pdbFL prot2_pdbFL
#1. The coordinates of prot2 is kept unchanged.
#2. If prot1 contains multiple chains, TMalign automatically combine them into one chain and use chain ID A in the output.
#Same for prot2, chain ID B is used in the output.
#3. TMalign output file has nothing after the cordinates. It lost the information about ocupancy,  B-factor, etc.
#4. TMalign output file will lose the HETATM section if the original pdb file has this section.
#5. TMalign does not change the atomResNum.
#6. If a non-amino acid is written as ATOM instead of HETATM, TMalign will write this chemical as X in the alignment file.
#7. If an atom has alternative location (pos. 17), TMalign will only use the first location and discard other locations.

#Usage:
#perl genSuperimposePDBFL.pl caseID template_rec_pdb template_lig_pdb outputDIR
#perl genSuperimposePDBFL.pl 1AHW ../data/Cluspro2_BM3_decoy/predicted_int_ClusPro2_EngyFun000_5angstrom_SeqIdenThr90_strict_newPSHomPPI/supUnbound2TemplatePDB
#perl genSuperimposePDBFL.pl 1CGI ../data/Cluspro2_BM3_decoy/predicted_int_ClusPro2_EngyFun000_5angstrom_SeqIdenThr90_strict_newPSHomPPI/supUnbound2TemplatePDB

    use strict;
    use File::Basename;
    use File::Copy;
    use File::Path;

    #--input
    if ( scalar @_ ne 5 ) {
        die("Number of input ARGV is wrong:$!");

    }

    my $outputDIR = shift @_;    #put all the output of the superimposition here
    my $template_r_PDBFL  = shift @_;    #dir/template1.rec.1m9xA.pdb
    my $template_l_PDBFL  = shift @_;    #dir/template1.lig.1m9xA.pdb
    my $PDBFL_unbound_rec = shift @_;
    my $PDBFL_unbound_lig = shift @_;

    #--
    my ($template_pdbID_rec) = $template_r_PDBFL =~ /.(\w+).pdb$/;
    my ($template_pdbID_lig) = $template_l_PDBFL =~ /.(\w+).pdb$/;

    #--
    my $TMalignOutputDIR = "$outputDIR/original_TMalign_output";

    #--output files
    my $basename1 = basename( $PDBFL_unbound_rec, '.pdb' );
    my $supFL_r_name =
      "$basename1.sup";    #the basic output file name of TMalign for receptor
    my $basename2 = basename( $PDBFL_unbound_lig, '.pdb' );
    my $supFL_l_name =
      "$basename2.sup";    #the basic output file name of TMalign for ligand

    my $supFL_RecTem_pdb =
      "$TMalignOutputDIR/$supFL_r_name\_all_atm"
      ;                    #the pdb file of supimposed rec and template
    my $supFL_LigTem_pdb =
      "$TMalignOutputDIR/$supFL_l_name\_all_atm"
      ;                    #the pdb file of supimposed rec and template

    my $supFL_Rec_pdb =
      "$outputDIR/$supFL_r_name.pdb"
      ;                    #the pdb file of superimposed rec and template
    my $supFL_Lig_pdb =
      "$outputDIR/$supFL_l_name.pdb"
      ;                    #the pdb file of superimposed rec and template

    #--
    rmtree($outputDIR)       if ( -d $outputDIR );
    mkdir($outputDIR)        if ( !-d $outputDIR );
    mkdir($TMalignOutputDIR) if ( !-d $TMalignOutputDIR );

    #-- clean template pdb file

    my $basename             = basename($template_r_PDBFL);
    my $template_r_PDBFL_new = "$outputDIR/$basename";
    $basename = basename($template_l_PDBFL);
    my $template_l_PDBFL_new = "$outputDIR/$basename";

    &cleanPDBFL( $template_r_PDBFL, $template_r_PDBFL_new );
    &cleanPDBFL( $template_l_PDBFL, $template_l_PDBFL_new );

    #--align rec to template
    print "\n\n*** Superimpose rec and lig to the template... ***\n\n";
    print "\n\nAlign rec to template ...\n\n";
    &callTMalign(
        $PDBFL_unbound_rec, $template_r_PDBFL_new,
        $TMalignOutputDIR,  $supFL_r_name
    );

    print "\n\nAlign lig to template ...\n\n";
    &callTMalign(
        $PDBFL_unbound_lig, $template_l_PDBFL_new,
        $TMalignOutputDIR,  $supFL_l_name
    );
    print "\n\ncallTMalign done.\n\n";

    #--get the 3D coordinates of aligned rec and output PDB file $supFL_Rec_pdb
    my $newChnID_rec = 'A';
    &parseSupFL( $supFL_RecTem_pdb, $supFL_Rec_pdb, $newChnID_rec,
        $PDBFL_unbound_rec, "$PDBFL_unbound_rec.HETATM" );

  #../data/Cluspro2_BM3_decoy/benchmark3/originalPdb/1AHW_r_u.pdb has two chains

#--get the 3D coordinates of aligned rec and output PDB file $supFL_Rec_pdb, and set the chain ID to B
    my $newChnID_lig = 'B';
    &parseSupFL( $supFL_LigTem_pdb, $supFL_Lig_pdb, $newChnID_lig,
        $PDBFL_unbound_lig, "$PDBFL_unbound_lig.HETATM" );

    #-- combine the superimposed rec and lig pdb files into one pdb file.

    print "\nCombine the coordinates of lig and rec into one PDB file...\n";
    my $final_supFL =
      "$outputDIR/finalSup.$template_pdbID_rec\_$template_pdbID_lig.pdb";
    &combineTwoPDBFL( $supFL_Rec_pdb, $supFL_Lig_pdb, $final_supFL );

}

#------------------------
#
sub renumberPDBFL {

 #-- renumber the at)mResNum so that it starts from 1
 #-- so that it is easier to know map the alignment file with the pdb file later
 #-- pdb_reres does not work with long dir
    my $pdbFL      = shift @_;
    my $currentDIR = getcwd;
    my $dataDIR    = dirname($pdbFL);
    my $filename   = basename($pdbFL);
    chdir $dataDIR;
    system("pdb_reres $filename -1 > $filename.tmp; mv $filename.tmp $filename")
      == 0
      or die("Cannot renumber the atomResNum in $pdbFL:$!");
    chdir $currentDIR;
    print "$pdbFL is renumbered to starting from 1.\n";

}

sub changChnID {
    use File::Copy;

    my $pdbFL    = shift @_;
    my $newChnID = shift @_;
    my $r        = rand(10);
    my $tempFL   = "$r.pdb";

    unlink $tempFL if ( -e $tempFL );
    open( OUTPUT, ">>$tempFL" ) || die("Cannot open $tempFL:$!");
    open( INPUT,  "<$pdbFL" )   || die("Cannot open $pdbFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^(ATOM|HETATM)/) {

            #http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
            substr( $_, 21, 1 ) = $newChnID;
            print OUTPUT "$_\n";

        }
        else {

            print OUTPUT "$_\n";
        }

    }
    close INPUT;
    close OUTPUT;

    unlink $pdbFL;
    copy( $tempFL, $pdbFL ) || die("Cannot rename $tempFL to $pdbFL:$!");
    unlink $tempFL;

}


sub callTMalign {
    our $TMalign;
    our $pdb2atomresnumPL;

    #--
    my $mobilePDBFL   = shift @_;
    my $templatePDBFL = shift @_;
    my $outputDIR     = shift @_;
    my $outputName    = shift @_;

    if ( !-e $mobilePDBFL || !-e $templatePDBFL ) {
        die("$mobilePDBFL or $templatePDBFL does not exist:$!");
    }
    my $command =
"$TMalign $mobilePDBFL $templatePDBFL -o $outputDIR/$outputName  > $outputDIR/$outputName.align";
    print "\n$command\n";
    system($command) == 0
      or die("TMalign does not work: $command : $!")
      ;    #the coordinates of the 2nd proteins will not change.

    #---------------------------
    #-- a side-step:
    #-- format the alignment file generated by TMalign so that it can be used later for predicted CA-CA distance calculation
    #-- When predicting CA-CA interface distances later, only aligned residues (corresponding the alignment symbol *** (aligned near) and ** (aligned far) ) are used.

    my $aligned_resiNumFL = "$outputDIR/$outputName.aligned_resiNum";
    &formatAlignFL( "$outputDIR/$outputName.align", $aligned_resiNumFL );

    #--replace the sequence resinum in $aligned_resiNumFL with atom resinum
    my $basename_Q       = basename( $mobilePDBFL, ".pdb" );
    my $dirname_Q        = dirname($mobilePDBFL);
    my $mapFL_atom2seq_Q = "$dirname_Q/$basename_Q.atomResNum";

    my $basename_T       = basename( $templatePDBFL, ".pdb" );
    my $dirname_T        = dirname($templatePDBFL);
    my $mapFL_atom2seq_T = "$dirname_T/$basename_T.atomResNum";

    system("$pdb2atomresnumPL $mobilePDBFL > $mapFL_atom2seq_Q") == 0
      or die("Cannot extract atomRestNum from $mobilePDBFL:$!");
    print "$mapFL_atom2seq_Q generated\n\n";

    system("$pdb2atomresnumPL $templatePDBFL > $mapFL_atom2seq_T") == 0
      or die("Cannot extract atomRestNum from $templatePDBFL:$!");
    print "$mapFL_atom2seq_T generated\n\n";

    #-replace seqResNum with atomResNum
    system(
"perl mapAtomResNum_4alnSymbolFL.pl $aligned_resiNumFL $mapFL_atom2seq_Q  $mapFL_atom2seq_T > $aligned_resiNumFL.tmp && mv $aligned_resiNumFL.tmp $aligned_resiNumFL "
      ) == 0
      or die(
"mapAtomResNum_4alnSymbolFL.pl failed to replace seqResNum in $aligned_resiNumFL with atomResNum in $mapFL_atom2seq_Q and $mapFL_atom2seq_T:$!"
      );
    print "seqResNum in $aligned_resiNumFL is replaced with atomResNum.\n";


    #- add chain IDs to $aligned_resiNumFL
    my @chnIDs_mobilePDB = @{&getChnID($mobilePDBFL)};
    my @chnIDs_templatePDB = @{&getChnID($templatePDBFL)};
    if (scalar @chnIDs_mobilePDB ne 1){
        die("The query pdb file $mobilePDBFL should have one and only one chain ID:$!");
    }

    &addChnID($aligned_resiNumFL,$chnIDs_mobilePDB[0], $chnIDs_templatePDB[0]);


    #-- side step done. ----------------



    #---------------------------
    #-- copy final output file to ..
    #
    my $upperDIR = dirname($outputDIR);
    system("mv $aligned_resiNumFL $upperDIR") == 0
      or die("Cannot move $aligned_resiNumFL to $upperDIR:$!");

}

#sub callTMalign {
#    our $TMalign;
#    our $pdb2atomresnumPL;
#
#    #--
#    my $mobilePDBFL   = shift @_;
#    my $templatePDBFL = shift @_;
#    my $outputDIR     = shift @_;
#    my $outputName    = shift @_;
#
#    if ( !-e $mobilePDBFL || !-e $templatePDBFL ) {
#        die("$mobilePDBFL or $templatePDBFL does not exist:$!");
#    }
#    my $command =
#"$TMalign $mobilePDBFL $templatePDBFL -o $outputDIR/$outputName  > $outputDIR/$outputName.align";
#    print "\nTMALIGN COMMAND: $command\n\n";
#    system($command) == 0
#      or die("TMalign does not work: $command : $!")
#      ;    #the coordinates of the 2nd proteins will not change.
#
##-- a side-step:
##-- format the alignment file generated by TMalign so that it can be used later for predicted CA-CA distance calculation
##-- When predicting CA-CA interface distances later, only well-aligned residues (corresponding the alignment symbol *** ) are used.
#    my $aligned_resiNumFL = "$outputDIR/$outputName.aligned_resiNum";
#    &formatAlignFL( "$outputDIR/$outputName.align", $aligned_resiNumFL );
#
#    #--replace the sequence resinum in $aligned_resiNumFL with atom resinum
#    print
#"\n\n-- replace the seqResNum in the alignment file generated by TMalign with atomResNum ...\n\n";
#
#    my $basename       = basename( $mobilePDBFL, ".pdb" );
#    my $dirname        = dirname($mobilePDBFL);
#    my $mapFL_atom2seq = "$dirname/$basename.atomResNum";
#
#    $command = "perl $pdb2atomresnumPL $mobilePDBFL > $mapFL_atom2seq";
#    print "COMMAND: $command\n\n";
#    system("$command") == 0
#      or die("Cannot extract atomRestNum from $mobilePDBFL:$!");
#
#    $command =
#"perl mapAtomResNum_4alnSymbolFL.pl $aligned_resiNumFL $mapFL_atom2seq > $aligned_resiNumFL.tmp && mv $aligned_resiNumFL.tmp $aligned_resiNumFL";
#    print "COMMAND: $command\n\n";
#    system("$command") == 0
#      or die(
#"\nERROR:mapAtomResNum_4alnSymbolFL.pl failed to replace seqResNum in $aligned_resiNumFL with atomResNum in $mapFL_atom2seq:$!");
#    print "seqResNum in $aligned_resiNumFL is replaced with atomResNum.\n\n";
#
#    #-- copy final output file to ..
#    my $upperDIR = dirname($outputDIR);
#    system("mv $aligned_resiNumFL $upperDIR") == 0
#      or die("Cannot move $aligned_resiNumFL to $upperDIR:$!");
#
#}
sub addChnID{
    # add chain ID to the input file
    # Input file format:
    ##Qry_aa Qry_seqResNum Template_aa Template_seqResNum alignment_symbol
    #Q   1   1   M   1   1   ***
    #T   2   2   V   2   2   ***

    my $aligned_resiNumFL =shift @_;
    my $qry_chnID = shift @_;
    my $template_chnID = shift @_;
    my $tmpFL = "$aligned_resiNumFL.tmp";
    unlink $tmpFL if (-e $tmpFL);
    open (OUTPUT, ">>$tmpFL") or die ("Cannot open $tmpFL:$!");

    open(INPUT, "<$aligned_resiNumFL") or die ("Cannot open $aligned_resiNumFL:$!");
    while(<INPUT>){
        s/[\n\r]//mg;

        if (/^#{0,}Qry/i){
            my @header_tmp = qw (Qry_chnID Qry_aa Qry_seqResNum Qry_atomResNum Template_chnID Template_aa Template_seqResNum Template_atomResNum Symbol);
            my $header = join("\t", @header_tmp);
            print OUTPUT "$header\n";
            next;
        }
        if (/^#/){
            print OUTPUT "$_\n";
            next;
        }
        if (/^\w{1}\s+/){
            my @a=split(/\s+/, $_);
            splice (@a, 0, 0, $qry_chnID);
            splice (@a, 4, 0, $template_chnID);
#            my  ($aa_Q, $seqResNum_Q,$atomResNum_Q, $aa_T, $seqResNum_T , $atomResNum_T, $modes) =split(/\s+/,$_);
            my $newLine = join("\t", @a);
            print OUTPUT "$newLine\n";

        }
    }
    close INPUT;
    close OUTPUT;

    move($tmpFL, $aligned_resiNumFL) or die ("$!");

    print "chnIDs added to $aligned_resiNumFL\n";
}

sub getChnID{
    my $pdbFL = shift @_;
    my @chnIDs;

    open (INPUT, "<$pdbFL") or die ("Cannot open $pdbFL:$!");
    while(<INPUT>){
        s/[\n\r]//mg;
        if (/^(ATOM|HETATM)/) {
           my $chnID = substr($_,21,1);
           push @chnIDs, $chnID;
        }
    }
    close INPUT;

    if (!@chnIDs){
        die("No chain IDs read from $pdbFL:$!");
    }

    my @final_chnIDs = @{&Unique(\@chnIDs)};

    return \@final_chnIDs;

}

sub Unique {
    my @a = @{shift @_} ;
    my %seen;
    @seen{@a}=1 x scalar @a;
    my @a_new = keys %seen;
    return \@a_new;
}


sub addTER {

    my $pdbFL_ori = shift @_;
    my $pdbFL_new = shift @_;

    #    print "\nAdding TER to the end of chains in $pdbFL_new ...\n";

    unlink $pdbFL_new if ( -e $pdbFL_new );
    open( OUTPUT, ">>$pdbFL_new" ) || die("Cannot open $pdbFL_new:$!");

    my $chainID_prev;

    open( INPUT, "<$pdbFL_ori" ) || die("Cannot open $pdbFL_ori:$!");
    while (<INPUT>) {

        s/[\n\r]//gm;

        if (/^TER/) {

#the original file has TER already: /data/benchmark/docking-benchmark4/runs-cmrestraints/1Z0K/ana_scripts/protein1.pdb
            next;
        }

        if (/^(ATOM|HETATM)/) {

            #http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
            my $chnID = substr( $_, 21, 1 );

            if ( !defined $chainID_prev ) {

                #first line of the input pdb file
                $chainID_prev = $chnID;
                print OUTPUT "$_\n";
                next;
            }
            else {

                if ( $chnID eq $chainID_prev ) {
                    print OUTPUT "$_\n";
                    $chainID_prev = $chnID;
                    next;
                }
                else {
                    print "previous chain ID: $chainID_prev\n";
                    print "new chain ID: $chnID\n";

                    print OUTPUT "TER\n";
                    print OUTPUT "$_\n";

                    $chainID_prev = $chnID;

                }

            }

        }
        else {

            #not ATOM line
            print OUTPUT "$_\n";
        }

    }
    close INPUT;
    close OUTPUT;

    #    print "$pdbFL_new generated.\n\n";
    print "Add TER done.\n";

}

sub combineTwoPDBFL {

    my $pdbFL1   = shift @_;
    my $pdbFL2   = shift @_;
    my $outputFL = shift @_;

    unlink $outputFL if ( -e $outputFL );

    system("cat $pdbFL1> $outputFL");
    system("echo \"TER\">>$outputFL");
    system("cat $pdbFL2>>$outputFL");
    system("echo \"TER\">>$outputFL");

    print "$outputFL generated. Chn A is rec, and Chn B is lig.\n";

}

sub parseSupFL {

    my $TMalignOutputFL = shift @_; #input
    my $PDBFL_supModel  = shift @_; #output
    my $supModel_chnID  = shift @_; #set the chain ID for the superimposed model
    my $PDBFL_qry_unbound = shift @_;    #used to add tail
    my $PDBFL_qry_HETATM  = shift @_;    #used  to add HETATM to the sup model

    &getCoor4supModel( $TMalignOutputFL, $PDBFL_supModel );

    #--
    #both $supFL_Lig_pdb and $supFL_Rec_pdb have chain ID A .
    #change $supFL_Lig_pdb to chain B
    &changChnID( $PDBFL_supModel, $supModel_chnID );

#--
# Add back the occupancy and B-factor from the original pdb file to the pdb file generated by TMalign
# NOTE: remove altLoc not solved yet. So this step cannot be on.
#    &addTail( $PDBFL_supModel, $PDBFL_qry_unbound );

#--
# Add back the HETATM section from the original pdb file to the pdb file generated by TMalign
    if ( -e $PDBFL_qry_HETATM ) {
        &addHETATM( $PDBFL_supModel, $PDBFL_qry_HETATM );
    }

#	#--
#	# Check the number of rows
#	if(!&equalRowNum($supFL_Lig_pdb,$PDBFL_unbound_lig_new)){
#			die("superimposed file do not have the same number of rows as the original pdb file. Check $supFL_Lig_pdb and $PDBFL_unbound_lig_new :$!");
#		})
#
#	if(!&equalRowNum($supFL_Rec_pdb,$PDBFL_unbound_rec_new)){
#			die("superimposed file do not have the same number of rows as the original pdb file. Check $supFL_Rec_pdb and $PDBFL_unbound_rec_new :$!");
#		}
#
#
#-- renumber the atomResNum so that it starts from 1
#    &renumberPDBFL($PDBFL_supModel);

}

sub getCoor4supModel {

#parse the output all_atm PDB file of TMalign. This function extract the first chain (chain ID A) and write it into a separate pdb file.

#	>TMalign prot1_pdbFL prot2_pdbFL
#The coordinates of prot2 is kept unchanged.
#If prot1 contains multiple chains, TMalign automatically combine them into one chain and use chain ID A in the output.
#Same for prot2, chain ID B is used in the output.

#Input file:
#
#	load inline
#	select *A
#	color blue
#	select *B
#	color red
#	select all
#	cartoon
#	exit
#	REMARK TM-align Version 20120707
#	REMARK Chain 1:../data/Cl  Size= 428
#	REMARK Chain 2:../data/Cl  Size= 304 (TM-score is normalized by  304, d0=  6.40)
#	REMARK Aligned length= 155, RMSD=  5.90, TM-score=0.32056, ID=0.065
#	ATOM      1  N   ASP A   1      24.492 -11.544  15.114
#	ATOM      2  CA  ASP A   1      24.632 -13.020  15.231
#	ATOM      3  C   ASP A   1      23.710 -13.671  14.207
#	ATOM      4  O   ASP A   1      24.167 -14.493  13.421
#	ATOM      5  CB  ASP A   1      24.325 -13.497  16.658
#	ATOM      6  CG  ASP A   1      24.472 -15.016  16.831
#	ATOM      7  OD1 ASP A   1      25.434 -15.607  16.286
#

    my $TMalignOutputFL = shift @_;
    my $outputFL        = shift @_;

    print
"\nExtracting the new coordinates of rec/lig PDB FL from TMalign output ....\n";

    unlink $outputFL if ( -e $outputFL );

    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
    open( INPUT, "<$TMalignOutputFL" )
      || die("Cannot open $TMalignOutputFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^TER/) {
            last;

        }
        if (/^(ATOM|HETATM)/) {

            print OUTPUT "$_\n";
        }

    }
    close INPUT;
    close OUTPUT;

    print "$outputFL generated.\n";

}

#---------------------

sub TemplateUsedByPSHomPPI_oneCase_old {

#Li Xue
#Jan 14th, 2013
#
#Input files: homologs_used_in_predictions.lst for each chain of a query complex
#Output format:
#	rec	rec	    lig
#	A	B		C
#	1awcA	1awcB	1awcE
#
#Usage:
#perl TemplateUsedByPSHomPPI.pl caseID PSHomPPIoutputDIR outputDIR_of_this_script

    use strict;
    use PSHomPPI_resiPairs;

    my $caseID = shift @_;    #'1CGI';    #
    my $PSHomPPIoutputDIR =
      shift @_
      ; #$PSHomPPIoutputDIR='../data/Cluspro2_BM3_decoy/predicted_int_ClusPro2_EngyFun000_5angstrom_SeqIdenThr90_strict_newPSHomPPI';
    my $outputDIR      = "$PSHomPPIoutputDIR/TemplatesUsed";
    my $outputDIR_case = "$outputDIR/$caseID";

    #--
    mkdir($outputDIR)      if ( !-d $outputDIR );
    mkdir($outputDIR_case) if ( !-d $outputDIR_case );

    #--
    my $case_predictedIntDIR = "$PSHomPPIoutputDIR/$caseID/predicted_int";

    my $recLigLst = "$case_predictedIntDIR/$caseID.lst";
    my ( $receptors, $ligands ) = &readRecLigLst($recLigLst);

    if ( scalar @$receptors > 1 || scalar @$ligands > 1 ) {

        print "\nWarning: This case is not dimmer. Skip.\n\n";
        return;
    }

    my %homologsUsedInPredictionFLs_rec;
    my %homologsUsedInPredictionFLs_lig;
    foreach my $rec (@$receptors) {

        foreach my $lig (@$ligands) {
            my $pair = "$rec:$lig";
            print "$pair\n";
            my ( $pairDIR, $chn1, $chn2 ) =
              &getPairDIR2( $pair, $case_predictedIntDIR );

            my $homologsUsedInPredictionFL_rec =
              "$pairDIR/$rec/homointerologs_used_in_prediction.txt";
            push @{ $homologsUsedInPredictionFLs_rec{$rec} },
              $homologsUsedInPredictionFL_rec;
            my $homologsUsedInPredictionFL_lig =
              "$pairDIR/$lig/homointerologs_used_in_prediction.txt";
            push @{ $homologsUsedInPredictionFLs_lig{$lig} },
              $homologsUsedInPredictionFL_lig;
        }

    }

    print
      "File for rec1: @{$homologsUsedInPredictionFLs_rec{$$receptors[0]}}\n";
    print "File for lig: @{$homologsUsedInPredictionFLs_lig{$$ligands[0]}}\n";

    #--

    foreach my $rec ( @{$receptors} ) {
        my $homologFL = "$outputDIR_case/homologsUsed_rec_Chn$rec.lst";
        &transformHomointerologFL2HomologFL(
            $homologsUsedInPredictionFLs_rec{$rec}, $homologFL )
          ;    #output: $homologFL

    }
    foreach my $lig ( @{$ligands} ) {
        my $homologFL = "$outputDIR_case/homologsUsed_lig_Chn$lig.lst";
        &transformHomointerologFL2HomologFL(
            $homologsUsedInPredictionFLs_lig{$lig}, $homologFL )
          ;    #output: $homologFL

    }

    #--
    my $homologFL_rec =
      "$outputDIR_case/homologsUsed_rec_Chn$receptors->[0].lst";
    my $homologFL_lig = "$outputDIR_case/homologsUsed_lig_Chn$ligands->[0].lst";

    if ( !-e $homologFL_rec || !-e $homologFL_lig ) {

        print
"Warning: No homologs found for rec or lig. homologFL_rec ($homologFL_rec) or homologFL_lig ($homologFL_lig) does not exist! No template file can be written!\n";

        return;

    }
    my @pdbIDs = @{
        &compareHomlogLsts_new( $homologFL_rec, $homologFL_lig,
            $outputDIR_case )
    };

    #--download pdb file of the templates
    my $pdbDIR = "$outputDIR_case/pdb";
    mkdir $pdbDIR if ( !-d $pdbDIR );

    foreach my $templatePDBID (@pdbIDs) {
        $templatePDBID = lc($templatePDBID);

        if ( !-e "$pdbDIR/$templatePDBID.pdb" ) {

            system("cp $unboundPDB_DIR/pdb$templatePDBID.ent.gz $pdbDIR") == 0
              or die("Cannot copy pdb file from ProtInDB to $pdbDIR:$!");
            system("gunzip $pdbDIR/pdb$templatePDBID.ent.gz") == 0
              or die("Cannot unzip pdb$templatePDBID.ent.gz:$!");
            move( "$pdbDIR/pdb$templatePDBID.ent",
                "$pdbDIR/$templatePDBID.pdb" );

            #            system(
            #"wget  http://www.rcsb.org/pdb/files/$templatePDBID.pdb -P $pdbDIR"
            #            );

        }
    }

    #--extract chains used as the templates from their pdb file
    my $templateFL = "$outputDIR_case/templates.lst";

    if ( !-e $templateFL ) {

        print
"No hom-complexexes found by PS-HomPPI. No template file can be written!\n";
        return;
    }
    &prepareRecLigTemplatePDB( $templateFL, $pdbDIR );
}

##----------------------------
#sub readTemplateFL{
##generated by genSuperimposePDBFL.pl
##pdbID  rec:lig
##3rdz    B:D
##3rdz    A:C
#
#my $templateFL = shift @_;
#my @templates;
#
#open(INPUT, "<$templateFL")||die("Cannot open $templateFL:$!");
#while(<INPUT>){
#
#	s/[\n\r]//gm;
#
#	if(/^\w{4}[\s\t]+\w+:\w+/){
#
#		push @templates, $_;
#	}
#
#}
#close INPUT;
#return \@templates;
#
#}
#
#
#------------------------
#sub extractPartPDBFL {
#
#	#extract the part of $chains out of $pdbFL_ori
#
#
#
#	my $pdbFL_ori = shift @_;
#	my $chains    = shift @_;    #'AC'
#	my $pdbFL_new = shift @_;
#
#
#	print "Extract the part of chain $chains out of $pdbFL_ori ...\n";
#
#	if($chains !~/^[a-zA-Z\d]+$/){
#		die("chains = $chains:$!");
#
#	}
#
#	my $flag =0; #0: the new pdb file is empty.
#
#	unlink $pdbFL_new if ( -e $pdbFL_new );
#	open( OUTPUT, ">>$pdbFL_new" ) || die("Cannot open $pdbFL_new:$!");
#
#	open( INPUT, "<$pdbFL_ori" ) || die("Cannot open $pdbFL_ori:$!");
#	while (<INPUT>) {
#		s/[\n\r]//mg;
#		if (/^(ATOM|HETATM)/) {
#
#			#http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
#			my $chnID = substr( $_, 21, 1 );
#
#			if ( $chnID =~ /^$chains/ ) {
#
#				print OUTPUT "$_\n";
#				$flag=1;
#			}
#
#		}
#
#	}
#	close INPUT;
#	close OUTPUT;
#
#
#
#	print
#"All template chains $chains from $pdbFL_ori are written into $pdbFL_new.\n";
#
#if($flag ==0){
#
#	die("$pdbFL_new is empty:$!");
#}
#
#}

#sub compareHomlogLsts_new {
#
#	#!usr/bin/perl -w
#	#Li Xue
#	#1/15/2013
#	#
#
##compare two homolog lists
##if they have the same pdbID and different chainID
##then this pdb complex is written to the output file
## perl ./compareHomlogLsts.pl HomologLst1 HomologLst2
##
##INPUT: two lists of homologs of two binding partners
##OUTPUT: a list of homo-complexes "../data/templates.lst"
##
##note: if one or two binding partners cannot find homologs, no output file will be generated.
##
##Usage example:
## perl ./compareHomlogLsts.pl ../data/1auiA_A_B/homologsOfquerySeq.lst ../data/1auiB_B_A/homologsOfquerySeq.lst
## perl ./compareHomlogLsts.pl ../data/1dkg_A_D/1dkgD_D_A/homologsOfquerySeq.lst ../data/1dkg_A_D/1dkgA_A_D/homologsOfquerySeq.lst
## perl ./compareHomlogLsts.pl ../data/trans208/2prg/2prg_B_C/2prgB_B_C/homologsOfquerySeq.lst ../data/trans208/2prg/2prg_B_C/2prgC_C_B/homologsOfquerySeq.lst
#	use strict;
#	use File::Basename;
#
#	my $inputFL1  = shift @_;
#	my $inputFL2  = shift @_;
#	my $outputDIR = shift @_;
#
#	my $query1 = basename( dirname($inputFL1) );
#	my $query2 = basename( dirname($inputFL2) );
#
#	my @Qry1homologs_in_homComplexes
#	  ;    #homologs of Qry chain 1 that forms homComplexes
#	my @Qry2homologs_in_homComplexes;
#	my @homComplexes;
#	my $homComplexesFL = "$outputDIR/templates.lst";
#	my @commonPDBIDs;
#
#	#----------------
#
#	print
#"\n\nCompare two homolog lists ($inputFL1 and $inputFL2) and write homo-interolog list....\n\n";
#
#	my $homologLst1 =
#	  &readPDBHomologLstFL($inputFL1);    #hash ref: $list1->{$pdbID} = @chainID
#	my $homologLst2 = &readPDBHomologLstFL($inputFL2);    #%hash ref
#
#	unlink($homComplexesFL) if ( -e $homComplexesFL );
#
#	if ( $homologLst1 && $homologLst2 ) {
#
#		#if both binding partners have homologs
#		my @pdbIDlst1 = keys(%$homologLst1);
#
#		print "\nPDB IDs in $inputFL1: @pdbIDlst1\n";
#
#		my @pdbIDlst2 = keys(%$homologLst2);
#		print "\nPDB IDs in $inputFL2: @pdbIDlst2\n";
#
#		my %seen;
#		if ( scalar @pdbIDlst1 >= scalar @pdbIDlst2 ) {
#			@seen{@pdbIDlst1} = (1) x $#pdbIDlst1;
#			@commonPDBIDs = grep { $seen{$_} } @pdbIDlst2;
#		}
#		else {
#			@seen{@pdbIDlst2} = (1) x $#pdbIDlst2;
#			@commonPDBIDs = grep { $seen{$_} } @pdbIDlst1;
#
#		}
#
#		if (@commonPDBIDs) {
#			print "\nCommon PDB IDs are: @commonPDBIDs\n";
#
#			foreach my $pdbid (<@commonPDBIDs>) {
#				my @homologChainIDs_Qry1 =
#				  &unique( @{ $homologLst1->{$pdbid} } );
#				my @homologChainIDs_Qry2 =
#				  &unique( @{ $homologLst2->{$pdbid} } );
#
#				foreach my $chainid1 (@homologChainIDs_Qry1) {
#					foreach my $chainid2 (@homologChainIDs_Qry2) {
#						if ( $chainid1 ne $chainid2 ) {
#
#							push @homComplexes,
#							  "$pdbid$chainid1:$pdbid$chainid2";
#
#							push @Qry1homologs_in_homComplexes,
#							  "$pdbid$chainid1\|$chainid1:$chainid2";
#							push @Qry2homologs_in_homComplexes,
#							  "$pdbid$chainid2\|$chainid2:$chainid1";
#
#						}
#					}
#				}
#
#			}
#
#			@Qry1homologs_in_homComplexes =
#			  &unique(@Qry1homologs_in_homComplexes);
#			@Qry2homologs_in_homComplexes =
#			  &unique(@Qry2homologs_in_homComplexes);
#			@homComplexes = @{ &uniquePair( \@homComplexes ) };
#
#			#2greB:2greO -> 2gre B:O
#			@homComplexes = @{ &reformHomComplexPairs( \@homComplexes ) };
#
#			&writeArrayIntoFL( \@homComplexes, $homComplexesFL );
#
#			print
#"compareHomlogLsts_new() is finished. $homComplexesFL is generated.  \n\n";
#
#		}
#		else {
#			my $qryComplex = basename( basename( dirname($inputFL1) ) );
#			print
#"compareHomlogLsts_new() is finished. $qryComplex does not have homoComplexes.\n";
#
#		}
#
#	}
#	else {
#
#		#at least one of the binding partners cannot find homologs
#		if ( !$homologLst1 ) {
#
#			#		../data/1i4fA_A_C/homologsOfquerySeq.lst
#			my ($partner) = $inputFL1 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
#			print
#"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";
#
#		}
#		if ( !$homologLst2 ) {
#
#			#		../data/1i4fC_C_A/homologsOfquerySeq.lst
#			my ($partner) = $inputFL2 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
#			print
#"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";
#
#		}
#	}
#
#	print "\n";
#
#	return \@commonPDBIDs;
#
#}

#sub transformHomointerologFL2HomologFL {
#
#	#Input (HomointerologFL):
#
##Homologs used in final prediction.
##	HomologID	logEval_H	PositiveS_H	Frac_Q	Frac_H	SimilarityScore
##	1tgsZ|Z:I	-108.578174314659	60.09	0.942857142857143	0.978165938864629	0.563818855368153
##	1hjaB|B:I	-189.505124806072	100	0.53469387755102	1	0.435053174254006
##	1choF|F:I	-189.505124806072	100	0.53469387755102	1	0.424073590325056
#
#	#Output(homologFL):
#	#	1tgsZ
#	#	1hjaB
#	#	1choF
#
#	my @homointerologFLs = @{ shift @_ };
#	my $outputFL         = shift @_;
#
#	my @homologs;
#
#	foreach my $homointerologFL (@homointerologFLs) {
#
#		if ( !-e $homointerologFL ) {
#			print
#"Warning: $homointerologFL does not exist. No template file can be generated.\n";
#			next;
#
#		}
#
#		open( INPUT, "<$homointerologFL" )
#		  || die("Cannot open $homointerologFL:$!");
#		while (<INPUT>) {
#
#			s/[\n\r]//mg;
#
#			if (/^(\w{5})\|\w+:\w+/) {
#
##	1tgsZ|Z:I	-108.578174314659	60.09	0.942857142857143	0.978165938864629	0.563818855368153
#				push @homologs, $1;
#			}
#
#		}
#		close INPUT;
#
#		@homologs = &unique(@homologs);
#
#		#--write to output file
#		unlink $outputFL if ( -e $outputFL );
#		open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
#
#		foreach (@homologs) {
#			print OUTPUT "$_\n";
#
#		}
#
#	}
#	close OUTPUT;
#
#	if ( -e $outputFL ) {
#
#		print "$outputFL generated.\n";
#	}
#	else {
#		print "No template used by PS-HomPPI. $outputFL cannot be generated.\n";
#	}
#
#}
#
#sub readRecLigLst {
#
#	#Input file:
#	#rec:lig pairs
#	#A:B
#	#C:B
#
#	my $recLigLst = shift @_;
#	my %receptors;
#	my %ligands;
#
#	open( INPUT, "<$recLigLst" ) || die("Cannot open $recLigLst:$!");
#	while (<INPUT>) {
#		s/[\n\r]//mg;
#
#		if (/^[\s\t]*(\w+:\w+)/) {
#
#			my ( $rec, $lig ) = split( /:/, $1 );
#			$receptors{$rec} = 1;
#			$ligands{$lig}   = 1;
#		}
#
#	}
#	close(INPUT);
#
#	my @rec = keys %receptors;
#	my @lig = keys %ligands;
#
#	if(!@rec ||!@lig){
#
#		die("No rec or lig read from $recLigLst:$!");
#
#	}
#
#	return ( \@rec, \@lig );
#
#}

#-----------------------------------------------
sub TemplateUsedByPSHomPPI_oneCase_new {

# Li Xue
# Nov. 11th, 2016
#
# Input files: statistics_wo_sameProt.txt
# OUTPUT 1:
# 	rec	rec	    lig
# 	A	B		C
# 	1awcA	1awcB	1awcE
#
#  OUTPUT 2: TemplatesUsed folder
#
#  NOTE: TemplateUsedByPSHomPPI_oneCase_new() and TemplateUsedByPSHomPPI_oneCase_old() differs
#  in the way to generate $templateFL = "$outputDIR_case/templates.lst".
#
#

    use strict;
    use PSHomPPI_resiPairs;
    use File::Copy;
    use File::Path;

    my $PSHomPPIoutputDIR_oneQryPair = shift @_;
    my $outputDIR_oneQryPair = "$PSHomPPIoutputDIR_oneQryPair/TemplatesUsed";
    my $pdbDIR = "$PSHomPPIoutputDIR_oneQryPair/TemplatesUsed/pdb";

    #--
    rmtree($outputDIR_oneQryPair) if (-d $outputDIR_oneQryPair);
    mkdir($outputDIR_oneQryPair) if ( !-d $outputDIR_oneQryPair );

    #--
    my $qry_pair = basename($PSHomPPIoutputDIR_oneQryPair);    #$pair='A:B'

    my $case_predictedIntDIR = dirname($PSHomPPIoutputDIR_oneQryPair);

    #-- for haddock dataset, $receptors is always chain A, and $ligands is always chain B.

    #--
    my $statFL = "$case_predictedIntDIR/statistics_wo_sameProt.txt"
    ; # This file will be used to extract sequence identity info of templates.

    if ( !-e $statFL ) {
        $statFL = "$case_predictedIntDIR/statistics.txt";
    }

    if (! -e $statFL){
        print "\nWARNING: $statFL does not exist. No homologs found for this case.\n";
        return;
    }

    #--
    my ($templateStats, $headers) =  &readStatisticsFL_all($statFL);

    my @templateIDs = keys % { $templateStats->{"$qry_pair"} };
    #-- @templateIDs = ('1ak4D,80:1ak4A,16', '1ak4D,80:1ak4A,37', ...)
    if (!@templateIDs){
        die("qry $qry_pair does not have templates in $statFL:$!");
    }

    #-- generate TemplatesUsed/templates_used_in_prediction.stat. Scan $statFL and only keep the templates in $templates_used_in_predictionFL.
    my $templates_used_in_predictionFL = "$PSHomPPIoutputDIR_oneQryPair/templates_used_in_prediction.lst";
    my $finalHomcomplexStatFL = "$outputDIR_oneQryPair/templates_used_in_prediction.stat";
    my @templateIDs_sorted = &genFinalHomcomplexStatFL( $qry_pair, $statFL, $templates_used_in_predictionFL, $finalHomcomplexStatFL ); #- templateIDs sorted by pred_CC


    #-- prepare template pdb files
    &prepareTemplatePDBfls ( \@templateIDs_sorted, $pdbDIR);
}

sub prepareTemplatePDBfls{
    #&prepareTemplatePDBfls ( \@templateIDs, $pdbDIR);


    use strict;
    my @templateIDs_sorted = @{shift @_}; # ('4lqwB,1:4lqwC,4', '1ak4B,16:1ak4C,51') They are sorted by pred_CC
    my $pdbDIR = shift @_; # output dir; $pdbDIR = "$outputDIR_oneQryPair/pdb";

    #-- remove group IDs from @templateIDs: '4lqwB,1:4lqwC,4' => '4lqwB:4lqwC'
    @templateIDs_sorted  = map {s/,\d+//g; $_} @templateIDs_sorted;
    my @templateIDs = @{&unique_keepOrder(\@templateIDs_sorted)};


    #-- extract template pdb IDs
    my @pdbIDs;
    foreach my $templateID (@templateIDs){
        # $templateID = '4lqwB,1:4lqwC,4'
        push @pdbIDs, substr($templateID, 0, 4);
    }
    @pdbIDs = &unique(@pdbIDs);

    if ( !@pdbIDs || scalar @pdbIDs == 0 ) {
        print "Warning: no homologs used in interface prediction. \n";
        return;
    }

    #-- 1. download pdb file of the templates
    mkdir $pdbDIR if ( !-d $pdbDIR );

    foreach my $templatePDBID (@pdbIDs) {

        $templatePDBID = lc($templatePDBID);

        if ( !-e "$pdbDIR/$templatePDBID.pdb" ) {
            my $command = "cp $unboundPDB_DIR/pdb$templatePDBID.ent.gz $pdbDIR";
            system("$command") == 0
              or die("Cannot copy pdb file $templatePDBID from ProtInDB to $pdbDIR. Command used: $command :$!");
            system("gunzip $pdbDIR/pdb$templatePDBID.ent.gz") == 0
              or die("Cannot unzip pdb$templatePDBID.ent.gz:$!");
            move( "$pdbDIR/pdb$templatePDBID.ent",
                "$pdbDIR/$templatePDBID.pdb" );

#                        system(
#            "wget  http://www.rcsb.org/pdb/files/$templatePDBID.pdb -P $pdbDIR"
#                          ) == 0
#                          or die("Cannot download $templatePDBID.pdb:$!");
#                        sleep 1;
        }
    }

    #-- 2. extract chains used as the templates from their pdb file

    for (my $i =0; $i< scalar @templateIDs; $i++){

        #-- @templateIDs = ('1m9xC:1m9xB', '1m9xG:1m9xF')
        #
        my $templateID = $templateIDs[$i]; # $templateID = '1m9xC,9:1m9xB,16'
        my ($prot1, $prot2) = split(/[,:]/, $templateID);
        my $pdbID = substr($prot1,0,4);
        my $chnID1 = substr($prot1, 4,1);
        my $chnID2 = substr($prot2, 4, 1);
        my $pdbFL = "$pdbDIR/$pdbID.pdb";
        my $n = $i+1;
        my $outPDBFL_rec = "$pdbDIR/template$n.rec.$prot1.pdb";
        my $outPDBFL_lig = "$pdbDIR/template$n.lig.$prot2.pdb";

        my $command ="pdb_selchain.py -$chnID1 $pdbFL > $outPDBFL_rec";
        system($command) ==0 or die ("FAILED: $command:$!");

         $command ="pdb_selchain.py -$chnID2 $pdbFL > $outPDBFL_lig";
        system($command) ==0 or die ("FAILED: $command:$!");

    }

}

#sub TemplateUsedByPSHomPPI_oneCase_new_ori {
#
##Li Xue
##Jan 14th, 2013
##
##Input files: homologs_used_in_predictions.lst for each chain of a query complex
##Output format:
##	rec	rec	    lig
##	A	B		C
##	1awcA	1awcB	1awcE
##
##
## TemplateUsedByPSHomPPI_oneCase_new() and TemplateUsedByPSHomPPI_oneCase_old() differs
## in the way to generate $templateFL = "$outputDIR_case/templates.lst".
##
##
##Usage:
##perl TemplateUsedByPSHomPPI.pl caseID PSHomPPIoutputDIR outputDIR_of_this_script
#
#    use strict;
#    use PSHomPPI_resiPairs;
#
#    my $PSHomPPIoutputDIR_oneQryPair = shift @_;
#    my $outputDIR_oneQryPair = "$PSHomPPIoutputDIR_oneQryPair/TemplatesUsed";
#
#    #--
#    mkdir($outputDIR_oneQryPair) if ( !-d $outputDIR_oneQryPair );
#
#    #--
#    my $pair = basename($PSHomPPIoutputDIR_oneQryPair);    #$pair='A_B'
#    my ( $rec, $lig ) = split( /:/, $pair );
#
#    my $case_predictedIntDIR = dirname($PSHomPPIoutputDIR_oneQryPair);
#
##-- for haddock dataset, $receptors is always chain A, and $ligands is always chain B.
#
#    #--
#    my $statFL = "$case_predictedIntDIR/statistics.txt"
#      ;    #statistics\_$caseID\_wo_sameProt.txt"
#    ; #statistics_1AVX_wo_sameProt.txt. This file will be used to extract sequence identity info of templates.
#
#    if ( !-e $statFL ) {
#        $statFL = "$case_predictedIntDIR/statistics_wo_sameProt.txt";
#
#    }
#
#    my %homologsUsedInPredictionFLs_rec;
#    my %homologsUsedInPredictionFLs_lig;
#    my @pdbIDs;
#    print "Qry ID pair: $pair\n\n";
#    my $pairDIR = $PSHomPPIoutputDIR_oneQryPair;
#
#    my $homologsUsedInPredictionFL_rec =
#      "$pairDIR/$rec/homointerologs_used_in_prediction.txt";
#    push @{ $homologsUsedInPredictionFLs_rec{$rec} },
#      $homologsUsedInPredictionFL_rec;
#    my $homologsUsedInPredictionFL_lig =
#      "$pairDIR/$lig/homointerologs_used_in_prediction.txt";
#    push @{ $homologsUsedInPredictionFLs_lig{$lig} },
#      $homologsUsedInPredictionFL_lig;
#
#    #--
#    my $homComplexUsedInPredictionFL =
#      "$pairDIR/templates_used_in_prediction.lst"
#      ;    #output of genHomcomplexFL
#    @pdbIDs = &genHomcomplexFL_2(
#        $rec, $lig,
#        $homologsUsedInPredictionFL_rec,
#        $homologsUsedInPredictionFL_lig,
#        $statFL, $homComplexUsedInPredictionFL
#    );
#
#    if ( -e $homComplexUsedInPredictionFL ) {
#        mkdir($outputDIR_oneQryPair) if ( !-d $outputDIR_oneQryPair );
#
#        copy( $homComplexUsedInPredictionFL, $outputDIR_oneQryPair )
#          || die(
#"Cannot move $homComplexUsedInPredictionFL to $outputDIR_oneQryPair:$!"
#          );
#        print
#          "$homComplexUsedInPredictionFL is copied to $outputDIR_oneQryPair.\n";
#    }
#
#    #--
#    if ( !@pdbIDs || scalar @pdbIDs == 0 ) {
#
#        print "Warning: no homologs used in interface prediction. \n";
#        return;
#    }
#
#    #--download pdb file of the templates
#    my $pdbDIR = "$outputDIR_oneQryPair/pdb";
#    mkdir $pdbDIR if ( !-d $pdbDIR );
#
#    foreach my $templatePDBID (@pdbIDs) {
#
#        #templatePDBID=4j2yB:4j2yA
#
#        $templatePDBID = lc($templatePDBID);
#
#        if ( !-e "$pdbDIR/$templatePDBID.pdb" ) {
#            my $command = "cp $unboundPDB_DIR/pdb$templatePDBID.ent.gz $pdbDIR";
#            system("$command") == 0
#              or die("Cannot copy pdb file $templatePDBID from ProtInDB to $pdbDIR. Command used: $command :$!");
#            system("gunzip $pdbDIR/pdb$templatePDBID.ent.gz") == 0
#              or die("Cannot unzip pdb$templatePDBID.ent.gz:$!");
#            move( "$pdbDIR/pdb$templatePDBID.ent",
#                "$pdbDIR/$templatePDBID.pdb" );
#
##                        system(
##            "wget  http://www.rcsb.org/pdb/files/$templatePDBID.pdb -P $pdbDIR"
##                          ) == 0
##                          or die("Cannot download $templatePDBID.pdb:$!");
##                        sleep 1;
#        }
#    }
#
#    #--extract chains used as the templates from their pdb file
#    my $templateFL =
#      "$outputDIR_oneQryPair/templates_used_in_prediction.stat";
#
#    if ( !-e $templateFL ) {
#
#        print
#"No hom-complexexes found by PS-HomPPI. No template file can be written!\n";
#        return;
#    }
#    &prepareRecLigTemplatePDB( $templateFL, $pdbDIR );
#}

sub genFinalHomcomplexStatFL {
    #-- scan $statFL and only keep the rows with homologs in P:Q/templates_used_in_prediction.lst
    #-- output: TemplatesUsed/templates_used_in_prediction.stat (templates sorted by pred_CC)
    #
    # INPUT 1 (statistics.txt):
    #
    #    #qryPairID homologPair pred_CC len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2    EVal1   EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL_Q1    LAL_H1  LAL_Q2  LAL_H2  frac_LAL1   frac_LAL2   TP_1    TN_1    FP_1    FN_1CC_1    specificity_1   sensitivity_1   TP_2 TN_2    FP_2    FN_2    CC_2    specificity_2   sensitivity_2
    #    1m9xE:1m9xH 4lqwA,1:4lqwC,26    0.575793151365371   165 146 177 146 14  9   9   11  1.5E-47 2.3E-72 66.2    97.9    88.9    98.6   3.4 4.3 100.0   100.0   163 163 146 146 0.91    1.002142    7   12  0.12    0.22    0.14    0   115 11  9   -0.08   0.00    0.00
    #
    # INPUT 2 (templates_used_in_prediction.lst):
    #
    #    #Template_ID predicted_CC
    #    2x83B,13:2x83A,26   0.993
    #    1m9xA,50:1m9xD,86   0.952
    #
    # Output format is the same as $statFL:
    #
    #   #qryPairID homologPair pred_CC len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2    EVal1   EVal2   SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL_Q1    LAL_H1  LAL_Q2  LAL_H2  frac_LAL1   frac_LAL2   TP_1    TN_1    FP_1    FN_1CC_1    specificity_1   sensitivity_1   TP_2 TN_2    FP_2    FN_2    CC_2    specificity_2   sensitivity_2
    #   1m9xE:1m9xH 4lqwA,1:4lqwC,26    0.575793151365371   165 146 177 146 14  9   9   11  1.5E-47 2.3E-72 66.2    97.9    88.9    98.6   3.4 4.3 100.0   100.0   163 163 146 146 0.91    1.002142    7   12  0.12    0.22    0.14    0   115 11  9   -0.08   0.00    0.00


    use PSHomPPI_resiPairs;

    my $qryPairID = shift @_;
    my  $statFL= shift @_;
    my  $templates_used_in_predictionFL = shift @_;
    my  $outFL = shift @_; #-- .../TemplatesUsed/templates_used_in_prediction.stat

    #---- read template IDs in the templates_used_in_prediction.lst
    my $predCCs = &readTemplates_used_in_predictionFL($templates_used_in_predictionFL);
    if (!defined $predCCs){
        die("nothing read from $templates_used_in_predictionFL:$!");
    }
    my @templateIDs_sorted = sort { $predCCs->{$b} cmp $predCCs->{$a} } keys %{$predCCs};


    #--- scan $statFL and only keep the rows with homologs in $templates_used_in_predictionFL
    my ($templateStats, $headers) = &readStatisticsFL_all($statFL);
    #-  $templateStats->{$qryPairID}->{$templateID}->{$header} = $stats[$i];


    unlink $outFL if (-e $outFL);
    open (OUTPUT, ">> $outFL") or die ("Cannot open $outFL:$!");

    #-- print header line
    my $tmp = join("\t", @{$headers});
    print OUTPUT "#qryPairID\ttemplateID\t$tmp\n";

    #-- print content
    foreach my $templateID (@templateIDs_sorted){

        my @line ;
        push @line, "$qryPairID";
        push @line, "$templateID";
        foreach my $header (@$headers){
            push @line, $templateStats->{$qryPairID}->{$templateID}->{$header};
        }
        my $LINE = join ("\t", @line);
        print OUTPUT "$LINE\n";
    }
    close OUTPUT;

    print "$outFL generated.\n";

    return @templateIDs_sorted;

}

#sub readTemplates_used_in_predictionFL{
#
#  # INPUT: templates_used_in_prediction.lst
#  #
#  #    Template_ID  predicted_CC
#  #    1m9xB,50:1m9xC,9    0.998
#  #    4lqwA,1:4lqwD,77    0.986
#
#  my $homcomplexes_used_in_predictionFL = shift @_;
#  my $predCCs;
#
#  open (INPUT, "<$homcomplexes_used_in_predictionFL") or die ("Cannot open $homcomplexes_used_in_predictionFL:$!");
#  while (<INPUT>){
#      s/[\n\r]//gm;
#
#      if (/^\w+/){
#          my ($templateID, $pred_CC) = split( /\s+/, $_ );
#          $predCCs->{$templateID} = $pred_CC;
#      }
#  }
#  close INPUT;

#  if (!defined $predCCs){
#      die("nothing read from $homcomplexes_used_in_predictionFL:$!");
#  }

#  return $predCCs;
#}

#sub genHomcomplexFL_3 {
#
##    @pdbIDs = &genHomcomplexFL_3(
##        $statFL, $homComplexUsedInPredictionFL
##    );
#
#    # Read all the A':B' in the statFL and use them for later predictions.
#    #
#    # INPUT 1 ( $statFL):
#    #
#    # #Homologs listed in ../../uploadData/test/delete.lst are removed from ../../uploadData/test/statistics.txt
#    # qryPairID  homologPair len_Qry1    len_Qry2    len_H1  len_H2  numInt_Q1   numInt_Q2   numInt_Homolog1 numInt_Homolog2 EVal1   EVal2  SID1    SID2    Positives1  Positives2  Similarity1 Similarity2 hhpredProb1 hhpredProb2 LAL1    LAL2    frac_LAL1 frac_LAL2   TP_1    TN_1    FP_1    FN_1    CC_1    specificity_1sensitivity_1  TP_2    TN_2    FP_2    FN_2    CC_2 specificity_2   sensitivity_2
#    # Q:P 2x2dC,42:2x2dD,69   164 137 165 147 Nan Nan 14  8   1.1E-45 3E-41   99.3    98.599.3    100.0   4.7 4.4 100.0   100.0   164    136 0.99    0.92
#    # Q:P 2x2dC,42:2x2dD,38   164 137 165 147 Nan Nan 14  8   1.1E-45 3.7E-68 99.3    98.599.3    100.0   4.7 4.4 100.0   100.0   164    136 0.99    0.92
#
#    #
#    # OUTPUT:
#    #	#pdbID	rec:lig
#    #	3rdz	B:D
#    #	3rdz	A:C
#
#    my (
#        $statFL, $outputFL
#    ) = @_;
#
#    print "Write hom-complexes used in predictions into FL ...\n\n";
#
#    my @pdbIDs;
#
#    if ( !-e $statFL ) {
#        print "\nWarning: $statFL does not exist. No homologs found for this case.\n";
#        return @pdbIDs;
#    }
#
#    #-- read statFL
#
#    my (
#        $Num_Int,      $Length1,     $Length2,      $logEval_H,
#        $logEval_HP,   $PositiveS_H, $PositiveS_HP, $IdentityS_H,
#        $IdentityS_HP, $alignLen_H,  $alignLen_HP,  $Frac_Q,
#        $Frac_H,       $Frac_P,      $Frac_HP,
#      )
#      = &readStatisticsFL_all($statFL)
#      ;    #-note: @{ $PositiveS_H->{$queryID}->{$homologID} }
#
#    #--write $homComplex to output file
#    @homcomplexes = &unique(@homcomplexes);
#
#    my $Len_recA = $Length1->{$recID};
#    my $Len_ligB = $Length1->{$ligID};
#
#    unlink $outputFL if ( -e $outputFL );
#    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
#    print OUTPUT "#generated by genSuperimposePDBFL.pl\n";
#    print OUTPUT "#A: rec. B: lig. A': homolog of A. B': homolog of B.\n";
#
##print OUTPUT "#pdbID\trecA:ligB\tlocal_IdentityS_AA'\tAlocal_IdentityS_BB'\tLALAA'\tLALBB'\tLenA\tLenB\tglobal_IdentitySAA'\tglobal_IdentitySBB'\tave_similarity_score\n";
#    print OUTPUT
#"# global_IdentityS_rec = ( num_identicalResiNum_AAprime / Len_recA + num_identicalResiNum_AAprime / Len_Aprime ) / 2 * 100\n";
#    print OUTPUT
#"# global_IdentityS_lig = ( num_identicalResiNum_BBprime / Len_ligB + num_identicalResiNum_BBprime / Len_Bprime ) / 2 * 100\n";
#    print OUTPUT
#"#pdbID\trecA:ligB\tlocal_IdentityS_AA'\tlocal_IdentityS_BB'\tglobal_IdentityS_AA'\tglobal_IdentityS_BB'\n";
#
#    #note:A-rec, B-lig
#    #
#
#    my %seen;
#
#    foreach my $templateID (@$homcomplexes_rec) {
#
#        #$tempalteID= 'pdbID	rec:lig'
#        $seen{$templateID} = 1;
#
#        print "Template ID for rec: $templateID.\n";
#
#        my ( $pdbID, $homologChnID_rec, $homologChnID_lig ) =
#          split( /[\s\t:]/, $templateID );
#
#        push @pdbIDs, $pdbID;
#
##my $aveSimilarityScore = &getAveSimilarityScore($templateScore,$homologChnID_rec, $homologChnID_lig);
#        my $queryID = "$recID|$recID:$ligID";
#
#        my $homolog_rec =
#          "$pdbID$homologChnID_rec|$homologChnID_rec:$homologChnID_lig";
#
#        my $local_IdentityS_rec =
#          &MeanArray( $IdentityS_H->{$queryID}->{$homolog_rec} );
#        my $local_IdentityS_lig =
#          &MeanArray( $IdentityS_HP->{$queryID}->{$homolog_rec} );
#
#        #print "local_IdentityS_rec: $local_IdentityS_rec\n";
#        #print "local_IdentityS_lig:$local_IdentityS_lig\n";
#
#        #note:A-rec, B-lig
#        my $LAL_AAprime = &MeanArray( $alignLen_H->{$queryID}->{$homolog_rec} )
#          ;    #LAL between rec A and A'
#        my $LAL_BBprime =
#          &MeanArray( $alignLen_HP->{$queryID}->{$homolog_rec} );
#
#        #note:A-rec, B-lig
#        my $Len_Aprime = $Length1->{"$pdbID$homologChnID_rec"};
#        my $Len_Bprime = $Length1->{"$pdbID$homologChnID_lig"};
#
#        print "LAL_BBprime: $LAL_BBprime\n";
#
#        if ( !defined $Len_Bprime ) {
#
#            die(
#"length of homolog of ligand \(B = $queryID, B\' = $pdbID$homologChnID_lig \) is not in $statFL:$!"
#            );
#        }
#
#        #note:A-rec, B-lig
#        my $num_identicalResiNum_AAprime =
#          $local_IdentityS_rec * $LAL_AAprime * 0.01;
#        my $num_identicalResiNum_BBprime =
#          $local_IdentityS_lig * $LAL_BBprime * 0.01;
#
#        #note:A-rec, B-lig
#        my $global_IdentityS_rec = ( $num_identicalResiNum_AAprime / $Len_recA +
#              $num_identicalResiNum_AAprime / $Len_Aprime ) / 2 * 100;
#        my $global_IdentityS_lig = ( $num_identicalResiNum_BBprime / $Len_ligB +
#              $num_identicalResiNum_BBprime / $Len_Bprime ) / 2 * 100;
#
##print "global_IdentityS_lig = ($num_identicalResiNum_BBprime/$Len_ligB + $num_identicalResiNum_BBprime/$Len_Bprime)/2\n";
#
##print OUTPUT "$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$LAL_AAprime\t$LAL_BBprime\t$Len_recA\t$Len_ligB\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#        print OUTPUT
#"$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#    }
#    close OUTPUT;
#
#    foreach my $templateID (@$homcomplexes_lig) {
#
#        #$tempalteID= 'pdbID	lig:rec'
#
#        print "Template ID for lig: $templateID.\n";
#
#        my ( $pdbID, $homologChnID_lig, $homologChnID_rec ) =
#          split( /[\s\t:]/, $templateID );
#
#        if (   $seen{$templateID}
#            || $seen{"$pdbID\t$homologChnID_rec:$homologChnID_lig"} )
#        {
#
##The information of this template has been written in the output file using the for loop of @$homcomplexes_rec
#            print
#"The information of this template has been written in the output file using the for loop of \@\$homcomplexes_rec\n";
#            next;
#        }
#
#        push @pdbIDs, $pdbID;
#
##my $aveSimilarityScore = &getAveSimilarityScore($templateScore,$homologChnID_rec, $homologChnID_lig);
#        my $queryID = "$ligID|$ligID:$recID";
#
#        my $homolog_lig =
#          "$pdbID$homologChnID_lig|$homologChnID_lig:$homologChnID_rec";
#
#        my $local_IdentityS_lig =
#          &MeanArray( $IdentityS_H->{$queryID}->{$homolog_lig} );
#        my $local_IdentityS_rec =
#          &MeanArray( $IdentityS_HP->{$queryID}->{$homolog_lig} );
#
#        #note:A-rec, B-lig
#        my $LAL_AAprime = &MeanArray( $alignLen_HP->{$queryID}->{$homolog_lig} )
#          ;    #LAL between rec A and A'
#        my $LAL_BBprime = &MeanArray( $alignLen_H->{$queryID}->{$homolog_lig} );
#
#        #note:A-rec, B-lig
#        my $Len_Aprime = $Length1->{"$pdbID$homologChnID_rec"};
#        my $Len_Bprime = $Length1->{"$pdbID$homologChnID_lig"};
#
#        #note:A-rec, B-lig
#        my $num_identicalResiNum_AAprime =
#          $local_IdentityS_rec * $LAL_AAprime * 0.01;
#        my $num_identicalResiNum_BBprime =
#          $local_IdentityS_lig * $LAL_BBprime * 0.01;
#
#        #note:A-rec, B-lig
#        my $global_IdentityS_rec = ( $num_identicalResiNum_AAprime / $Len_recA +
#              $num_identicalResiNum_AAprime / $Len_Aprime ) / 2 * 100;
#        my $global_IdentityS_lig = ( $num_identicalResiNum_BBprime / $Len_ligB +
#              $num_identicalResiNum_BBprime / $Len_Bprime ) / 2 * 100;
#
##print "global_IdentityS_lig = ($num_identicalResiNum_BBprime/$Len_ligB + $num_identicalResiNum_BBprime/$Len_Bprime)/2\n";
#
##print OUTPUT "$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$LAL_AAprime\t$LAL_BBprime\t$Len_recA\t$Len_ligB\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#        print OUTPUT
#"$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#    }
#    close OUTPUT;
#
#    if ( -e $outputFL ) {
#        print "$outputFL generated.\n";
#    }
#    else {
#        die("$outputFL does not exist:$!");
#    }
#
#    #--
#    @pdbIDs = &unique(@pdbIDs);
#    return @pdbIDs;
#
#}
#
#
#sub genHomcomplexFL_2 {
#
#    #Input ($homologsUsedInPredictionFL):
#    #HomologID	logEval_H	PositiveS_H	Frac_Q	Frac_H	SimilarityScore
#    #3rdzB|B:D	-107.52835219016	59.91	0.938775510204082	1	0.390479103445849
#    #3rdzA|A:C	-107.52835219016	59.91	0.938775510204082	1	0.390479103445849
#
#    #Input $statFL:
#    #>QUERY PDBID + CHAINID|CHAINID:INTERACTING CHAINID
#    ##HOMOLOG-PDBID+CHAINID|CHAINID:INTERACTING CHAINID	Num_Int	num_residue1	num_residue2	num_residue1_P	Bit_score	Eval_H	Eval_HP	PositiveS_H	PositiveS_HP	IdentityS_H	IdentityS_HP	alignLen_H	alignLen_HP	Frac_Q	Frac_H	Frac_P	Frac_HP	TP	TN	FP	FN	CC	Specificity	Sensitivity	Accuracy
##>A|A:B	33	223	223	169
##2qyiC|C:	37	223	223	183	349	1e-123	2e-30	93.27	58.29	81.17	46.29	223	175	1	1	0.96	0.91	29	182	8	4	0.798614482572727	0.783783783783784	0.878787878787879	0.946188340807175
##2qyiA|A:B	40	223	223	183	349	1e-123	2e-30	93.27	58.29	81.17	46.29	223	175	1	1	0.96	0.91	30	180	10	3	0.79265640773292	0.75	0.909090909090909	0.94170403587444
##
##
##Output:
##	#pdbID	rec:lig
##	3rdz	B:D
##	3rdz	A:C
#
#    my $recID = shift @_;
#    my $ligID = shift @_;
#    my (
#        $homologsUsedInPredictionFL_rec,
#        $homologsUsedInPredictionFL_lig,
#        $statFL, $outputFL
#    ) = @_;
#
#    print "Write hom-complexes used in predictions into FL ...\n\n";
#    print "homologsUsedInPredictionFL_rec:$homologsUsedInPredictionFL_rec\n";
#    print "homologsUsedInPredictionFL_lig:$homologsUsedInPredictionFL_lig\n\n";
#
#    my @pdbIDs;
#
#    if ( !-e $homologsUsedInPredictionFL_rec ) {
#
#        print
#"\nWarning: $homologsUsedInPredictionFL_rec does not exist. No homologs used in interface prediction.\n";
#        return @pdbIDs;
#    }
#
#    if ( !-e $homologsUsedInPredictionFL_lig ) {
#
#        print
#"\nWarning: $homologsUsedInPredictionFL_lig does not exist. No homologs used in interface prediction.\n";
#        return @pdbIDs;
#    }
#
#    if ( !-e $statFL ) {
#
#        print
#"\nWarning: $statFL does not exist. No homologs found for this case.\n";
#        return @pdbIDs;
#    }
#
#    my $homComplex
#      ; #$homComplex->{homologChn_of_rec}->{homologChn_of_lig} = "pdbID\thomologChn_of_rec:homologChn_of_lig";
#    my $RecTemplateScore;
#    my $LigTemplateScore;
#
#    #--
#
#    my ( $homcomplexes_rec, $templatePDBID_rec ) =
#      &readHomologsUsedInPredictionFL($homologsUsedInPredictionFL_rec)
#      ;    #templates used for predicting int of rec
#    my ( $homcomplexes_lig, $templatePDBID_lig ) =
#      &readHomologsUsedInPredictionFL($homologsUsedInPredictionFL_rec)
#      ;    #templates used for predicting int of lig
#           #--note: push @homcomplexes, "$pdbID\t$chn1:$chn2";
#
#    #--combine @homcomplexes_rec and @homcomplexes_lig
#    my @homcomplexes = @$homcomplexes_rec;
#
#    foreach my $ID (@$homcomplexes_lig) {
#
#        #$ID = "1acb	homolog_of_lig:homolog_of_rec
#        my ( $pdbID, $homologChnID_lig, $homologChnID_rec ) =
#          split( /[\t:]/, $ID );
#        push @homcomplexes, "$pdbID\t$homologChnID_rec:$homologChnID_lig";
#
#    }
#
#    #-- read statFL
#
#    my (
#        $Num_Int,      $Length1,     $Length2,      $logEval_H,
#        $logEval_HP,   $PositiveS_H, $PositiveS_HP, $IdentityS_H,
#        $IdentityS_HP, $alignLen_H,  $alignLen_HP,  $Frac_Q,
#        $Frac_H,       $Frac_P,      $Frac_HP,
#      )
#      = &readStatisticsFL_all($statFL)
#      ;    #-note: @{ $PositiveS_H->{$queryID}->{$homologID} }
#
#    #--write $homComplex to output file
#    @homcomplexes = &unique(@homcomplexes);
#
#    my $Len_recA = $Length1->{$recID};
#    my $Len_ligB = $Length1->{$ligID};
#
#    unlink $outputFL if ( -e $outputFL );
#    open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
#    print OUTPUT "#generated by genSuperimposePDBFL.pl\n";
#    print OUTPUT "#A: rec. B: lig. A': homolog of A. B': homolog of B.\n";
#
##print OUTPUT "#pdbID\trecA:ligB\tlocal_IdentityS_AA'\tAlocal_IdentityS_BB'\tLALAA'\tLALBB'\tLenA\tLenB\tglobal_IdentitySAA'\tglobal_IdentitySBB'\tave_similarity_score\n";
#    print OUTPUT
#"# global_IdentityS_rec = ( num_identicalResiNum_AAprime / Len_recA + num_identicalResiNum_AAprime / Len_Aprime ) / 2 * 100\n";
#    print OUTPUT
#"# global_IdentityS_lig = ( num_identicalResiNum_BBprime / Len_ligB + num_identicalResiNum_BBprime / Len_Bprime ) / 2 * 100\n";
#    print OUTPUT
#"#pdbID\trecA:ligB\tlocal_IdentityS_AA'\tlocal_IdentityS_BB'\tglobal_IdentityS_AA'\tglobal_IdentityS_BB'\n";
#
#    #note:A-rec, B-lig
#    #
#
#    my %seen;
#
#    foreach my $templateID (@$homcomplexes_rec) {
#
#        #$tempalteID= 'pdbID	rec:lig'
#        $seen{$templateID} = 1;
#
#        print "Template ID for rec: $templateID.\n";
#
#        my ( $pdbID, $homologChnID_rec, $homologChnID_lig ) =
#          split( /[\s\t:]/, $templateID );
#
#        push @pdbIDs, $pdbID;
#
##my $aveSimilarityScore = &getAveSimilarityScore($templateScore,$homologChnID_rec, $homologChnID_lig);
#        my $queryID = "$recID|$recID:$ligID";
#
#        my $homolog_rec =
#          "$pdbID$homologChnID_rec|$homologChnID_rec:$homologChnID_lig";
#
#        my $local_IdentityS_rec =
#          &MeanArray( $IdentityS_H->{$queryID}->{$homolog_rec} );
#        my $local_IdentityS_lig =
#          &MeanArray( $IdentityS_HP->{$queryID}->{$homolog_rec} );
#
#        #print "local_IdentityS_rec: $local_IdentityS_rec\n";
#        #print "local_IdentityS_lig:$local_IdentityS_lig\n";
#
#        #note:A-rec, B-lig
#        my $LAL_AAprime = &MeanArray( $alignLen_H->{$queryID}->{$homolog_rec} )
#          ;    #LAL between rec A and A'
#        my $LAL_BBprime =
#          &MeanArray( $alignLen_HP->{$queryID}->{$homolog_rec} );
#
#        #note:A-rec, B-lig
#        my $Len_Aprime = $Length1->{"$pdbID$homologChnID_rec"};
#        my $Len_Bprime = $Length1->{"$pdbID$homologChnID_lig"};
#
#        print "LAL_BBprime: $LAL_BBprime\n";
#
#        if ( !defined $Len_Bprime ) {
#
#            die(
#"length of homolog of ligand \(B = $queryID, B\' = $pdbID$homologChnID_lig \) is not in $statFL:$!"
#            );
#        }
#
#        #note:A-rec, B-lig
#        my $num_identicalResiNum_AAprime =
#          $local_IdentityS_rec * $LAL_AAprime * 0.01;
#        my $num_identicalResiNum_BBprime =
#          $local_IdentityS_lig * $LAL_BBprime * 0.01;
#
#        #note:A-rec, B-lig
#        my $global_IdentityS_rec = ( $num_identicalResiNum_AAprime / $Len_recA +
#              $num_identicalResiNum_AAprime / $Len_Aprime ) / 2 * 100;
#        my $global_IdentityS_lig = ( $num_identicalResiNum_BBprime / $Len_ligB +
#              $num_identicalResiNum_BBprime / $Len_Bprime ) / 2 * 100;
#
##print "global_IdentityS_lig = ($num_identicalResiNum_BBprime/$Len_ligB + $num_identicalResiNum_BBprime/$Len_Bprime)/2\n";
#
##print OUTPUT "$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$LAL_AAprime\t$LAL_BBprime\t$Len_recA\t$Len_ligB\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#        print OUTPUT
#"$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#    }
#    close OUTPUT;
#
#    foreach my $templateID (@$homcomplexes_lig) {
#
#        #$tempalteID= 'pdbID	lig:rec'
#
#        print "Template ID for lig: $templateID.\n";
#
#        my ( $pdbID, $homologChnID_lig, $homologChnID_rec ) =
#          split( /[\s\t:]/, $templateID );
#
#        if (   $seen{$templateID}
#            || $seen{"$pdbID\t$homologChnID_rec:$homologChnID_lig"} )
#        {
#
##The information of this template has been written in the output file using the for loop of @$homcomplexes_rec
#            print
#"The information of this template has been written in the output file using the for loop of \@\$homcomplexes_rec\n";
#            next;
#        }
#
#        push @pdbIDs, $pdbID;
#
##my $aveSimilarityScore = &getAveSimilarityScore($templateScore,$homologChnID_rec, $homologChnID_lig);
#        my $queryID = "$ligID|$ligID:$recID";
#
#        my $homolog_lig =
#          "$pdbID$homologChnID_lig|$homologChnID_lig:$homologChnID_rec";
#
#        my $local_IdentityS_lig =
#          &MeanArray( $IdentityS_H->{$queryID}->{$homolog_lig} );
#        my $local_IdentityS_rec =
#          &MeanArray( $IdentityS_HP->{$queryID}->{$homolog_lig} );
#
#        #note:A-rec, B-lig
#        my $LAL_AAprime = &MeanArray( $alignLen_HP->{$queryID}->{$homolog_lig} )
#          ;    #LAL between rec A and A'
#        my $LAL_BBprime = &MeanArray( $alignLen_H->{$queryID}->{$homolog_lig} );
#
#        #note:A-rec, B-lig
#        my $Len_Aprime = $Length1->{"$pdbID$homologChnID_rec"};
#        my $Len_Bprime = $Length1->{"$pdbID$homologChnID_lig"};
#
#        #note:A-rec, B-lig
#        my $num_identicalResiNum_AAprime =
#          $local_IdentityS_rec * $LAL_AAprime * 0.01;
#        my $num_identicalResiNum_BBprime =
#          $local_IdentityS_lig * $LAL_BBprime * 0.01;
#
#        #note:A-rec, B-lig
#        my $global_IdentityS_rec = ( $num_identicalResiNum_AAprime / $Len_recA +
#              $num_identicalResiNum_AAprime / $Len_Aprime ) / 2 * 100;
#        my $global_IdentityS_lig = ( $num_identicalResiNum_BBprime / $Len_ligB +
#              $num_identicalResiNum_BBprime / $Len_Bprime ) / 2 * 100;
#
##print "global_IdentityS_lig = ($num_identicalResiNum_BBprime/$Len_ligB + $num_identicalResiNum_BBprime/$Len_Bprime)/2\n";
#
##print OUTPUT "$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$LAL_AAprime\t$LAL_BBprime\t$Len_recA\t$Len_ligB\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#        print OUTPUT
#"$templateID\t$local_IdentityS_rec\t$local_IdentityS_lig\t$global_IdentityS_rec\t$global_IdentityS_lig\n";
#    }
#    close OUTPUT;
#
#    if ( -e $outputFL ) {
#        print "$outputFL generated.\n";
#    }
#    else {
#        die("$outputFL does not exist:$!");
#    }
#
#    #--
#    @pdbIDs = &unique(@pdbIDs);
#    return @pdbIDs;
#
#}
#
#----------------------------

#sub readHomologsUsedInPredictionFL {
#
#    my $homologsUsedInPredictionFL = shift @_;
#
#    my @homcomplexes;
#    my @pdbIDs;
#
#    open( INPUT, "<$homologsUsedInPredictionFL" )
#      || die("Cannot open $homologsUsedInPredictionFL:$!");
#    while (<INPUT>) {
#
#        s/[\n\r]//mg;
#
#        if (/^(\w{5})\|(\w+):(\w+)/) {
#
##	1tgsZ|Z:I	-108.578174314659	60.09	0.942857142857143	0.978165938864629	0.563818855368153
#            my $homolog = $1;
#            my $chn1    = $2;                         #homologChn_of_rec
#            my $chn2    = $3;                         #homologChn_of_lig
#            my $pdbID   = substr( $homolog, 0, 4 );
#
#            my ( $template, $logEval, $positiveS, $FracQ, $fracH,
#                $similarityScore )
#              = split( /[\s\t]+/, $_ );
#
##$homComplex->{$pdbID}->{$chn1}->{$chn2} =
##  "$pdbID\t$chn1:$chn2"
##  ; #always put the homologs of the rec ($chn1) before the homologs of lig ($chn2)
##$RecTemplateScore->{$template}->{similarityScore} =
##  $similarityScore
##  ; #always put the homologs of the rec ($chn2) before the homologs of lig ($chn1)
##$RecTemplateScore->{$template}->{positiveS} = $positiveS;
##$RecTemplateScore->{$template}->{FracQ}     = $FracQ;
#            push @homcomplexes, "$pdbID\t$chn1:$chn2";
#            push @pdbIDs,       $pdbID;
#
#        }
#
#    }
#    close INPUT;
#
#    if ( !@homcomplexes || !@pdbIDs ) {
#        die("no templates read from $homologsUsedInPredictionFL:$!");
#    }
#
#    return ( \@homcomplexes, \@pdbIDs );
#
#}
#
sub readTemplateFL {

    #generated by genSuperimposePDBFL.pl
    #pdbID  rec:lig
    #3rdz    B:D
    #3rdz    A:C

    my $templateFL = shift @_;
    my @templates;

    open( INPUT, "<$templateFL" ) || die("Cannot open $templateFL:$!");
    while (<INPUT>) {

        s/[\n\r]//gm;

        if (/^(\w{4}[\s\t]+\w+:\w+)/) {

            push @templates, $1;
        }

    }
    close INPUT;
    return \@templates;

}

#------------------------
sub prepareRecLigTemplatePDB {

    #input template FL
    #1hja	AC:IM
    #1ldt	TS:L

    #write chains on one side of the colon into one file
    #Left side of the colon are template chains for receptors
    #Right side for ligands.

    print "\n\n\nprepare PDB file of templates ...\n\n";

    my ( $templateFL, $pdbDIR ) = @_;

    my @templates = @{ &readTemplateFL($templateFL) };

    if ( scalar @templates == 0 ) {
        die("No templates read from $templateFL:$!");

    }
    my $templateNo = 0;

    foreach my $template (@templates) {

        #	$template = '1hja	AC:IM'

        print "prepare PDB file of template $template ...\n";

        $templateNo++;

        my ( $pdbID, $rec_template, $lig_template ) =
          split( /[\s\t:]+/, $template );

        my $pdbFL_ori = "$pdbDIR/$pdbID.pdb";
        my $pdbFL_new =
          "$pdbDIR/template$templateNo.rec.$pdbID$rec_template.pdb";
        &extractPartPDBFL( $pdbFL_ori, $rec_template, $pdbFL_new );

        $pdbFL_new = "$pdbDIR/template$templateNo.lig.$pdbID$lig_template.pdb";
        &extractPartPDBFL( $pdbFL_ori, $lig_template, $pdbFL_new );

    }
    print "\nTemplate Number is based on the order in $templateFL\n";

}

sub extractPartPDBFL {

    #extract the part of $chains out of $pdbFL_ori

    my $pdbFL_ori = shift @_;
    my $chains    = shift @_;    #'AC'
    my $pdbFL_new = shift @_;

    print "Extract the part of chain $chains out of $pdbFL_ori ...\n";

    if ( $chains !~ /^[a-zA-Z\d]+$/ ) {
        die("chains = $chains:$!");

    }

    my $flag = 0;                #0: the new pdb file is empty.

    unlink $pdbFL_new if ( -e $pdbFL_new );
    open( OUTPUT, ">>$pdbFL_new" ) || die("Cannot open $pdbFL_new:$!");

    open( INPUT, "<$pdbFL_ori" ) || die("Cannot open $pdbFL_ori:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;
        if (/^(ATOM|HETATM)/) {

            #http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
            my $chnID = substr( $_, 21, 1 );

            if ( $chnID =~ /^$chains/ ) {

                print OUTPUT "$_\n";
                $flag = 1;
            }

        }

    }
    close INPUT;
    close OUTPUT;

    print
"All template chains $chains from $pdbFL_ori are written into $pdbFL_new.\n";

    if ( $flag == 0 ) {

        die("chain $chains does not exist in $pdbFL_ori:$!");
    }

}

sub compareHomlogLsts_new {

    #!usr/bin/perl -w
    #Li Xue
    #1/15/2013
    #

#compare two homolog lists
#if they have the same pdbID and different chainID
#then this pdb complex is written to the output file
# perl ./compareHomlogLsts.pl HomologLst1 HomologLst2
#
#INPUT: two lists of homologs of two binding partners
#OUTPUT: a list of homo-complexes "../data/templates.lst"
#
#note: if one or two binding partners cannot find homologs, no output file will be generated.
#
#Usage example:
# perl ./compareHomlogLsts.pl ../data/1auiA_A_B/homologsOfquerySeq.lst ../data/1auiB_B_A/homologsOfquerySeq.lst
# perl ./compareHomlogLsts.pl ../data/1dkg_A_D/1dkgD_D_A/homologsOfquerySeq.lst ../data/1dkg_A_D/1dkgA_A_D/homologsOfquerySeq.lst
# perl ./compareHomlogLsts.pl ../data/trans208/2prg/2prg_B_C/2prgB_B_C/homologsOfquerySeq.lst ../data/trans208/2prg/2prg_B_C/2prgC_C_B/homologsOfquerySeq.lst
    use strict;
    use File::Basename;

    my $inputFL1  = shift @_;
    my $inputFL2  = shift @_;
    my $outputDIR = shift @_;

    my $query1 = basename( dirname($inputFL1) );
    my $query2 = basename( dirname($inputFL2) );

    my @Qry1homologs_in_homComplexes
      ;    #homologs of Qry chain 1 that forms homComplexes
    my @Qry2homologs_in_homComplexes;
    my @homComplexes;
    my $homComplexesFL = "$outputDIR/templates.lst";
    my @commonPDBIDs;

    #----------------

    print
"\n\nCompare two homolog lists ($inputFL1 and $inputFL2) and write homo-interolog list....\n\n";

    my $homologLst1 =
      &readPDBHomologLstFL($inputFL1);    #hash ref: $list1->{$pdbID} = @chainID
    my $homologLst2 = &readPDBHomologLstFL($inputFL2);    #%hash ref

    unlink($homComplexesFL) if ( -e $homComplexesFL );

    if ( $homologLst1 && $homologLst2 ) {

        #if both binding partners have homologs
        my @pdbIDlst1 = keys(%$homologLst1);

        print "\nPDB IDs in $inputFL1: @pdbIDlst1\n";

        my @pdbIDlst2 = keys(%$homologLst2);
        print "\nPDB IDs in $inputFL2: @pdbIDlst2\n";

        my %seen;
        if ( scalar @pdbIDlst1 >= scalar @pdbIDlst2 ) {
            @seen{@pdbIDlst1} = (1) x $#pdbIDlst1;
            @commonPDBIDs = grep { $seen{$_} } @pdbIDlst2;
        }
        else {
            @seen{@pdbIDlst2} = (1) x $#pdbIDlst2;
            @commonPDBIDs = grep { $seen{$_} } @pdbIDlst1;

        }

        if (@commonPDBIDs) {
            print "\nCommon PDB IDs are: @commonPDBIDs\n";

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
            @homComplexes = @{ &uniquePair( \@homComplexes ) };

            #2greB:2greO -> 2gre B:O
            @homComplexes = @{ &reformHomComplexPairs( \@homComplexes ) };

            &writeArrayIntoFL( \@homComplexes, $homComplexesFL );

            print
"compareHomlogLsts_new() is finished. $homComplexesFL is generated.  \n\n";

        }
        else {
            my $qryComplex = basename( basename( dirname($inputFL1) ) );
            print
"compareHomlogLsts_new() is finished. $qryComplex does not have homoComplexes.\n";

        }

    }
    else {

        #at least one of the binding partners cannot find homologs
        if ( !$homologLst1 ) {

            #		../data/1i4fA_A_C/homologsOfquerySeq.lst
            my ($partner) = $inputFL1 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
            print
"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";

        }
        if ( !$homologLst2 ) {

            #		../data/1i4fC_C_A/homologsOfquerySeq.lst
            my ($partner) = $inputFL2 =~ /\/(\w+)\/homologsOfquerySeq.lst$/;
            print
"compareHomlogLsts.pl is finished. $partner does not have homologs.\n";

        }
    }

    print "\n";

    return \@commonPDBIDs;

}

sub transformHomointerologFL2HomologFL {

    #Input (HomointerologFL):

#Homologs used in final prediction.
#	HomologID	logEval_H	PositiveS_H	Frac_Q	Frac_H	SimilarityScore
#	1tgsZ|Z:I	-108.578174314659	60.09	0.942857142857143	0.978165938864629	0.563818855368153
#	1hjaB|B:I	-189.505124806072	100	0.53469387755102	1	0.435053174254006
#	1choF|F:I	-189.505124806072	100	0.53469387755102	1	0.424073590325056

    #Output(homologFL):
    #	1tgsZ
    #	1hjaB
    #	1choF

    my @homointerologFLs = @{ shift @_ };
    my $outputFL         = shift @_;

    my @homologs;

    foreach my $homointerologFL (@homointerologFLs) {

        if ( !-e $homointerologFL ) {
            print
"Warning: $homointerologFL does not exist. No template file can be generated.\n";
            next;

        }

        open( INPUT, "<$homointerologFL" )
          || die("Cannot open $homointerologFL:$!");
        while (<INPUT>) {

            s/[\n\r]//mg;

            if (/^(\w{5})\|\w+:\w+/) {

#	1tgsZ|Z:I	-108.578174314659	60.09	0.942857142857143	0.978165938864629	0.563818855368153
                push @homologs, $1;
            }

        }
        close INPUT;

        @homologs = &unique(@homologs);

        #--write to output file
        unlink $outputFL if ( -e $outputFL );
        open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");

        foreach (@homologs) {
            print OUTPUT "$_\n";

        }

    }
    close OUTPUT;

    if ( -e $outputFL ) {

        print "$outputFL generated.\n";
    }
    else {
        print "No template used by PS-HomPPI. $outputFL cannot be generated.\n";
    }

}

sub readRecLigLst {

    #Input file:
    #rec:lig pairs
    #A:B
    #C:B

    my $recLigLst = shift @_;
    my %receptors;
    my %ligands;

    open( INPUT, "<$recLigLst" )
      || die("Cannot open rec:lig file $recLigLst:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^[\s\t]*(\w+:\w+)/) {

            my ( $rec, $lig ) = split( /:/, $1 );
            $receptors{$rec} = 1;
            $ligands{$lig}   = 1;
        }

    }
    close(INPUT);

    my @rec = keys %receptors;
    my @lig = keys %ligands;

    if ( !@rec || !@lig ) {

        die("No rec or lig read from $recLigLst:$!");

    }

    return ( \@rec, \@lig );

}

sub getCaseIDs_new {

    #		1a2k_AB_C
    #		1a2y_AB_C
    #		1akj_AB_DE
    #		1avw_A_B
    #		1bth_LH_P

    use strict;

    my $caseIDFL = shift @_;
    my @caseIDs;

    open( INPUT, "<$caseIDFL" ) || die("Cannot open $caseIDFL:$!");

    foreach (<INPUT>) {
        s/[\n\r]//mg;

        if (/^(\w+)/) {

            #		1a2k_AB_C
            #		1a2y_AB_C
            #		1akj_AB_DE
            #		1avw_A_B
            #		1bth_LH_P
            my $caseID = $1;
            push @caseIDs, $caseID;

        }
    }

    if ( !@caseIDs ) {
        die("No caseIDs read from $caseIDFL:$!");
    }

    my $num = scalar @caseIDs;

    print "Totally $num cases read from $caseIDFL\n";

    return \@caseIDs;

}

#sub readStatisticsFL_all {
#
#    #read statisticsFL
#
#    #A:B and A':B' are homo-interolog
#
##Input $statisticsFL example:
##>QUERY PDBID + CHAINID|CHAINID:INTERACTING CHAINID
##HOMOLOG-PDBID+CHAINID|CHAINID:INTERACTING CHAINID	Num_Int	num_residue1	num_residue2	num_residue1_P	Bit_score	Eval_H	Eval_HP	PositiveS_H	PositiveS_HP	IdentityS_H	IdentityS_HP	alignLen_H	alignLen_HP	Frac_Q	Frac_H	Frac_P	Frac_HP	TP	TN	FP	FN	CC	Specificity	Sensitivity	Accuracy
##	>1jpsL|1jpsL:1jpsT	2	213	213	200
##	1uj3A|A:C	4	215	214	205	420	1e-118	3e-112	99	97	96	97	213	205	1	0.990697674418605	0.99	1	2	209	2	0	0.703747585080836	0.5	1	0.990610328638498
##	1ahwD|D:F	4	214	214	219	315	9e-087	8e-114	84	97	70	97	213	207	1	0.995327102803738	1	0.945205479452055	1	208	3	1	0.345139365888927	0.25	0.5	0.981220657276995
##	1jrhL|L:I	3	213	167	108	284	2e-077	0.059	78	48	66	30	212	89	0.995305164319249	0.990610328638498	0.42	0.805555555555556	1	163	2	0	0.573840508584859	0.333333333333333	1	0.987951807228916
#
#    use strict;
#
#    my $statisticsFL = shift @_;
#
#    my $queryID;    #	$queryID = '1avxA|A:B'
#    my $Num_Int;
#    my $Length1;    #length of A or A'
#    my $Length2;    #length of B or B'
#    my $logEval_H;
#    my $logEval_HP;
#    my $PositiveS_H;
#    my $PositiveS_HP;
#    my $IdentityS_H;
#    my $IdentityS_HP;
#    my $alignLen_H;
#    my $alignLen_HP;
#    my $Frac_Q;     #(endQ-startQ)/lenQ
#    my $Frac_H;     #(endH-startH)/lenH
#    my $Frac_P;     #(endP-startP)/lenP
#    my $Frac_HP;    #(endHP-startHP)/lenHP
#    my $homologID;
#    my $prevID;
#    my @temp;
#
#    my $flag          = 0;
#    my $flag_homologs = 1
#      ; #0: this query has at least one homolog in $statisticsFL. 1: this query has no homolog in $statisticsFL.
#
#    open( INPUT, "<$statisticsFL" ) || die("Cannot open $statisticsFL:$!");
#    foreach (<INPUT>) {
#        s/[\n\r]//mg;
#        if (/^>((\w+)\|\w+:\w+)\t/) {
#            $flag    = 1;
#            $queryID = $1;    #'A|A:B'
#
#            (
#                my $tmp, $Num_Int->{$queryID}, $Length1->{$2}, my $tmpQ,
#                $Length2->{$2}
#            ) = split( /\t/, $_ );
#            next;
#        }
#
#        if ( $flag == 1 ) {
#            if (/^((\w{5})\|\w{1}:\w{1})/) {
#
##2gvdA|A:C	3	225	190	394	77.8	1e-17	0.0	52	94	29	94	159	357	0.892045454545455	0.675555555555556	1	0.906091370558376
#                @temp = split( /\t/, $_ );
#                $homologID = $1;    #1fs1D|D:C
#                my ( $homolog, $chn1, $chn2 ) =
#                  split( /[:\|]/, $homologID );    #$homolog=1fs1D
#                my $pdbID = substr( $homolog, 0, 4 );
#
#                $Num_Int->{$homologID}    = $temp[1];
#                $Length1->{"$pdbID$chn1"} = $temp[2];
#                $Length2->{"$pdbID$chn1"} = $temp[4];
#
#                $Length1->{"$pdbID$chn2"} = $temp[4];
#
#                $flag_homologs =
#                  0;   #0: this query has at least one homolog in $statisticsFL.
#
#                $prevID = $homologID;
#
#            }
#            elsif (/^[\t\s]+\d+/) {
#
##      377                     5e-38           49              33              337
#                @temp = split( /\t/, $_ );
#                $homologID = $prevID;    #1fs1D|D:C
#
#                next;    #skip the 2nd alignment
#            }
#
#            #			if($temp[5] eq '' || $temp[6] eq ''){
#            #				print "$homologID\n";#xue
#            #			}
#
#            if ( $temp[ 5 + 1 ] ne '' ) {
#                push @{ $logEval_H->{$queryID}->{$homologID} },
#                  &LogEval( $temp[ 5 + 1 ] );
#            }
#
#            if ( $temp[ 6 + 1 ] ne '' ) {
#
#                push @{ $logEval_HP->{$homologID} }, &LogEval( $temp[ 6 + 1 ] );
#            }
#
#            if ( !defined $homologID ) {
#
#                die(
#"homologID not defined. Current line: $_\nCheck $statisticsFL:$!"
#                );
#            }
#            push @{ $PositiveS_H->{$queryID}->{$homologID} },  $temp[ 7 + 1 ];
#            push @{ $PositiveS_HP->{$queryID}->{$homologID} }, $temp[ 8 + 1 ];
#            push @{ $IdentityS_H->{$queryID}->{$homologID} },  $temp[ 9 + 1 ];
#            push @{ $IdentityS_HP->{$queryID}->{$homologID} }, $temp[ 10 + 1 ];
#            push @{ $alignLen_H->{$queryID}->{$homologID} },   $temp[ 11 + 1 ];
#            push @{ $alignLen_HP->{$queryID}->{$homologID} },  $temp[ 12 + 1 ];
#            push @{ $Frac_Q->{$queryID}->{$homologID} },       $temp[ 13 + 1 ];
#            push @{ $Frac_H->{$queryID}->{$homologID} },       $temp[ 14 + 1 ];
#            push @{ $Frac_P->{$queryID}->{$homologID} },       $temp[16];
#            push @{ $Frac_HP->{$queryID}->{$homologID} },      $temp[17];
#
#            $prevID = $homologID;
#            next;
#        }
#
#    }
#
#    if ( !defined $queryID || !defined $PositiveS_H ) {
#        die("flag=$flag. Nothing read from $statisticsFL:$!");
#    }
#
#    return (
#
#        $Num_Int,    $Length1,     $Length2,      $logEval_H,
#        $logEval_HP, $PositiveS_H, $PositiveS_HP, $IdentityS_H,
#        $IdentityS_HP, $alignLen_H, $alignLen_HP, $Frac_Q, $Frac_H, $Frac_P,
#        $Frac_HP,
#
#    );
#
#}
#
sub addTail {

#-- pdb files generated by TMalign lost the information after the cooridinates
#-- here this information is extracted from the original pdb file and added to the TMalign output pdb file

    my ( $pdbFL_sup, $pdbFL_ori ) = @_;
    print
"\nAdding tail info (such as occupancy and B-factor from $pdbFL_ori) to $pdbFL_sup... \n";

    my $command =
"perl -ne 'if (/^ATOM/) {\$a = substr(\$_,54,length(\$_)-54);print \$a;}' $pdbFL_ori > tailInfo";

    system($command) == 0 or die("Cannot extract info from $pdbFL_ori:$!");

    #-- check row numbers before vertically concatenate them

    if ( !&equalRowNum( $pdbFL_sup, 'tailInfo' ) ) {
        die(
"$pdbFL_sup  has different row numbers as tailInfo. TailInfo file is from $pdbFL_ori :$!"
        );
    }

    #-- vertically concatenate two files
    system("paste -d '' $pdbFL_sup tailInfo > $pdbFL_sup.tmp") == 0
      or die("Cannot vertically concatenate $pdbFL_sup and $pdbFL_ori:$!");
    unlink $pdbFL_sup;

    copy( "$pdbFL_sup.tmp", $pdbFL_sup )
      or die("Cannot rename $pdbFL_sup.tmp to $pdbFL_sup:$!");
    unlink "tailInfo";
    unlink "$pdbFL_sup.tmp";
    return;

}    ## --- end sub addTail

sub addHETATM {
    my ( $pdbFL_sup, $pdbFL_ori ) = @_;
    print "\nAdd HETATM to $pdbFL_sup\n";

    system("sed -n '/^HETATM/p'  $pdbFL_ori >> $pdbFL_sup") == 0
      or die("Cannot add HETATM section. Check $pdbFL_sup and $pdbFL_ori:$!");

    return;
}    ## --- end sub addHETAM

sub equalRowNum {
    my ( $file1, $file2 ) = @_;
    my $ans;

    my $num_row1 = `wc -l $file1`;
    my $num_row2 = `wc -l $file2`;
    ($num_row1) = $num_row1 =~ /^(\d+)/;
    ($num_row2) = $num_row2 =~ /^(\d+)/;

    if ( $num_row1 == $num_row2 ) {
        $ans = 1;
    }
    else {
        $ans = 0;
    }

    return $ans;
}    ## --- end sub equalRowNum

#sub formatAlignFL {
#
##format the alignment file generated by TMalign
##INPUT:
##
##(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
##---MELITILEKTVSPDRL-ELEAAQKFLERAAVENLPTFLVELSRVLANPGNSQVARVAAGLQIKNSLTSKDPDIKAQYQQRWLAIDAN-AR
##   ::::::::::::.::: .:::::::::::::::::::::::::::::...:::::::::::::::::::......:::::::::::: ::
##   MSTAEFAQLLENSILSPDQNIRLTSETQLKKLSNDNFLQFAGLSSQVLIDENTKLEGRILAALTLKNELVSKDSVKTQQFAQRWITQVSPEAK
##
##OUTPUT:
## 1 *
## 2 *
## 3 ***
## 4 ***
##
## Frist column: sequence residue number of query structual sequence
## *: not aligned
## **: aligned but >= 5 angstroms
## ***: aligned and < 5 angstroms
##
#    my $alignFL = shift @_;
#    my $outFL   = shift @_;
#
#    my $flag = 1;
#
#    #-- read alignment
#    my $qrySeq;
#    my $align_mode;
#    open( INPUT, "<$alignFL" ) or die("Cannot open $alignFL:$!");
#    while (<INPUT>) {
#        s/[\n\r]//mg;
#
#        if (/^[\-A-Za-z]+$/) {
#            $qrySeq = $_;
#            next;
#        }
#        if (/^[\s:\.]+$/) {
#
#            #   :: :::::::::: . ::
#            $align_mode = $_;
#            $flag       = 0;
#            last;
#        }
#    }
#
#    close INPUT;
#
#    if ( $flag == 1 || !defined $qrySeq ) {
#        die("Nothing read from $alignFL:$!");
#    }
#
#    #--write output file
#    unlink $outFL if ( -e $outFL );
#    open( OUTPUT, ">>$outFL" ) or die("Cannot open $outFL:$!");
#
#    print OUTPUT "# generated by formatAlignFL{}\n";
#    print OUTPUT "# FORMAT:\n";
#    print OUTPUT "# Frist column: query residue\n";
#    print OUTPUT
#      "# Second column: sequence residue number of query structual sequence\n";
#    print OUTPUT "# *: not aligned with the template\n";
#    print OUTPUT "# **: aligned but >= 5 angstroms\n";
#    print OUTPUT "# ***: aligned and < 5 angstroms\n";
#
#    my $resiNum_qry = 0;
#    my @residues    = split( //, $qrySeq );
#    my @modes       = &translateAlignMode($align_mode);
#
#    for ( my $i = 0 ; $i < scalar @residues ; $i++ ) {
#
#        if ( $residues[$i] eq '-' ) {
#            next;
#        }
#        $resiNum_qry = $resiNum_qry + 1;
#        print OUTPUT "$residues[$i]\t$resiNum_qry\t$modes[$i]\n";
#    }
#    close OUTPUT;
#
#    print "\nReformat alignment file done: $outFL generated.\n";
#
#}
sub formatAlignFL {

#format the alignment file generated by TMalign
#INPUT:
#
#(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
#---MELITILEKTVSPDRL-ELEAAQKFLERAAVENLPTFLVELSRVLANPGNSQVARVAAGLQIKNSLTSKDPDIKAQYQQRWLAIDAN-AR
#   ::::::::::::.::: .:::::::::::::::::::::::::::::...:::::::::::::::::::......:::::::::::: ::
#   MSTAEFAQLLENSILSPDQNIRLTSETQLKKLSNDNFLQFAGLSSQVLIDENTKLEGRILAALTLKNELVSKDSVKTQQFAQRWITQVSPEAK
#
# OUTPUT:
#
# Qry_aa Qry_seqResNum Template_aa Template_seqResNum alignment_symbol
# T 1 T 9 ***
# E 2 T 10 ***
#
# Frist column: sequence residue number of query structual sequence
# *: not aligned
# **: aligned but >= 5 angstroms
# ***: aligned and < 5 angstroms
#
    my $alignFL = shift @_;
    my $outFL   = shift @_;

    #-- read alignment
    my @seqs;
    my $align_mode;
    open( INPUT, "<$alignFL" ) or die("Cannot open $alignFL:$!");
    while (<INPUT>) {
        s/[\n\r]//mg;

        if (/^[\-A-Za-z]+$/) {
            push @seqs, $_;
            next;
        }
        if (/^[\s:\.]+$/) {

            #   :: :::::::::: . ::
            $align_mode = $_;
            next;
        }
    }

    close INPUT;
    my $qrySeq = $seqs[0];
    my $temSeq = $seqs[1];

    #--write output file
    unlink $outFL if ( -e $outFL );
    open( OUTPUT, ">>$outFL" ) or die("Cannot open $outFL:$!");

    print OUTPUT "# generated by formatAlignFL{}\n";
    print OUTPUT "# FORMAT:\n";
    print OUTPUT "# *: not aligned with the template\n";
    print OUTPUT "# **: aligned but >= 5 angstroms\n";
    print OUTPUT "# ***: aligned and < 5 angstroms\n";
    print OUTPUT "#Qry_aa Qry_seqResNum Template_aa Template_seqResNum alignment_symbol\n";

    my $resiNum_Q = 0;
    my $resiNum_T = 0;
    my @residues_Q    = split( //, $qrySeq );
    my @residues_T    = split( //, $temSeq );
    my @modes       = &translateAlignMode($align_mode);
    for ( my $i = 0 ; $i < scalar @residues_Q ; $i++ ) {


        if ($residues_Q[$i] ne '-'){
          $resiNum_Q = $resiNum_Q + 1;
        }
        if ($residues_T[$i] ne '-'){
            $resiNum_T = $resiNum_T + 1;
        }

        if ( $residues_Q[$i] eq '-' || $residues_T[$i] eq '-' ) {
            next;
        }

        print OUTPUT "$residues_Q[$i]\t$resiNum_Q\t$residues_T[$i]\t$resiNum_T\t$modes[$i]\n";
    }
    close OUTPUT;

    print "\nReformat alignment file done. $outFL generated.\n";

}


sub translateAlignMode {

    #use *, **, *** to replace: space, : , .
    #
    #
    my $align_mode = shift @_;    #'   :::::...:::    ....::::'
    my @align_mode_new;

    foreach ( split( //, $align_mode ) ) {
        if (/\s+/) {
            push @align_mode_new, '*';
        }
        elsif (/:/) {
            push @align_mode_new, '***';
        }
        elsif (/\./) {
            push @align_mode_new, '**';
        }
        else {
            die("Unrecognized alignment symbol: $_ :$!");
        }

    }

    return @align_mode_new;

}

sub changeATOM2HETATM {
    use File::Copy;

    my @aa = qw(Ala
      Arg
      Asn
      Asp
      Cys
      Glu
      Gln
      Gly
      His
      Ile
      Leu
      Lys
      Met
      Phe
      Pro
      Ser
      Thr
      Trp
      Tyr
      Val
      CYM
      CSP
      CYF
      NEP
      ALY
      MLZ
      MLY
      M3L
      HYP
      PTR
      SEP
      TOP
      TYP
      TYS);
    @aa = map { uc($_) } @aa;

    my $pdbFL = shift @_;

    my @HETATMresi;
    my $flag      = 0;              #$flag =1 meaning there is cofactor in ATM
    my $pdbFL_tmp = "$pdbFL.tmp";
    unlink $pdbFL_tmp if ( -e $pdbFL_tmp );
    open( OUTPUT, ">>$pdbFL_tmp" ) or die("Cannot open $pdbFL_tmp:$!");

    open( INPUT, "<$pdbFL" ) || die("Cannot open $pdbFL:$!");
    while (<INPUT>) {
        s/[\n\r]//gm;

        if (/^ATOM.+/) {
            my $resi = substr( $_, 17, 3 );
            my $line = $_;

            if ( !( $resi ~~ @aa ) ) {
                $flag = 1;
                substr( $line, 0, 6, 'HETATM' );
                push @HETATMresi, $resi;
            }
            print OUTPUT "$line\n";
        }
    }
    close INPUT;
    close OUTPUT;

    unlink $pdbFL;
    move( $pdbFL_tmp, $pdbFL ) or die("Cannot rename $pdbFL_tmp to $pdbFL:$!");

#    if ( $flag == 1 ) {
#        @HETATMresi = &unique(@HETATMresi);
#        print
#"$pdbFL had non-amino acid (@HETATMresi)in ATOM section. They are moved to HETATM section.\n";
#    }
#    else {
#
#        print "$pdbFL has no cofactors. No atoms need to be moved to HETATM.\n";
#    }
#
    print "Check HETATM atoms done.\n";

    #    my $command = "sed -r -i 's/^ATOM(.+$resi.+\$)/HETATM\\1/g' $pdbFL";
    #    print "$command\n";
    #    system($command) == 0
    #      or die("Cannot change GDP in ATOM to HETATM from $pdbFL:$!");
}

sub cleanPDBFL {

    #1. add TER to the end of each chain
    #2. change non-aa (such as, GDP) in the ATOM section to HETATM
    #3. TODO: remove additional AltLoc

    #1.-add TER to the end of each chain (required by TMalign)
    my $pdbFL_ori = shift @_;
    my $pdbFL_new = shift @_;

    print "\n\nCleaning PDB file: $pdbFL_ori ..\n\n";
    unlink $pdbFL_new if ( -e $pdbFL_new );

    &addTER( $pdbFL_ori, $pdbFL_new );

#2.- change non-aa (such as, GDP) in ATOM section to HETATM in pdb file (tmalign will turn them into X in the alignment file. Also, haddock cannot process it, 1GRN, 1F6M)

    &changeATOM2HETATM($pdbFL_new);

#3. separate hetatm section from the original pdb (TMalign output lose the information of hetatm)

    system("egrep '^HETATM' $pdbFL_new > $pdbFL_new.HETATM");
    system("egrep '^ATOM' $pdbFL_new > $pdbFL_new.tmp ;mv $pdbFL_new.tmp $pdbFL_new") == 0
      or die("Failed.  NO ATOM found in $pdbFL_new:$!");

    print "\nCleaning PDB file done.\n\n";

}

sub preparePDBFLs4qry {

# When user chose not to upload pdb files, so we do not know the structure of query proteins.
# To cluster the templates, we need to generate superimposed model based on the query protein pdb
#
# Here, we  treat Template 1 as the structure of the query and superimpose it to other templates

    # 1. Change the chain IDs
    # 2. Clean the PDB file
    # 3. Output: $outputDIR/protein1.pdb;  $outputDIR/protein2.pdb
    #

    our $pdb_chainPY;
    my $rec =
      'A';    #these are the chain IDs that will be used in superimposed models
    my $lig =
      'B';    #these are the chain IDs that will be used in superimposed models

    my $rec_pdbFL_template1 = shift @_;
    my $lig_pdbFL_template1 = shift @_;
    my $outputDIR           = shift @_;

    my $unbound_r_PDBFL = "$outputDIR/protein1.pdb";
    my $unbound_l_PDBFL = "$outputDIR/protein2.pdb";

    copy( $rec_pdbFL_template1, $unbound_r_PDBFL )
      or die("Cannot copy $rec_pdbFL_template1 to $unbound_r_PDBFL:$!");
    copy( $lig_pdbFL_template1, $unbound_l_PDBFL )
      or die("Cannot copy $lig_pdbFL_template1 to $unbound_l_PDBFL:$!");

    &cleanPDBFL( $unbound_l_PDBFL, "$unbound_l_PDBFL.tmp" );
    &cleanPDBFL( $unbound_r_PDBFL, "$unbound_r_PDBFL.tmp" );
    move( "$unbound_r_PDBFL.tmp", $unbound_r_PDBFL ) or die("Failed:$!");
    move( "$unbound_l_PDBFL.tmp", $unbound_l_PDBFL );
    move( "$unbound_r_PDBFL.tmp.HETATM", "$unbound_r_PDBFL.HETATM" )
      if ( -e "$unbound_r_PDBFL.tmp.HETATM" );
    move( "$unbound_l_PDBFL.tmp.HETATM", "$unbound_l_PDBFL.HETATM" )
      if ( -e "$unbound_l_PDBFL.tmp.HETATM" );

    #--add chain ID
    my $currentDIR         = getcwd;
    my $dirname            = dirname($unbound_r_PDBFL);
    my $basename_u_r_PDBFL = basename($unbound_r_PDBFL);
    my $basename_u_l_PDBFL = basename($unbound_l_PDBFL);
    system(
"cd $dirname && python $pdb_chainPY -$rec $basename_u_r_PDBFL  > tmp && mv tmp $basename_u_r_PDBFL && cd $currentDIR"
      ) == 0
      or die("Cannot add chain ID to $unbound_r_PDBFL:$!");
    system(
"cd $dirname && python $pdb_chainPY -$lig $basename_u_l_PDBFL  > tmp && mv tmp $basename_u_l_PDBFL && cd $currentDIR"
      ) == 0
      or die("Cannot add chain ID to $unbound_l_PDBFL:$!");

    print "pdb files for qry prepared: $unbound_r_PDBFL, $unbound_l_PDBFL\n";
    return ( $unbound_r_PDBFL, $unbound_l_PDBFL );

}

sub preparePDBFLs4qry2 {

# To cluster the templates, we need to generate superimposed model based on the query protein pdb
#

    # 2. Clean the PDB file
    # 3. Output: $outputDIR/protein1.pdb;  $outputDIR/protein2.pdb
    #

    my $rec =
      'A';    #these are the chain IDs that will be used in superimposed models
    my $lig =
      'B';    #these are the chain IDs that will be used in superimposed models

    my $rec_pdbFL_template1 = shift @_;
    my $lig_pdbFL_template1 = shift @_;
    my $outputDIR           = shift @_;

    my $unbound_r_PDBFL = "$outputDIR/protein1.pdb";
    my $unbound_l_PDBFL = "$outputDIR/protein2.pdb";

    copy( $rec_pdbFL_template1, $unbound_r_PDBFL )
      or die("Cannot copy $rec_pdbFL_template1 to $unbound_r_PDBFL:$!");
    copy( $lig_pdbFL_template1, $unbound_l_PDBFL )
      or die("Cannot copy $lig_pdbFL_template1 to $unbound_l_PDBFL:$!");

    &cleanPDBFL( $unbound_l_PDBFL, "$unbound_l_PDBFL.tmp" );
    &cleanPDBFL( $unbound_r_PDBFL, "$unbound_r_PDBFL.tmp" );
    move( "$unbound_r_PDBFL.tmp", $unbound_r_PDBFL ) or die("Failed:$!");
    move( "$unbound_l_PDBFL.tmp", $unbound_l_PDBFL );
    move( "$unbound_r_PDBFL.tmp.HETATM", "$unbound_r_PDBFL.HETATM" )
      if ( -e "$unbound_r_PDBFL.tmp.HETATM" );
    move( "$unbound_l_PDBFL.tmp.HETATM", "$unbound_l_PDBFL.HETATM" )
      if ( -e "$unbound_l_PDBFL.tmp.HETATM" );

#    #--add chain ID
#    my $currentDIR         = getcwd;
#    my $dirname            = dirname($unbound_r_PDBFL);
#    my $basename_u_r_PDBFL = basename($unbound_r_PDBFL);
#    my $basename_u_l_PDBFL = basename($unbound_l_PDBFL);
#    system(
#"cd $dirname && python $pdb_chainPY -$rec $basename_u_r_PDBFL  > tmp && mv tmp $basename_u_r_PDBFL && cd $currentDIR"
#      ) == 0
#      or die("Cannot add chain ID to $unbound_r_PDBFL:$!");
#    system(
#"cd $dirname && python $pdb_chainPY -$lig $basename_u_l_PDBFL  > tmp && mv tmp $basename_u_l_PDBFL && cd $currentDIR"
#      ) == 0
#      or die("Cannot add chain ID to $unbound_l_PDBFL:$!");
#
    print "pdb files for qry prepared: $unbound_r_PDBFL, $unbound_l_PDBFL\n";
    return ( $unbound_r_PDBFL, $unbound_l_PDBFL );

}
