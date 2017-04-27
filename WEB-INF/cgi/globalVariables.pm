#===============================================================================
#
#         FILE: globalVariables.pm
#
#  DESCRIPTION: This file contains all the global variables for PS-HomPPI server v2.0
#
#       AUTHOR: Li Xue, L.Xue@uu.nl
#      CREATED: 11/30/2016 05:52:29 PM
#===============================================================================

use strict;
use warnings;

our $intDef = 'alphaCarbonDistances'
  ;    #atomDistances  vanDerWaalsDistances alphaCarbonDistances

our $SERVERurl = 'http://ailab1.ist.psu.edu/PSHOMPPIv2.0/';
our $SERVER  = 'smtp.psu.edu';    #'mailhub.iastate.edu';
our $FROM    = 'L.Xue@uu.nl';     #'lxue@ist.psu.edu';    #'lixue@iastate.edu';
our $version = 'v2.0';
our $serverName = "PSHOMPPI$version";
my $serverDIR = "/data/web_servers/lxue/$serverName";
our $HOMPPIcgibin = "$serverDIR/WEB-INF/cgi";

our $dataDIR = "$serverDIR/uploadData";        #store all the user-uploaded data for HOMPPI
our $logDIR = "$dataDIR/LOGs";
our $safe_filename_characters = "a-zA-Z0-9_\.-";
our $mailprog =  "/usr/sbin/sendmail";    # Location of sendmail; may vary on your system

our $intNum_cutoff = 4; #-- if a homolog has <= 4 interfacial residues, this homolog will not be used in predictions.

#------------------- on hanavar server --------------------------#
#----------------------------------------------------------------#

# our $searchNewPPIDBDIR = '/data/web_servers/lxue/searchNewPPIDB';
#
# #our $BlastDIR    = '/data/web_servers/lxue/ncbi-blast-2.2.29+';
# our $blastBinDIR = "$BlastDIR/bin";
#
# our $BLASTDB_nr = "$BlastDIR/data/nr/nr";
# our $BlastDataset =
#   "$BlastDIR/data/nr_pdbaa_s2c/nr_pdbaa_s2c"
#   ;    #non-redundant protein dataset from S2C fasta file
# our $searchNewPPIDB_PS_PL =
#   "$searchNewPPIDBDIR/searchNewPPIDB_parnterSpecific.pl";
#
# #-------
# #our $BlastDataset =
# #  "$BlastDIR/data/scop_v175_prot/scop_v175_prot"
# #  ;    # => searchNewPPIDB_parnterSpecific_AtomRes.pl
# #our $searchNewPPIDB_PS_PL =
# #  "$searchNewPPIDBDIR/searchNewPPIDB_parnterSpecific_AtomRes.pl";
#
# our $pdb_chain_taxonomyLst = '/data/web_servers/lxue/pdb_chain_taxonomy.lst';
#
# #-------
# our $BlastP   = "$BlastDIR/bin/psiblast";    #not used.
# our $scoreThr = 0.5;

#------------------- on alcazar ---------------------------------#
#----------------------------------------------------------------#


#our $perl5LibDIR =  '/home/lixue/perl5/lib/perl5';
#
#------------------ hhpred related parameters -------------------#
# See example configure file: /home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/configExample.txt
#
#our $PYTHON = '/usr/bin/python';
#our $callHHpredPY =
#'/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/versionToWorkWith22-12-2016/runquery_abs.py';
#our $flag_hhpredrun = 'TRUE'; #if false, abort on error. if true, ignore and attempt to continue
##our $hhpredConfigTemplateFL = '/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/configExample.txt';
#our $addssPL = '/home/mwalter/hh/scripts/addss.pl'; #location of ssscript file
#our $uniprotDB_dir='/home/mwalter/project/v0.0.02/projectfiles/databases/uniprotHMMs'; #location of UNIPROT DB folder
#our $pdb70HMM_dir = '/home/mwalter/project/v0.0.02/projectfiles/databases/pdb70HMMs'; #location of PDB70 DB folder
#our $uniprotDB_name = 'uniprot20_2016_02'; #name of UNIPROT DB file
#our $pdb70DB_name = 'pdb70'; #name of PDB70 DB file
#our $pickle_dir = '/home/mwalter/project/v0.0.02/projectfiles/pickle.pickle'; #location of startscript pickle file
#our $queryName = 'query.fa'; #for now unsafe to modify, need to also modify HHR1OUTPUTLOCATION
#our $numItr_uniprotDB = 1; #number of hhsearch iterations against the UNIPROT DB
#our $numItr_PDB70 = 1; #number of hhsearch iterations against the PDB70 DB
#our $HHR1OUTPUTLOCATION = "%s/%s.hhr"; #contains QUERY1NAME and jobFolder, DO NOT MODIFY/UNCOMMENT
#our $lookuptable_dir = '/home/mwalter/project/v0.0.02/projectfiles/databases/lookuptable'; # location of the folder with the lookup table created by startscript.py
#our $PPIDBseqFL = '/home/mwalter/project/v0.0.02/projectfiles/databases/PDBallFASTA/pdb_seqres.txt'; # location of PDB all sequences fasta file
#our $evalThr_uniprot =  0.0001; ##both the -e -E options in the HHsearch against UNIPROT, E-value cutoff
##HHRPDBMINSID = 5 #-qid, dont use for now
#our $miniCov_hhpredPDB70  =  5; # the -cov option in the HHsearch against pdb70, min coverage
#our  $maxSID_hhpredPDB70  =  100; #the -id option in the HHsearch against pdb70, max SID cutoff (for less redundancy)
#our $EvalThr_hhpredPDB70 = 1; #both the -e -E options in the HHsearch against pdb70, E-value cutoff
#our $PvalThr_hhpredPDB70 = 40; #min P-value cutoff for the HHsearch against PDB70
#our $minResNum_in_alignment = 1; #min number of residues in alignment. must be >0 or results include start = stop residues
#our $BlastDIR =
#  '/home/sbgrid/programs/x86_64-linux/blastplus/2.3.0';    # on alcazar
# our $Rscript = "/home/lixue/tools/Anaconda3/bin/Rscript";

#------------------- on honavar01 ---------------------------------#
#----------------------------------------------------------------#


our $perl5LibDIR =  '/home/lxue/perl5/lib/perl5';
our $profit = '/home/lxue/tools/profit';

#------------------ hhpred related parameters -------------------#
# See example configure file: /home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/configExample.txt
#
our $PYTHON = '/usr/bin/python';
our $callHHpredPY ='/data/web_servers/lxue/call_hhpred/MICKSCODE/runquery_abs_report_all.py';
our $flag_hhpredrun = 'TRUE'; #if false, abort on error. if true, ignore and attempt to continue
#our $hhpredConfigTemplateFL = '/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/configExample.txt';
our $addssPL = '/data/web_servers/lxue/call_hhpred/SOFTWARE/hhsuite/hhsuite-3.0.1-Source/scripts/addss.pl'; #location of ssscript file
our $uniprotDB_dir='/data/web_servers/lxue/call_hhpred/DATABASES/uniprotHMMs'; #location of UNIPROT DB folder
our $pdb70HMM_dir = '/data/web_servers/lxue/call_hhpred/DATABASES/pdb70HMMs'; #location of PDB70 DB folder
our $uniprotDB_name = 'uniprot20_2016_02'; #name of UNIPROT DB file
our $pdb70DB_name = 'pdb70'; #name of PDB70 DB file
#our $pickle_dir = '/home/mwalter/project/v0.0.02/projectfiles/pickle.pickle'; #location of startscript pickle file
our $queryName = 'query.fa'; #for now unsafe to modify, need to also modify HHR1OUTPUTLOCATION
our $numItr_uniprotDB = 1; #number of hhsearch iterations against the UNIPROT DB
our $numItr_PDB70 = 1; #number of hhsearch iterations against the PDB70 DB
our $HHR1OUTPUTLOCATION = "%s/%s.hhr"; #contains QUERY1NAME and jobFolder, DO NOT MODIFY/UNCOMMENT
our $lookuptable_dir = '/data/web_servers/lxue/call_hhpred/DATABASES/lookuptable'; # location of the folder with the lookup table created by startscript.py
our $PPIDBseqFL = '/data/web_servers/lxue/call_hhpred/DATABASES/PDBallFASTA/pdb_seqres.txt';# location of PDB all sequences fasta file
our $evalThr_uniprot =  0.0001; ##both the -e -E options in the HHsearch against UNIPROT, E-value cutoff
#HHRPDBMINSID = 5 #-qid, dont use for now
our $miniCov_hhpredPDB70  =  5; # the -cov option in the HHsearch against pdb70, min coverage
our  $maxSID_hhpredPDB70  =  100; #the -id option in the HHsearch against pdb70, max SID cutoff (for less redundancy)
our $EvalThr_hhpredPDB70 = 1; #both the -e -E options in the HHsearch against pdb70, E-value cutoff
our $PvalThr_hhpredPDB70 = 40; #min P-value cutoff for the HHsearch against PDB70
our $minResNum_in_alignment = 1; #min number of residues in alignment. must be >0 or results include start = stop residues
our $BlastDIR ='/data/web_servers/lxue/ncbi-blast-2.2.29+';
our $Rscript = "/home/lxue/tools/Anaconda3/bin/Rscript";

#----------------------------------------------------------------------------------
#

our $searchNewPPIDBDIR = 'searchNewPPIDB';

our $blastBinDIR = "$BlastDIR/bin";
our $BLASTDB_nr  = "$BlastDIR/data/nr/nr";
our $BlastDataset =
  "$BlastDIR/data/nr_pdbaa_s2c/nr_pdbaa_s2c"
  ;    #non-redundant protein dataset from S2C fasta file
our $searchNewPPIDB_PS_PL =
  "$searchNewPPIDBDIR/searchNewPPIDB_parnterSpecific.pl";


#-------
our $BlastP   = "$BlastDIR/bin/psiblast";    #not used.
our $scoreThr = 0.5;

#-- variable used for mapping interfacial pairs from A':B' to A:B

our $s2cDIR = '/data/web_servers/lxue/S2C';
our $precomputedDIR =  '/data/isubk/einstein/home/ppidev/ppidev/ppiInfoNew/asymmetric/precomputed';

#our $CaDistThr = 8.5; #15;                             #8.5
#; #Ca-Ca distance threshold used for mapping template interating residue distances to query


#-------------------------- for genSeqIntFL4unboundQry.pl

#--- on lrs-honavar01.ist.psu server
#
our $unboundPDB_DIR =
  '/data/isubk/einstein/home/ppidev/ppidev/ppiInfoNew/asymmetric/pdb'
  ;    #--to get pdb files for templates
our $pdb_chainPY      = '/home/lxue/tools/pdb-tools/pdb_chain.py'; #set chain ID for a pdb file
our $TMalign          = '/home/lxue/tools/TMtool/TMalign';
our $pdb2atomresnumPL = '/home/lxue/tools/pdb-tools/PDB2AtomResNum.pl';

#- on alcazar
#our $unboundPDB_DIR = '/data/lixue/DBs/ppidev/ppiInfoNew/asymmetric/pdb';    #--to get pdb files for templates
#our $pdb_chainPY      = '/home/lixue/tools/pdb-tools/pdb_chain.py'; #set chain ID for a pdb file
#our $TMalign          = '/home/lixue/tools/TMtool/TMalign';
#our $pdb2atomresnumPL = '/home/lixue/tools/PDB2AtomResNum.pl';

#- tools

our $toolDIR ='/home/lxue/tools';
our $pdbtoolDIR = '/home/lxue/tools/pdb-tools';
our $perlDIR = '/usr/bin/perl';
our $path = "/home/lxue/tools/Anaconda3/bin:/usr/bin:/usr/local/bin:/bin";
our $PSHomPPI_path="$toolDIR:$pdbtoolDIR:$perlDIR:$path";


1;
