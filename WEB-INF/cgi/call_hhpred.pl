use strict;
use warnings;

my $filename = 'settings.config';
open(my $fh, ">$filename") or die ("Could not open file $filename: $!");
print $fh "E-value=0.0001\n";
close $fh;

my $FastaFileName = 'query1.fasta';
open(my $fh, '>', $FastaFileName) or die "Could not open file '$FastaFileName' $!";
print $fh ">example fasta file";
close $fh;
print "done\n";



&callpsiBlast_old( $BlastDataset, $FastaFileName, $Qry1BlastFL, $EvalThr );    #Blast the query directly against PDB database

#system("command", "arg1", "arg2", "arg3");
#YOUR_FOLDER_PATH=‘../../../../../../../home/mwalter/custom/location/YOURFOLDERNAME’
#qsub -q short -F YOUR_FOLDER_PATH /home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/runquery.py



my $FOLDER_PATH = '../../../../../../../'.$CurrentFolder;
system("qsub", "-q short", "-F ".$FOLDER_PATH, "/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/runquery.py");
