
Before Run
=======
Open searchNewPPIDB.pl and change "my $searchNewPPIDBDIR=" according to where searchNewPPIDB is extracted.
For example, 
my $searchNewPPIDBDIR='D:/lixue/courses/!BCB569/final_project/HomPPI/searchNewPPIDB';

To Run
=======

Use the following command:

perl searchNewPPIDB.pl <inputFile> [configurationFile]

The inputFile is mandatory where as the configurationFile is optional and can be used to overwrite default arguments.
Appropriate Error message is flagged



OutPut File
=============
The output is written to a file named seq_interface_$inputFile$ where $inputFile$ is replaced by the inputFile.
You should have permissions to write to the  directory (appropriate error message will be flagged).
If a value does not exist for an input in PPIDB it is replaced by null



Structure of Input File
=======================
The input file  should be in the format of pdbID,chainID.
The name is a sequentially increasing value. If it is repeated only one of them will be considered.
Comma or a space can be used as separators. Lines starting with '#' will be considered as comments
1a81,A
1a81,E
1a81,G
1a81,I
1a81,K
1m61,A
1csyA
1k9aA
1k9aB
1a9a C




ConfigurationFile
================
The following are the default values used by the program. If you want to overwrite anything, set
the correct value in the configuration file

factoryURI=http://hebb.cs.iastate.edu:8080/wsrf/services/ppi/ChainInterfaceQueryFactoryService
interfaceResidueType = ASA
interfaceResidueThreshold = 1	            
similarityType = ignore
similarityThreshold = 10
surfaceResidueThreshold = 30

#presence of lowercase=true will force the input(protien ID) into lowercase before using it
lowercase=false





