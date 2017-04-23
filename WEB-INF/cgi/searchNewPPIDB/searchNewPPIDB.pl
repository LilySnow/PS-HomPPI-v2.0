#!/usr/bin/perl -w
#
#Li Xue
#12/6/2010

#format of DistFL
# :A,:B,41:6.9002,,41:3.60848,,,,,,,
#  41 is the residue in chain B,  The residue in chain A is the first residue in the sequence file.
# ",," means that the distance is larger than 13 angstrom
#
# And if a chain has no int at all, this chain will be not appear in dist file.
#
#
#Check the new PPIDB. If the proteins in the input file exists in PPIDB, output their sequence and interface information. Input file is the homologs of the query seq, which is generated by BLASTp.
#Input: the homologs of query sequence, e.g, ./data/1lukA/PDBIDchainID_homologsOfquerySeq.lst
#Output: the sequence and interface info of the homologs, e.g., seq_interface_homologsOfquerySeq.lst
#
#
#   INPUT file:
#   2vl5, B
#   2vl5, D
#
# e.g.,
# perl searchNewPPIDB.pl ./example/test_NPS.lst conf/config
# perl searchNewPPIDB.pl ../data/nr2_test.txt conf/config
# perl searchNewPPIDB.pl ../data/all102853.lst

use strict;
use File::Basename;
use File::Path;

#---- new PPIDB data hierarchy
my $precomputedDIR =
'/data/isubk/einstein/home/ppidev/ppidev/ppiInfoNew/asymmetric/precomputed';

if ( !-d $precomputedDIR ) {
	die("$precomputedDIR does not exist:$!");
}

my $s2cDIR       = '/data/web_servers/lxue/S2C';

if (!-d $s2cDIR){
    die("$s2cDIR does not exist:$!");
}

#------


my $searchNewPPIDBDIR = dirname($0);
    #change this diretory according to where searchNewPPIDB is extracted.


  my $distThresh_Default = 4;    #4;
my $rasaThr_Default    = 5;    #5;

my $inputFL  = shift @ARGV;
my $configFL = shift @ARGV;
my $basename = basename( $inputFL, ( '.lst', '.txt' ) );

if ( !defined $inputFL ) {
	die("Please input an input file:$!");
}

my $sequencesDIR = $precomputedDIR . '/sequences';
my $rasaDIR      = $precomputedDIR . '/rasaBeforeComplexation';
my $distancesDIR;
my $distThresh;
my $rasaThr;         #[0-100]
my $outDIR;
my $interfaceDef;    # eg., $interfaceDef = atomDistances

#set parameters
if ( defined $configFL ) {

   #read config file to get interface definition parameters and output directory
	( $interfaceDef, $distThresh, $rasaThr, $outDIR ) =
	  &readConfigFL($configFL);

	if ( !defined $distThresh ) {
		print
"config file does not specify distance threshold. Now use the default value.\n";
		$distThresh = $distThresh_Default;
	}
	if ( !defined $rasaThr ) {
		print
"config file does not specify RASA threshold. Now use the default value.\n";
		$rasaThr = $rasaThr_Default;
	}
	if ( !defined $outDIR ) {
		$outDIR = "$searchNewPPIDBDIR/out";
	}
	if ( !defined $interfaceDef ) {
		$interfaceDef = 'atomDistances';
	}
}
else {

	#default parameter values

	$interfaceDef = 'atomDistances';
	$distThresh   = $distThresh_Default;
	$rasaThr      = $rasaThr_Default;
	$outDIR       = "$searchNewPPIDBDIR/out";

}

$distancesDIR = $precomputedDIR . '/distances' . "/$interfaceDef";

my $outputFL = "$outDIR/seq_int\_$basename.lst";

# process starts
my $subOutDIR;    #"$outDIR/2NZ8_A"
my $PDBID;
my $chainID;
my $seq;
my $interface;
my $ATOMindex
  ; #array ref. the index of aa in $seq of PPIDB. Some aa may be absent from PDB file.
my $totalNum = 0;   #The total number of protein chains in the input file
my $num1     = 0;   # The number of sequences in the output file
my $num2     = 0;   # The number of sequences that does not exist in NMR dataset
my $num_S2Cerr =0; #The number of queries that has abnormal s2c file.

open( INPUT, "<$inputFL" ) || die("Cannot open $inputFL:$!");

&header_outputFL( $outputFL, $interfaceDef, $distThresh, $rasaThr );

foreach (<INPUT>) {
	s/[\n\r]//mg;

	#read PDBID and chainID
	if (/^([\dA-Za-z]{4})[,\s]{0,}([A-Za-z\d]{1})/) {

		$totalNum++;

		$PDBID   = uc($1);
		$chainID = $2;

		print "Searching PPIDB for ......  $PDBID $chainID\n";

		#DIRs
		$PDBID = lc($PDBID);
		my $seqFL = $sequencesDIR . "/$PDBID.seq";

		if ( !defined $seqFL ) {

			#did not find the seq file for $pdbID
			$num2++;

			print "$PDBID $chainID does not exist in new PPIDB!\n";
			next;

		}

		if ( -e $seqFL ) {

			#extract sequence info from new PPIDB
			( $seq, $ATOMindex ) =
			  &extractSeqFromPPIDB( $sequencesDIR, $seqFL, $chainID );

              if (! defined $seq || ! defined $ATOMindex){
                  # chain ID does not have sequence info/AtomIndex info in the seqFL
                  next;
              }


			#read s2c files
			my $pdbID_tmp = lc($PDBID);
			my $s2cFL     = $s2cDIR . "/$pdbID_tmp.sc";

			if ( !-e $s2cFL ) {
				$num_S2Cerr ++;
				print
"$s2cFL does not exist. S2C DB needs to be updated. Skip this query.\n";
				next;
			}

			my ( $s2c_residues, $SEQRESresinum, $ATOMresinum ) =
			  &readS2Cfile( $PDBID, $chainID, $s2cDIR );

			  if(! @$s2c_residues){
			  print "$s2cFL cannot process this complex $PDBID. Next.";
			  $num_S2Cerr ++;
			  next;
			  }

			$num1++;

			my $fastaSeq = join( '', @$s2c_residues );

			#extract interface info from new PPIDB

			my @int =
			  &extractIntFromPPIDB( $distancesDIR, $rasaDIR, $rasaThr,
				$distThresh, $PDBID, $chainID );    #@int=(0,1,0,0,...);

            if (scalar @$s2c_residues < scalar @int){
                #s2c seq (ie., the whole seq) is shorter than the ppidb seq (ie, the structual seq) !!
                print "s2c seq (ie., the whole seq) is shorter than the ppidb seq (ie, the structual seq) !! Probably something wrong with this s2c entry, for example, 1qgc_4. Skip this query.\n";
				next;
            }

			if ( !@int ) {

#@int is empty. $chainID is not in $PDBID.dis file of PPIDB, because it has no int at all.
				$interface = '0' x length($fastaSeq);
			}
			else {

				#
				my %int_atomResiNum;
				@int_atomResiNum{@$ATOMindex} = @int;

				#
				my $newInt =
				  &mapInt2fasta( \%int_atomResiNum, $s2c_residues,
					$ATOMresinum );

				if ( !defined $newInt ) {

					die("Check $s2cFL and $seqFL:$!");
				}
				$interface = join( '', @$newInt );

			}

			#write %int into an output file

			&writeOutputFL( $interface, $fastaSeq, $PDBID, $chainID,
				$outputFL );

		}
		else {
			$num2++;

			print "$PDBID $chainID does not exist in new PPIDB!\n";
			next;
		}

	}

}
close(INPUT);

#
print
"\nSearching PPIDB is finished. The distance threshold: $distThresh, RASA threshold:	$rasaThr\%.\n";

print
"$outputFL is generated. Total $totalNum queries. Total $num1 in $outputFL, and $num2 queries do not exist in  PPIDB, and $num_S2Cerr queries have abnormal s2c files.\n\n\n";

#---------------------------

sub readConfigFL {

   #read config file to get interface definition parameters and output directory

	my $configFL = shift @_;
	my $interfaceDef;
	my $distThresh;
	my $outDIR;

	open( CONFIG, "<$configFL" ) || die("cannot open $configFL:$!");
	foreach (<CONFIG>) {
		if (/^#/) {
			next;
		}

		if (/interfaceDefinition[\s\t]{0,}=[\s\t]{0,}(\w+)/)

		  #interfaceDefinition = atomDistances
		{
			$interfaceDef = $1;
			next;
		}

		if (/distanceThreshold[\s\t]{0,}=[\s\t]{0,}([\d\.]+)/) {
			$distThresh = $1;
			next;
		}

		if (/surfaceResidueThreshold[\s\t]{0,}=[\s\t]{0,}([\d\.]+)/) {
			$rasaThr = $1;
			next;
		}
		if (/outDir=[\s\t]+(.+)/) {
			$outDIR = $1;
			next;
		}
	}
	close(CONFIG);
	return ( $interfaceDef, $distThresh, $rasaThr, $outDIR );
}

sub extractSeqFromPPIDB {

	#extract sequence info from new PPIDB
	my $sequencesDIR = shift @_;
	my $seqFL        = shift @_;
	my $chainID      = shift @_;

	my $seq;
	my $index
	  ;  #the index of aa in $seq of PPIDB. Some aa may be absent from PDB file.

	if ( -e $seqFL ) {
		open( SEQ, "<$seqFL" ) || die("Cannot open $seqFL:$!");
	}
	else {
		print
"$PDBID $chainID does not exist in new PPIDB!\n";
		return;
	}

	my $flag = 0;
	foreach my $line (<SEQ>) {
		$line =~ s/[\r\n]//mg;    #remove the new line symbol

		#		print "$line\n";#xue
		if ( $line =~ /\w{0,}:$chainID,([A-Za-z]+)/ ) {
			$seq  = $1;
			$flag = 1;
			next;
		}

		#		if ( $flag == 1 && $line =~ /^[\s\t\d,\-]+$/ ) {
		if ( $flag == 1 && $line =~ /^[\s\da-zA-Z\-]+,/ ) {
			$line =~ s/\s//g;
			@$index = split( /,/, $line );
			$flag = 0;    #reset $flag
			last;
		}
	}
	close(SEQ);

    if (!defined $seq){
        print ("\nWARNING: No sequence info extract for chain $chainID from $seqFL !!\n\n");
    }

    if (!defined $index){
        print ("\nWARNING: No index extract for chain $chainID from $seqFL !!\n\n");
    }

	return ( $seq, $index );
}

sub extractIntFromPPIDB {

	#read distance file and rasa file and extract interface info from new PPIDB

	my $distancesDIR = shift @_;
	my $rasaDIR      = shift @_;
	my $rasaThr      = shift @_;
	my $distThresh   = shift @_;
	my $PDBID        = shift @_;
	my $chainID      = shift @_;

	my @int;

	#read distance file and extract interface info from new PPIDB

	my @int_dist = &readDistFL( $distancesDIR, $distThresh, $PDBID, $chainID );

	if ( !@int_dist ) {

		#@int_dist is empty. $chainID is not in $PDBID.dis file of PPIDB.
		#print "$chainID does not exist in dist file of PPIDB of $PDBID.\n";
		return @int;
	}

	#read RASA files to determine Surface residues
	my @surfaces =
	  &extractSurfaceFromPPIDB( $rasaDIR, $rasaThr, $PDBID, $chainID );

	#get interfaces on the surfaces

	for ( my $i = 0 ; $i < scalar @surfaces ; $i++ ) {

		#						print "i: $i\n";
		#
		#						print "int_dist: $int_dist[$i]\n";
		#						print "surfaces: $surfaces[$i]\n";#xue
		#						print "---\n";#xue

		push @int, $int_dist[$i] & $surfaces[$i];
	}

	return @int;

}

sub extractSurfaceFromPPIDB {

	#read RASA files to determine Surface residues

	my $rasaDIR = shift @_;
	my $rasaThr = shift @_;
	my $PDBID   = lc( shift @_ );

	#	my $PDBID   = shift @_;
	my $chainID = shift @_;
	my $rasaFL  = $rasaDIR . "/$PDBID" . '.rsa';
	my @surfaces;

	open( RASAFL, "<$rasaFL" ) || die("Cannot open $rasaFL");
	foreach (<RASAFL>) {

		s/[\n\r]//;

		if (/^\w{0,}:$chainID,([\d\.\-,]+)/) {
			my @rasaLine = split( /,/, $1 );

			#
			foreach my $rasa (@rasaLine) {
				if ( $rasa >= $rasaThr )    #it is a surface residue
				{
					push @surfaces, 1;
				}
				else {
					push @surfaces, 0;
				}
			}
			last;

		}
	}

	return @surfaces;

}

sub readDistFL {

	#read distance file and extract interface info from new PPIDB

    #format of DistFL
    # :A,:B,41:6.9002,,41:3.60848,,,,,,,
    #  41 is the residue in chain B,  The residue in chain A is the first residue in the sequence file.
    # ",," means that the distance is larger than 13 angstrom

	my $distancesDIR = shift @_;
	my $distThresh   = shift @_;
	my $PDBID        = lc( shift @_ );
	my $chainID      = shift @_;
	my @temp;

	my $distFL = $distancesDIR . "/$PDBID" . '.dis';
	my @int_dist;    #interface determined by distance

	open( DISTFL, "<$distFL" ) || die("Cannot open $distFL:$!\n");
	foreach (<DISTFL>) {

		#		print "$_";
		s/[\n\r]//mg;

		if (/^\w{0,}:$chainID,\w{0,}:\w{1},([:\d\.\-,]+)$/) {
			my $distLine = $1;

			if ( substr( $distLine, length($distLine) - 1, 1 ) eq ',' ) {

            #if the last symbol of $1 is ',', for example, :A,:B,41:6.9002,,41:3.60848,
            #after 41:3.60848 there is another distance that is very large
            #
				$distLine = $distLine . ',9999';
				@temp = split( /,/, $distLine );
				pop @temp;    #get rid of the last element '9999'

			}
			else {

				@temp = split( /,/, $distLine );
			}
#            my $num_aa = scalar @temp;
#            print "The query protein $PDBID$chainID has $num_aa in the pdb file.\n";

			for ( my $i = 0 ; $i < scalar @temp ; $i++ ) {
                # the i-th residue of query chain

				# $temp[$i] = '137:3.25';

				my ($distance) = $temp[$i] =~ /:([\.\d]+)/;

				if ( !defined $distance ) {

					#the distance is very large

					if ( $int_dist[$i] ) {
						$int_dist[$i] =
						  $int_dist[$i] | 0
						  ; #get the union of interfaces if qry chain interacts with multiple chains
					}
					else {
						$int_dist[$i] = 0;    # 	$i is not an interface residue
					}

					next;
				}

				if ( $distance <= $distThresh ) {
					$int_dist[$i] = 1;        # 	$i is an interface residue
				}
				else {

					if ( $int_dist[$i] ) {
						$int_dist[$i] =
						  $int_dist[$i] | 0
						  ; #get the union of interfaces if qry chain interacts with multiple chains
					}
					else {
						$int_dist[$i] = 0;    # 	$i is not an interface residue
					}

				}
			}
		}

	}

	close(DISTFL);

	return (@int_dist);

}

#sub writeOutputFL {
#	my @ATOMindex = @{ shift @_ };
#	my @seq       = split( //, shift @_ );
#	my @int       = @{ shift @_ };
#	my $outputFL  = shift @_;
#
#	unlink($outputFL) if ( -e $outputFL );
#	open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL:$!");
#	for ( my $i = 0 ; $i < scalar @ATOMindex ; $i++ ) {
#
#		print OUTPUT "$ATOMindex[$i]\t$seq[$i],$int[$i]\n";
#	}
#	close OUTPUT;
#
#	#	print "$outputFL is generated.\n";
#
#}

sub writeOutputFL {

	#	(\%int,$PDBID,$chainID);
	#read the complete fasta file into a hash and compare with %int
	#for aa that are in %int, write the outputfile accordingly;
	#for aa that are not in %int, use "?" to denote the interface.

	my $interface = shift @_;
	my $seq       = shift @_;
	my $PDBID     = shift @_;
	my $chainID   = shift @_;
	my $outputFL  = shift @_;

	open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL!\n");
	print OUTPUT ">$PDBID$chainID\n";
	print OUTPUT "$seq\n";
	print OUTPUT "$interface\n";
	close(OUTPUT);

}

sub header_outputFL {
	my $outputFL = shift @_;
	my ( $interfaceDef, $distThresh, $rasaThr ) = @_;

	unlink($outputFL) if ( -e $outputFL );
	open( OUTPUT, ">>$outputFL" ) || die("Cannot open $outputFL!\n");
	print OUTPUT
"# =========================================================================================================\n";
	print OUTPUT "# auto generated file  by author lixue\@iastate.edu\n";
	print OUTPUT "# Format description:\n";
	print OUTPUT "# >Protein name(PDBid + chain ID)\n";
	print OUTPUT "# Amino acid sequence in new PPIDB.\n";
	print OUTPUT
"# (non)interface of the protein: 1s denote interface residues, 0s denote non-interface residue\n";
	print OUTPUT "# interface definition: $interfaceDef\n";
	print OUTPUT "# distThr: $distThresh\n";
	print OUTPUT "# rasaThr: $rasaThr\n";
	print OUTPUT
"# =========================================================================================================\n";
	close(OUTPUT);
}
sub readS2Cfile_old {

	#Note: Sometimes s2c cannot process a PDB file. For example, 1kld.
	#
	my $pdbID   = shift @_;
	my $chainID = shift @_;
	my $s2cDIR  = shift @_;
	my $s2cFL   = $s2cDIR . "/$pdbID.sc";
	my @residues;
	my @SEQRESresinum;
	my @ATOMresinum;

	open( S2C, "<$s2cFL" ) || die("Cannot open $s2cFL:$!");
	foreach (<S2C>) {
		s/[\n\r]//gm;

		if (/^SEQCRD $chainID/) {

			#			$chainID=substr($_,7,1);
			push @residues,      substr( $_, 9,  1 );
			push @SEQRESresinum, substr( $_, 19, 5 );

			my ($temp) = substr( $_, 25, 6 ) =~ /([^\s\t]+)/;
			push @ATOMresinum, $temp;

		}
	}
	close S2C;



	if(!@residues){
		print("s2c file error. Sometimes s2c cannot process a PDB file. For example, 1kld. Check $s2cFL.\n");
	}

	return ( \@residues, \@SEQRESresinum, \@ATOMresinum );
}


sub readS2Cfile {

	#Note: Sometimes s2c cannot process a PDB file. For example, 1kld.
	#
    # A s2c file:
    # SEQCRD    A R ARG ARG   122    201 C C 0
    # SEQCRD    A K LYS ---   123      - - - 0
    #
	my $pdbID   = shift @_;
	my $chainID = shift @_;
	my $s2cDIR  = shift @_;
	my $s2cFL   = $s2cDIR . "/$pdbID.sc";
	my @residues;
	my @SEQRESresinum;
	my @ATOMresinum;

	open( S2C, "<$s2cFL" ) || die("Cannot open $s2cFL:$!");
	while (<S2C>) {
		s/[\n\r]//gm;

		if (/^SEQCRD\s+$chainID/) {

			#			$chainID=substr($_,7,1);

            my @tmp = split(/\s+/, $_);
            my $aa = $tmp[2];
            my $seqResNum = $tmp[5];
            my $atomResNum = $tmp[6];
			push @residues,      $aa;
			push @SEQRESresinum, $seqResNum;
			push @ATOMresinum, $atomResNum ;

		}
	}
	close S2C;



	if(!@residues){
		print("s2c file error. Sometimes s2c cannot process a PDB file. For example, 1kld. Check $s2cFL.\n");
	}

	return ( \@residues, \@SEQRESresinum, \@ATOMresinum );
}

sub mapInt2fasta {
	my $int_atomResiNum_ref = shift @_;
	my $residues_ref        = shift @_;    # in the order of fasta.
	my $ATOMresinum_ref     = shift @_;
	my $int;

	for ( my $i = 0 ; $i < scalar @$residues_ref ; $i++ ) {

		# the residue is missing from the structure
		if ( $ATOMresinum_ref->[$i] eq '-' ) {
			$int->[$i] = '?';
			next;
		}

		#no interface info for this residue if it is not among 20 types aa
		elsif ( !exists $int_atomResiNum_ref->{ $ATOMresinum_ref->[$i] } ) {
			$int->[$i] = '?';
			next;
		}

		#there is interface info for this residue
		else {
			$int->[$i] = $int_atomResiNum_ref->{ $ATOMresinum_ref->[$i] };
		}

	}
	return $int;
}
