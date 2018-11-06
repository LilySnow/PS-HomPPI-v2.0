#!/usr/bin/perl -w
#Li Xue
#July 3rd, 2014
#
#Combine multiple precalculated distance files into one haddock CA-CA restraint file.
#Only the residue pairs with at least one CA-CA distances within the specified threshold are passed onto the final file.
#
#INPUT: the $dir which has multiple precalculatd restraint files. Each file has the folllowing format:
# 48 A 91 B 6.247917
# 49 A 90 B 6.427874
# 49 A 91 B 4.569639
# 50 A 90 B 5.799189
#
#OUTPUT 1:
# assign (name ca and segid A and resi 62) (name ca and segid B and resi 140) 8.836 0.986 0.276
# assign (name ca and segid A and resi 36) (name ca and segid B and resi 262) 8.719 0.346 0.246
# assign (name ca and segid A and resi 36) (name ca and segid B and resi 143) 5.805 0.113 0.148
#
#OUTPUT 2:
#A aa_A B aa_B mean min max
#A 128 P 90 14.437 13.6659 15.1642
#A 128 P 92 12.270 12.163 12.4211
#A 61 P 95 13.603 13.5638 13.6424
#
#OUTPUT 3 (pml for visualization).


use strict;
use List::Util qw(min max sum);

our $CaDistThr = shift @ARGV; #-- 8.5
my $dir = shift @ARGV;
    #$dir has one cluster of template-based restraint files. Each file has the format in the perl header.
my $outputFL = shift @ARGV;
&genCaDisRestraintFLs_oneCluster($dir,$outputFL);


sub genCaDisRestraintFLs_oneCluster {
    #$dir has one cluster of template-based restraint files. Each file has the folllowing format:
    #48 A 91 B 6.247917
    #49 A 90 B 6.427874
    #49 A 91 B 4.569639
    #50 A 90 B 5.799189

    our $CaDistThr;

    my $dir      = shift @_; #input dir
    my $outputFL = shift @_; #in haddock tbl format
    my $outputFL2 = "$dir/final_caca.txt"; # easy format


  #-- read Ca-Ca distances that are mapped from templates to the query sequences
    opendir( DIR, $dir ) or die("Cannot open dir:$dir:$!");
    my @distFLs = grep { /^.+.contact$/ } readdir(DIR);
    closedir(DIR);

    my $num_distFLs = scalar @distFLs;
    print "There are $num_distFLs restraint files under inputDIR\n";

    if ( !@distFLs ) {

        print
"No distanceMapped files exist. No interface residue pairs mapped to this query chain pair.\n\n";
        return;
    }

    my $dist_map;
    #- $dist_map->{Ca_of_residue_in_protA}->{Ca_of_residue_in_protB}=(6.1,10.2)

    foreach my $distFL (@distFLs) {
        $distFL = "$dir/$distFL";

        open( INPUT, "<$distFL" ) or die("Cannot open $distFL:$!");
        while (<INPUT>) {
            s/[\n\r]//gm;
            if (/^-{0,1}\d+/) {
                #  48 A 91 B 6.247917
                my ( $resiNum_A, $chnA, $resiNum_B,$chnB, $dist ) = split( /\s+/, $_ );

                if ($dist =~/^\s{0,}$/){
                    die("Input file format wrong. Check $distFL:$!");
                }

                my $resi_A="$resiNum_A:$chnA";
                my $resi_B="$resiNum_B:$chnB";
                push @{ $dist_map->{$resi_A}->{$resi_B} }, $dist;
            }
        }
        close INPUT;

    }

    #-- calcualte the min, max and mean Ca-Ca distantce for each residue pair
    &header1($outputFL);
    &header2($outputFL2);
    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    open( OUTPUT2, ">>$outputFL2" ) or die("Cannot open $outputFL2:$!");
#    print OUTPUT2 "#A aa_A B aa_B mean min max\n";
    print OUTPUT2 "ChnID1 aa1 ChnID2 aa2 mean min max\n";

    foreach my $resi_A ( keys %$dist_map ) {
        foreach my $resi_B ( keys %{ $dist_map->{$resi_A} } ) {
            my @distance = @{ $dist_map->{$resi_A}->{$resi_B} };

            #            print "$resi_A\t$resi_B\n";
            #            print "@distance\n";
            my ($resiNum_A, $chnA) = split(/:/, $resi_A);
            my ($resiNum_B, $chnB) = split(/:/, $resi_B);

            my $min  = sprintf("%.2f", min(@distance));
            my $max  = sprintf("%.2f", max(@distance));
            my $mean = sprintf( "%.2f", sum(@distance) / scalar @distance );

            if ( $min > $CaDistThr ) {
                next;
            }

            #--give 0.5 angstroms more flexibility to the restraints

            my $lower = $mean - $min + 0.5;
            my $upper = $max - $mean + 0.5;

                      printf OUTPUT
"assign (name ca and segid $chnA and resi %s) (name ca and segid $chnB and resi %s) %.3f %.3f %.3f\n",
              $resiNum_A, $resiNum_B, $mean, $lower, $upper;
              printf OUTPUT2 "$chnA $resiNum_A $chnB $resiNum_B $mean $min $max\n";

        }
    }
    close OUTPUT;
    close OUTPUT2;

    #--write $outputFL2 into pml file
    my $outputFL3=&writePML($outputFL2);
    print "\n$outputFL generated.\n";
    print "$outputFL2 generated.\n";
    print "$outputFL3 generated.\n";

}
sub   header1{
    #header for haddock restraint file
    our $CaDistThr;
    my $outputFL = shift @_;
    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    print OUTPUT "! generated by $0\n";
    print OUTPUT "! CA-CA cutoff = $CaDistThr\n";
    print OUTPUT "! 0.5 angstroms are added both sides of the distance restraints derived from templates\n";
    close OUTPUT;

}
sub   header2{
    #header for CA-CA file

    our $CaDistThr;
    my $outputFL = shift @_;
    unlink $outputFL if ( -e $outputFL );
    open( OUTPUT, ">>$outputFL" ) or die("Cannot open $outputFL:$!");
    print OUTPUT "# generated by $0\n";
    print OUTPUT "# CA-CA cutoff = $CaDistThr\n";
    close OUTPUT;

}
sub writePML{
    #Input:
    #
    #ChnID1 aa1 ChnID2 aa2 mean min max
    #A 128 B 24 9.01 8.45 9.93
    #A 40 B 41 8.01 7.49 8.58

    use File::Basename;

    my $inputFL = shift @_;

    my $dirname = dirname ($inputFL);
    my $filename = basename ($inputFL, (".txt", ".lst"));
    my $outputFL = "$dirname/$filename.pml";


    my $command = "echo \"#generated by $0\" > $outputFL";
    system($command ) ==0 or die ("FAILED: $command:$!");

    $command = "egrep '^\\w' $inputFL |awk '{print \"distance pair\" NR \", chain \" \$1 \" and name ca and resi \" \$2 \", chain \" \$3 \" and name ca and resi \" \$4 }' >> $outputFL";

#    print "\nWrite pml file:\n";
#    print "$command\n";

    system($command ) ==0 or die ("FAILED: $command:$!");

    return $outputFL;
}
