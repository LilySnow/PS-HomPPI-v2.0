# Li Xue
# July 2nd, 2014
#
# Given some templates,
# 1) calculate atom-atom pair distance. One template one file.
# 2) calculate the CA-CA distance for residue-residue pair. One template one file.
# 3) map qry atomResNums to templates' CA-CA files.
# 4) Add alignment symbols to the CA-CA files from step 3.
# 5) combine the output from step 4) into one file.

set -e

#----------global var
#contact='/home/lixue/tools/contact-chainID_new' # the cpp to calculate contacts
contact='/home/lxue/tools/contact-chainID_new' # the cpp to calculate contacts


if [ ! -e $contact ];then
    echo "$contact does not exist"
fi

#---------input

DIR=$1 #-- dir for a cluster of pdb files
atomDistThr=$2; #20 #-- atom distance threshold for calculating the initial contact file, which is later used to extract CA-CA distances
caDistThr=$3; #15 #-- only include contacts with CA-CA distance <= 8.5 angstroms into the final output file

scriptDIR=`pwd`


num_pdbFL=0




#-------------------------------------------------------
#-- step 1 and 2: Calculate CA-CA contacts for each template pdb file

for i in `ls $DIR/template*.pdb`;do
    printf "Calculate CA-CA contact file from $i \n"
    (( num_pdbFL= num_pdbFL+1 ))
    filename=`basename $i ".pdb"`
    $contact $i $atomDistThr > $DIR/$filename.contact.atom

    #-- format of $filename.contact.atom
    #HISD  84CC  ALAB 103CBB 11.4961
    #PROD  85NN  ALAB 103CBB 11.0693

    #-- if condition is important. Otherwise, when grep cannot find the match, it terminates the whole bash script!
    if  egrep -q ' CA .+ CA ' $DIR/$filename.contact.atom ; then
         egrep  ' CA .+ CA ' $DIR/$filename.contact.atom > $DIR/$filename.caca.contact;
    else
        printf "\n# WARNING: Cannot find CA-CA contact pairs in $DIR/$filename.contact.atom !!\n\n"
        continue
    fi


    sed -i 's/ CA / /g ' $DIR/$filename.caca.contact
    echo "$DIR/$filename.caca.contact generated."


    if [ ! -e $DIR/$filename.caca.contact ]; then
       printf "\nERROR: $DIR/$filename.caca.contact was not generated!!!"
       exit 1
    fi
done
printf "number of pdb files under $DIR is $num_pdbFL\n\n"
echo







#-------------------------------------------------------
#-- step 3 and 4: map qry atomResNum and alignment symbols to templates' CA-CA files; add alignment symbol.

for i in `ls $DIR/*contact`;do
    # $i = '$DIR/template10.rec_lig.3tcxO_3tcxP.caca.contact'

    printf "Add qry atomResNum to the template's CA-CA file: $i \n";
    templateID=`basename $i|perl -ne '@a=split(/\./, $_); END{print $a[0]}'`

    cat $DIR/protein1.$templateID.aligned_resiNum $DIR/protein2.$templateID.aligned_resiNum > $DIR/qry_$templateID.aligned_resiNum
#    echo "perl $scriptDIR/addQryAtomResNum2templateContactFL.pl $i $DIR/qry_$templateID.aligned_resiNum >  $i.AlnSymbol"
    perl $scriptDIR/addQryAtomResNum2templateContactFL.pl $i $DIR/qry_$templateID.aligned_resiNum >  $i.AlnSymbol
    #unlink protein2.$templateID.aligned_resiNum
    #unlink protein1.$templateID.aligned_resiNum

done
echo



#-- current format of *.AlnSymbol  (after running addQryAtomResNum2templateContactFL.pl). $i.AlnSymbol later is used to generate *.contact
#
# atomResNum1_qry chnID1_qry  atomResNum2_qry chnID2_qry  Templ_aa1   Templ_chnID1    Templ_atomResNum1   Templ_aa2   Templ_chnID2    Templ_atomResNum2   Templ_dist  align_symbol
# 2   A   36  B   ALA I      2        ARG E     36        19.5321 *** ***
# 2   A   38  B   ALA I      2        GLY E     38        18.6497 *** ***
#
#--------------------------



#-- only keep the relevant CA-CA distances


for i in `ls $DIR/*contact`;do


    #--check whether the contact file is empty

    if [ ! -s $i ]; then
        printf "# WARNING: contact file $DIR/$i is empty. The template and the unbound pdb files may have a bad alignment, or have few atom contacts and no CA-CA contacts !!\n\n\n"
        continue;
    fi

    #--------
    #-- modify the following line to keep the relavant CA-CA distances into $i
    #--


    #printf "\n\n----NOTE: ONly the aligned residues (aligned far and aligned near) are included in the haddock restraints!!\n\n";

    if egrep -q "(\*\*)\s+(\*\*)" $i.AlnSymbol; then
         egrep "(\*\*)\s+(\*\*)" $i.AlnSymbol > $i.tmp
         egrep "^#" $i.AlnSymbol > $i.header
         cat $i.header $i.tmp > $i
         rm $i.tmp $i.header
    else
        printf "# WARNING: No well aligned residue pairs found. $DIR/$i will be empty.\n"
        continue
    fi

    #echo "$i.AlnSymbol filtered based on the alignment symbols (for example, removing all the unaligned residues)"
    echo "$i.AlnSymbol filtered based on the alignment symbols"

    #-- check the resulting contact file $i whether it is empty or not. If empty, the template and the unbound pdb has a very bad alignment.

    if [ ! -s $i ]; then
        printf "\n\n# WARNING: contact file $DIR/$i is empty. The template and the unbound pdb files may have a bad alignment!!\n\n\n"
    fi

done

    rm $DIR/*.AlnSymbol
    rm $DIR/*contact.atom


echo



#-- current format of contact files
#
# atomResNum1_qry chnID1_qry  atomResNum2_qry chnID2_qry  Templ_aa1   Templ_chnID1    Templ_atomResNum1   Templ_aa2   Templ_chnID2    Templ_atomResNum2   Templ_dist  align_symbol
# 2   A   36  B   ALA I      2        ARG E     36        19.5321 *** ***
# 2   A   38  B   ALA I      2        GLY E     38        18.6497 *** ***
#
#--------------------------



#-- only keep the qry columns from the contact file
for i in `ls $DIR/*contact`;do
    mv $i $i.final
    egrep '^[^#]' $i.final | awk '{print $1 " " $2 " " $3 " " $4 " " $11 }' > $i
done

#-- current format of *.caca.concact file

# atomResNum1_qry chnID1_qry  atomResNum2_qry chnID2_qry  Templ_dist
# 2   A   36  B   19.5321
# 2   A   38  B   18.6497



#-------------------------------------------------------
#- step 5:  combine the distance restraints in all *.contact

echo "Combine distance restraints in all *.contact"
final_haddockRestraintFL="$DIR/restraints_4haddock.CaCa.txt"
perl $scriptDIR/combineDisRestraintFLs_oneCluster.pl  $caDistThr $DIR $final_haddockRestraintFL



#-- clean up
rm $DIR/*caca.contact

