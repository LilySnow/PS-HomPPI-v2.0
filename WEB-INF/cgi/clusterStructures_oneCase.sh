# Li Xue
# July 2nd, 2014
# Given a folder of superimposed models, calculate the L-RMSDs against each other, and cluster them.

set -e

#---------global var

#-- on honavar
#profit='/home/lxue/tools/ProFitV3.1/profit'
#cluster_struc='/home/lxue/tools/cluster_struc'

#-- on alcazar
#profit='/home/software/bin/profit'
#cluster_struc='/home/software/haddock/haddock2.2/tools/cluster_struc'

#if [ ! -e $profit ];then
#    echo "$profit does not exist!!";
#    exit 1
#fi
#
#if [ ! -e $cluster_struc ];then
#    echo "$cluster_struc does not exist!!"
#    exit 1
#fi

#---------input

DIR=$1
num_pdb=`ls $DIR/*pdb |wc -l`
echo "There are $num_pdb pdb files under $DIR"

if  [ -e $DIR/l-rmsd.lst ];then
    rm $DIR/l-rmsd.lst;
fi

scriptDIR=`pwd`

cd $DIR


#--- if there is only one structure under $DIR

if [ $num_pdb == 1 ];then
    rm -rf cluster1
    mkdir cluster1

    echo "Copy files for template1 into the folder of cluster1 "
    cp template_pdbs/template1.rec_lig*.pdb cluster1
    mv *template1.aligned_resiNum cluster1

    cd $scriptDIR
    exit 0
fi




#--- if there is more than one structure under $DIR

for (( i=1; i<=$num_pdb; i++)); do
    for (( j=i+1; j<= $num_pdb; j++  ));do
        refePDBfl=`ls *template${i}.pdb`
        mobiPDBfl=`ls *template${j}.pdb`

        if [ ! -e $refePDBfl ] || [ ! -e $mobiPDBfl ];then
            printf "\nError: $refePDBfl or $mobiPDBfl does not exist !!!\n\n"
            exit 1
        fi

        printf "$i $j" >> l-rmsd.lst

        profit << END |egrep RMS | tail -1 >> l-rmsd.lst
        refe $refePDBfl
        mobi $mobiPDBfl
        atom CA,C,N,O
        zone A*
        fit
        rzone B*
        quit
END
    done
done

if [ ! -e l-rmsd.lst ];then
    echo
    echo "ERROR: $DIR/l-rmsd.lst does not exist !!!"
    echo
    exit 1
fi

awk '{print $1 " " $2 " " $4}' l-rmsd.lst > x
mv x l-rmsd.lst

echo '--------'
echo "$DIR/l-rmsd.lst generated."


#--cluster the structures based on l-rmsd.lst
l_rmsd_Thr=5
minsize=1
cluster_struc l-rmsd.lst $l_rmsd_Thr $minsize > Clusters.lst


echo "$DIR/Clusters.lst generated."


#-- copy template structures to corresponding cluster folder
if [ -d cluster1 ];then
    rm -rf cluster*
fi

for clusterID in `awk '{print $2}' Clusters.lst`;do
    mkdir cluster$clusterID
    echo "cluster$clusterID generated"

    for templateID in `awk -v row=$clusterID '{if (NR==row) {$1=$2=$3="";print $0 } }' Clusters.lst`;do
        echo "Copy files for template$templateID into the folder of cluster$clusterID "
        cp template_pdbs/template$templateID.rec_lig*.pdb cluster$clusterID
        mv *template$templateID.aligned_resiNum cluster$clusterID
    done
done



cd $scriptDIR


