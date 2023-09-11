#!/bin/bash


# Path to the script
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
echo -n "path to script: $SCRIPTPATH"


pdb=$1
outpath=$2

if [ ! -d $outpath ]
then
    mkdir $outpath
fi

id=$(basename $pdb)
outfile=$outpath/${id%.pdb}_all_csym.dat

declare -a listrmsd=()

echo $id
echo $outpath

echo "symmetry av.rmsd clashscore" > $outfile

for sym in 2 3 4 5 6 7 8 9 10 11 12
do
    rmsd=`$ANANAS $pdb c$sym -C 100 | grep "Average RMSD" | awk '{ print $4 }' `
    lower=`awk -v a=$rmsd 'BEGIN { print (a <= 7) ? "YES" : "NO" }'`
    clashscore="NA"
    echo -e "\n*************************\nrmsd : $rmsd\n*************************\n"

    if [ $lower == "YES" ] && [ $sym == 2 ]
    then
        echo -e "\n*************************\nsymmetry c$sym has rmsd below 2.5A: $rmsd\n*************************\n"
        # The model has a C2 symmetry, no need to test other symmetries
        break
    
    # Is the rmsd below 2.5A threshold?
    elif [ $lower == "YES" ]
    then
        echo -e "\n*************************\nsymmetry c$sym has rmsd below 2.5A: $rmsd\n*************************\n"
        
        # check if the exact same rmsd is in values already computed. 
        # if so it means that we are just testing a multiple of a smaller symmetry 
        # (for ex if C3 real symmetry, the C6 and C9 have same rmsd but with clashes).
        if [[ ! " ${listrmsd[*]} " =~ " ${rmsd:0:6} " ]]
        then
            echo "symmetry c$sym has rmsd below 2.5A: $rmsd"
            
            # step 1: reconstruct the pdb of this specific symmetry
            $ANANAS $pdb c$sym --symmetrize $outpath/${id%.pdb}_c${sym}.pdb
            
            # step 2: assess clash score with phenix.clashscore
            clashscore=`$PHENIX_CLASHSCORE $outpath/${id%.pdb}_c${sym}.pdb | grep "clashscore" | cut -d'=' -f2`
            
            rm $outpath/${id%.pdb}_c${sym}.pdb # removing symmetrized pdb after evaluating clashes
        fi
    fi
    
    # append rmsd to the array
    listrmsd+=(${rmsd:0:6})

    echo c$sym $rmsd $clashscore >> $outfile
done
