#!/bin/bash

pdb=$1
outpath=$2

id=$(basename $pdb)
outfile=$outpath/${id%.pdb}_all_csym.dat

declare -a listrmsd=()

echo $id
echo $outpath

echo "symmetry av.rmsd clashscore" > $outfile

#for sym in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
for sym in 2 3 4 5 6 7 8 9 10 11 12
do
    rmsd=`/opt/ananas $pdb c$sym -C 100 | grep "Average RMSD" | awk '{ print $4 }' `
    #lower=`awk -v a=$rmsd 'BEGIN { print (a <= 2.5) ? "YES" : "NO" }'`
    lower=`awk -v a=$rmsd 'BEGIN { print (a <= 7) ? "YES" : "NO" }'`
    clashscore="NA"

    if [ $lower == "YES" & $sym == 2 ]
    then
        echo "symmetry c$sym has rmsd below 2.5A: $rmsd"
        clashscore=`/data4/00_BIN/phenix/phenix-1.20.1-4487/build/bin/phenix.clashscore $outpath/${id%.pdb}_c${sym}.pdb | grep "clashscore" | cut -d'=' -f2`
    
    # Is the rmsd below 2.5A threshold?
    elif [ $lower == "YES" ]
    then
        # check if the exact same rmsd is in values already computed. 
        # if so it means that we are just testing a multiple of a smaller symmetry 
        # (for ex if C3 real symmetry, the C6 and C9 have same rmsd but with clashes).
        if [[ ! " ${listrmsd[*]} " =~ " ${rmsd:0:6} " ]]
        then
            echo "symmetry c$sym has rmsd below 2.5A: $rmsd"
            # step 1: reconstruct the pdb of this specific symmetry
            /opt/ananas $pdb c$sym --symmetrize $outpath/${id%.pdb}_c${sym}.pdb
            # step 2: assess clash score with phenix.clashscore
            clashscore=`/data4/00_BIN/phenix/phenix-1.20.1-4487/build/bin/phenix.clashscore $outpath/${id%.pdb}_c${sym}.pdb | grep "clashscore" | cut -d'=' -f2`
            rm $outpath/${id%.pdb}_c${sym}.pdb # removing symmetrized pdb after evaluating clashes
        fi
    fi
    
    # append rmsd to the array
    listrmsd+=(${rmsd:0:6})

    #echo -e "\n\n*********************\n"
    #echo -e "${listrmsd[*]}"
    #echo -e "\n*********************\n\n"

    echo c$sym $rmsd $clashscore >> $outfile
done
