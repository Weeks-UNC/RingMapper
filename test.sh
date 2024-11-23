#!/bin/bash


function runDiff(){
if [[ ! -e $2 ]]
then
    echo "FAILED Test :: $3 output file $2 not generated"
    exit 1
fi

diff $1 $2 > d.txt 2>err.txt

if [[ -s err.txt ]] || [[ -s d.txt ]]
then
    echo "FAILED Test :: $3 output file $2 does not match reference"
    exit 1
fi

}


function checkError(){
if [[ -s $2 ]]
then
    echo "FAILED Test :: $1 exited with errors"
    exit 1
fi
}


cd testFiles

# clean up prior tests
if [[ -e err.txt ]]
then
rm err.txt 
fi

# clean up prior tests
if [[ -e d.txt ]]
then
    rm d.txt
fi


tar -xf referencefiles.tgz 

../ringmapper.py --fasta ref_addWTfull.fa --untreated ref_untreated.mut ref_modified.mut test_ring-corrs.txt --mincount 50 > test_ring.out 2>err.txt

checkError ringmapper err.txt
#runDiff ref_ring-corrs.txt test_ring-corrs.txt ringmapper
runDiff ref_ring-corrs.txt test_ring-corrs.txt_N1 ringmapper
runDiff ref_ring.out test_ring.out ringmapper


../pairmapper.py --profile ref_profile.txt --untreated ref_untreated.mut --mod ref_modified.mut --out test_pm --over --mincount 50 --secondary_react 0.5 --renorm > test_pm.out 2>err.txt


checkError pairmapper err.txt
#runDiff ref_pm.out test_pm.out pairmapper
runDiff ref_pm-allcorrs.txt test_pm-allcorrs.txt pairmapper
runDiff ref_pm-pairmap.txt test_pm-pairmap.txt pairmapper
runDiff ref_pm.dms test_pm.dms pairmapper
runDiff ref_pm.bp test_pm.bp pairmapper

echo 'All tests PASSED'

# clean up
rm ref_* err.txt d.txt

cd ../

