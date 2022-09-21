#!/bin/bash

homedir=$HOME/AllenModels
subdir=$1
iAmp=$2

cp ./commonfiles/Biophys1.hoc $subdir
cp -rf ./commonfiles/modfilesMMRT $subdir
cp -rf ./commonfiles/modfiles_orig $subdir
cp -rf ./commonfiles/*.py $subdir
rm -rf $subdir/sims

echo Running $subdir
./myNets.sh $subdir
./myEnvs.sh $subdir $iAmp
./myComp.sh $subdir
./myRuns.sh $subdir
