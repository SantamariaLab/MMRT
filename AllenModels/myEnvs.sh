#!/bin/bash

homedir=$HOME/AllenModels
subdir=$1
iAmp=$2
cd $homedir/$subdir

python myEnv.py 21 'Q23' $iAmp
python myEnv.py 21 'MMRT' $iAmp

python myEnv.py 34 'Q23' $iAmp
python myEnv.py 34 'MMRT' $iAmp

python myEnv.py 38 'Q23' $iAmp
python myEnv.py 38 'MMRT' $iAmp


