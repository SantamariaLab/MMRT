#!/bin/bash
homedir=$HOME/AllenModels
subdir=$1
cd $homedir/$subdir

python myRun.py 21 'Q23'
python myRun.py 21 'MMRT'

python myRun.py 34 'Q23'
python myRun.py 34 'MMRT'

python myRun.py 38 'Q23'
python myRun.py 38 'MMRT'



