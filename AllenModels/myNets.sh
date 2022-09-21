#!/bin/bash

homedir=$HOME/AllenModels
subdir=$1
cd $homedir/$subdir


python myNetwork.py 21 'Q23'
#python myNetwork.py 34 'Q23'
#python myNetwork.py 38 'Q23'
#
#python myNetwork.py 21 ''
#python myNetwork.py 34 ''
#python myNetwork.py 40 ''
#
