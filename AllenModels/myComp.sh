#!/bin/bash

homedir=$HOME/AllenModels
subdir=$1
cd $homedir/$subdir

cp fit_parameters.json sims/MMRT/components/biophysical_neuron_models
mkdir sims/MMRT/components/mechanisms/modfiles
cp modfilesMMRT/*mod sims/MMRT/components/mechanisms/modfiles
cd sims/MMRT/components/mechanisms/
rm -rf x86_64
nrnivmodl ./modfiles

cd $homedir/$subdir
cp Biophys1.hoc sims/MMRT/components/templates
cp reconstruction.swc sims/MMRT/components/morphologies


cd $homedir/$subdir
cp fit_parameters.json sims/Q23/components/biophysical_neuron_models
mkdir sims/Q23/components/mechanisms/modfiles
cp modfiles_orig/*mod sims/Q23/components/mechanisms/modfiles
cd sims/Q23/components/mechanisms/
rm -rf x86_64
nrnivmodl ./modfiles

cd $homedir/$subdir
cp Biophys1.hoc sims/Q23/components/templates
cp reconstruction.swc sims/Q23/components/morphologies

