#!/bin/bash

EOS=('v2cmf' 'ptv2' 'pt2v2')
ener=('123' '200' '400' '600' '1000')
path1='/scratch/deeplearn/savchuk/potential/Skyrme/'
for eos in ${!EOS[@]}
do
    for e in ${!ener[@]}
    do
	path2=\'$path1${EOS[eos]}'/'${ener[e]}\'
	echo 'sbatch dileptons.sh '$path2' 1 60'
    done
done