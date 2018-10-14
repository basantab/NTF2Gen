#!/bin/bash

module load intel-python3_2017
source /gscratch/baker/basantab/PyRosetta4.Release.python35.linux.release-195/setup.sh
export PYTHONPATH=$PYTHONPATH:/gscratch/baker/basantab/NTF2Gen

python3.5 ../DesignBeNTF2.py -database /gscratch/baker/basantab/NTF2Gen -input_pdb ./AADMDBFA_AHYDXLMV_0.pdb -nstruct 1
