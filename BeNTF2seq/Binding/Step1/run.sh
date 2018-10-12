#!/bin/bash
module load intel-python3_2017
source /gscratch/baker/basantab/PyRosetta4.Release.python35.linux.release-175/setup.sh
export PYTHONPATH=$PYTHONPATH:/gscratch/baker/basantab/NTF2_project/BeNTF2_Koga_Rules
cd /gscratch/baker/basantab/NTF2_project/20181005_design_HCY_docks_for_chip/Step1_dock_design
line_n=$1;
pdb=$(sed "${line_n}q;d" inputs.list );
mkdir $(printf %05d $1)
cd $(printf %05d $1)
cp $pdb .;
python3.5 ../../FixBeNTF2Fragments_no_pocket_positions.py -input_pdb $pdb;
(\time /gscratch/baker/basantab/precompiled_Rosetta/bin/rosetta_scripts.static.linuxgccrelease -s $pdb -parser:protocol ../design_files_penalty/ala_redes_pack_binder.xml -parser:script_vars SS=$(cat $pdb | grep SECSTRUCT | awk '{print $2}')L -parser:script_vars lig_pos=$( /gscratch/baker/basantab/utils/get_pos.sh $pdb ) resfile=$(ls *.resfile) -load_PDB_components false -extra_res_fa ../../HCY.params -mute all @../design_files_penalty/general_flags.txt) &> local.log
cd ..
