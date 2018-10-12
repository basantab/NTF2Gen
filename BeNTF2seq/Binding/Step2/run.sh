#!/bin/bash
module load intel-python3_2017
source /gscratch/baker/basantab/PyRosetta4.Release.python35.linux.release-175/setup.sh
export PYTHONPATH=$PYTHONPATH:/gscratch/baker/basantab/NTF2_project/BeNTF2_Koga_Rules
cd /gscratch/baker/basantab/NTF2_project/20181005_design_HCY_docks_for_chip/Step2_hybridize_hbnet
line_n=$1;
pdb=$(sed "${line_n}q;d" inputs.list );
mkdir $(printf %05d $1)
cd $(printf %05d $1) 
cat $pdb | grep -v HETATM | grep -v CONECT | grep -v HETNAM | grep -m1 -B 30000 TER > model_1.pdb;
cp $pdb .;
python3.5 ../../generate_bp_from_NTF2dict.py -database /gscratch/baker/basantab/NTF2_project/BeNTF2_Koga_Rules -input_pdb $(basename $pdb) -output_bp blueprint
python3.5 ../../FixBeNTF2Fragments_no_pocket_positions.py -input_pdb $pdb;
echo "0     A    LX    R" >> blueprint
(\time /gscratch/baker/basantab/precompiled_Rosetta/bin/rosetta_scripts.static.linuxgccrelease -s $pdb -parser:protocol ../design_files/diversify_and_HBNet.xml -parser:script_vars SS=$(cat $pdb | grep SECSTRUCT | awk '{print $2}')L -parser:script_vars lig_pos=$( /gscratch/baker/basantab/utils/get_pos.sh $pdb ) resfile=$(ls *.resfile) -parser:script_vars model_1=model_1.pdb bp=./blueprint -load_PDB_components false -extra_res_fa ../../HCY.params -out:nstruct 100 -maxruntime 7200 -mute all @../design_files/general_flags.txt) &> local.log
cd ..
