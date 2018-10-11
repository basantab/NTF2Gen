#!/software/miniconda3/envs/pyrosetta3/bin/python3

from BeNTF2_toolkit import *
import json
import argparse
import copy
from Blueprint import *

def get_SecBulgePos(bp,BeNTF2_obj):
	sheet_obj = BeNTF2_obj.ring_obj.sheet_obj
	E3 = bp.segment_dict['E3']
	if sheet_obj.ExtendedE4:
		main_bulge_pos_from_C = -1*(sheet_obj.base_width)-3
	else:
		main_bulge_pos_from_C = -1*(sheet_obj.base_width)-1
	second_bulge_position_from_C = main_bulge_pos_from_C - ( 2*sheet_obj.Second_b_place + 3 )
	return E3.bp_data[second_bulge_position_from_C][0]

########## Option system: ###########
argparser = argparse.ArgumentParser(description='Create resfile for optimal fragment quality in all places outside of the NTF2 pocket')
argparser.add_argument('-database', type=str, default='/home/basantab/scripts/BeNTF2_Koga_Rules' , help='Database path for accesory files')
argparser.add_argument('-input_pdb', type=str,help='Input NTF2 file with dict string')
args = argparser.parse_args()

db = args.database
PDB_fname = args.input_pdb
pdb_handle = open(PDB_fname,'r')
dict_line = [ line[:-1] for line in pdb_handle.readlines() if 'BENTF2DICT' in line ][0]
pdb_handle.close()
dict_string = ' '.join(dict_line.split()[1:])
BeNTF2dict = json.loads(dict_string)
BeNTF2_obj = CreateBasicNTF2fromDict(BeNTF2dict,NTF2_is_complete=True,db=db)

init()
POSE = pose_from_file( PDB_fname )

# ABEGO string from pose:
abego_str = core.sequence.get_abego(POSE,1)

# Complete blueprint:
bp_data = copy.copy(BeNTF2_obj.NTF2_bp.bp_data)
bp = Blueprint(data=bp_data)
bp.reindex_blueprint()

SS = ''.join([ pos[2][0] for pos in bp.bp_data ])
print("Pose from %s:"%PDB_fname)
print(SS)
print(abego_str)

for n,pos in enumerate(bp.bp_data):
	bp.bp_data[n][2] = pos[2][0]+abego_str[n+1]

# Resfile
resfile_handle = open('%s.resfile'%(PDB_fname.split('/')[-1]),'w')
resfile_handle.write('ALLAA\n')
resfile_handle.write('START\n')

# Bulges:

# Main bulge positions: B-A-B-B
E3 = bp.segment_dict['E3']
MainBulge_A = max([ pos[0] for pos in E3.bp_data if pos[2][1] == 'A' ]); MainBulge_A_aa = 'DEK' # Simple heuristic to get the main bulge, the main bulge is always the last one in E3
MainBulge_pre = MainBulge_A - 1; MainBulge_pre_aa = 'IV'
MainBulge_B = MainBulge_A + 1; MainBulge_B_aa = 'RES'
MainBulge_post = MainBulge_A + 2;  MainBulge_post_aa = 'VI'

resfile_handle.write( '%d A PIKAA %s\n'%(MainBulge_A,MainBulge_A_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(MainBulge_pre,MainBulge_pre_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(MainBulge_B,MainBulge_B_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(MainBulge_post,MainBulge_post_aa) )

# E6 bulge:

E6 = bp.segment_dict['E6']
E6Bulge_A = min([ pos[0] for pos in E6.bp_data if pos[2][1] == 'A' ]); E6Bulge_A_aa = 'E'
if BeNTF2_obj.ring_obj.sheet_obj.short_arm_l == 1: E6Bulge_A_aa = 'KR'
E6Bulge_pre = E6Bulge_A - 1; E6Bulge_pre_aa = 'IV'
E6Bulge_B = E6Bulge_A + 1; E6Bulge_B_aa = 'SE'
E6Bulge_post = E6Bulge_A + 2;  E6Bulge_post_aa = 'VI'

resfile_handle.write( '%d A PIKAA %s\n'%(E6Bulge_A,E6Bulge_A_aa) )
#resfile_handle.write( '%d A PIKAA %s\n'%(E6Bulge_pre,E6Bulge_pre_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(E6Bulge_B,E6Bulge_B_aa) )
#resfile_handle.write( '%d A PIKAA %s\n'%(E6Bulge_post,E6Bulge_post_aa) )

# Seconary bulge if present:
if BeNTF2_obj.ring_obj.sheet_obj.Second_bulge_E3:
	SecBulge_A = get_SecBulgePos(bp,BeNTF2_obj); SecBulge_A_aa = 'EDH'
	SecBulge_pre = SecBulge_A - 1; SecBulge_pre_aa = 'IV'
	SecBulge_B = SecBulge_A + 1; SecBulge_B_aa = 'EDH'
	SecBulge_post = SecBulge_A + 2; SecBulge_post_aa = 'IV'

	resfile_handle.write( '%d A PIKAA %s\n'%(SecBulge_A,SecBulge_A_aa) )
	#resfile_handle.write( '%d A PIKAA %s\n'%(SecBulge_pre,SecBulge_pre_aa) )
	resfile_handle.write( '%d A PIKAA %s\n'%(SecBulge_B,SecBulge_B_aa) )
	#resfile_handle.write( '%d A PIKAA %s\n'%(SecBulge_post,SecBulge_post_aa) )

# Loop #2:

loop2 = bp.segment_dict['L2']
loop2_C = loop2.bp_data[-1][0]; loop2_C_aa = 'D'
resfile_handle.write( '%d A PIKAA %s\n'%(loop2_C,loop2_C_aa) )

# Loop #3:

loop3 = bp.segment_dict['L3']
loop3_N = loop3.bp_data[0][0]; loop3_N_aa = 'LIMVA'
loop3_C = loop3_N + 3; loop3_C_aa = 'ND'
loop3_C_1 = loop3_C + 1; loop3_C_1_aa = 'AVT'

resfile_handle.write( '%d A PIKAA %s\n'%(loop3_C,loop3_C_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(loop3_C_1,loop3_C_1_aa) )

# Loop #5:

loop5 = bp.segment_dict['L5']
loop5_N = loop5.bp_data[0][0] - 1; loop5_N_aa = 'H'
loop5_C = loop5_N + 2; loop5_C_aa = 'RW'
loop5_H = loop5_N + 4; loop5_H_aa = 'E'

resfile_handle.write( '%d A PIKAA %s\n'%(loop5_N,loop5_N_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(loop5_C,loop5_C_aa) )
resfile_handle.write( '%d A PIKAA %s\n'%(loop5_H,loop5_H_aa) )

# 2-residue Strand-Strand loops:

for key in bp.segment_dict.keys():
	if (key[0] == 'L') and (len(key) == 2) :
		abego = ''.join([ pos[2][1] for pos in bp.segment_dict[key].bp_data ])
		if abego == 'GG':
			pos_1st = bp.segment_dict[key].bp_data[0][0]
			resfile_handle.write( '%d A PIKAA N\n'%(pos_1st) )
		elif abego == 'EA':
			pos_2nd = bp.segment_dict[key].bp_data[1][0]
			resfile_handle.write( '%d A PIKAA D\n'%(pos_2nd) )

# Frontal hairpin:

abego = ''.join([ pos[2][1] for pos in bp.segment_dict['L4'].bp_data ])
if abego == 'GG':
	loop4 = bp.segment_dict['L4']
	loop4_N = loop4.bp_data[0][0] - 1; loop4_N_aa = 'F'
	#resfile_handle.write( '%d A PIKAA %s\n'%(loop4_N,loop4_N_aa) )

	if BeNTF2_obj.ring_obj.hairpin_len == 4:
		HP_acid_1 = loop4_N - 1; HP_acid_1_aa = 'E'
		resfile_handle.write( '%d A PIKAA %s\n'%(HP_acid_1,HP_acid_1_aa) )
		HP_acid_2 = loop4_N + 4; HP_acid_2_aa = 'E'
		resfile_handle.write( '%d A PIKAA %s\n'%(HP_acid_2,HP_acid_2_aa) )
		HP_base = loop4_N + 3; HP_base_aa = 'RQK'
		resfile_handle.write( '%d A PIKAA %s\n'%(HP_base,HP_base_aa) )

# Short Arm:
L9 = bp.segment_dict['L9']
L9_1st_pos = L9.bp_data[0][0]
if BeNTF2_obj.ring_obj.sheet_obj.short_arm_l == 1:
	L9_base = L9_1st_pos + 2; L9_base_aa = 'RK'
	resfile_handle.write( '%d A PIKAA %s\n'%(L9_base,L9_base_aa) )
else:
	L9_base = L9_1st_pos + 2; L9_base_aa = 'RKQ'
	resfile_handle.write( '%d A PIKAA %s\n'%(L9_base,L9_base_aa) )
	L9_aro = L9_1st_pos - 1; L9_aro_aa = 'FY'
	resfile_handle.write( '%d A PIKAA %s\n'%(L9_aro,L9_aro_aa) )
	L9_leu = L9_1st_pos - 3; L9_leu_aa = 'LWF'
	resfile_handle.write( '%d A PIKAA %s\n'%(L9_leu,L9_leu_aa) )
	L9_aro = L9_1st_pos + 4; L9_aro_aa = 'IVF'
	resfile_handle.write( '%d A PIKAA %s\n'%(L9_aro,L9_aro_aa) )
	# Loop 7 K:
	L7 = bp.segment_dict['L7']
	L7_1st_pos = L7.bp_data[0][0]
	L7_base = L7_1st_pos + 2; L7_base_aa = 'K'
	resfile_handle.write( '%d A PIKAA %s\n'%(L7_base,L7_base_aa) )

# H3-E3 loops:
L6 = bp.segment_dict['L6']
L6_1st_pos = L6.bp_data[0][0]

#if BeNTF2_obj.ring_obj.connection_type=='BA':
# Sequences seem to be optimal already	

if BeNTF2_obj.ring_obj.connection_type=='GBA':
	L6_D = L6_1st_pos + 2; L6_D_aa = 'D'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_D,L6_D_aa) )
	L6_H = L6_1st_pos + 3; L6_H_aa = 'H'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_H,L6_H_aa) )

#if BeNTF2_obj.ring_obj.connection_type=='GB':
# Already optimal	

if BeNTF2_obj.ring_obj.connection_type=='ClassicDirect':
	L6_base = L6_1st_pos - 2; L6_base_aa = 'KRS'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_base,L6_base_aa) )
	L6_I = L6_1st_pos ; L6_I_aa = 'I'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_I,L6_I_aa) )

if BeNTF2_obj.ring_obj.connection_type=='BulgeAndB':
	L6_V = L6_1st_pos + 1; L6_V_aa = 'V'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_V,L6_V_aa) )
	L6_D = L6_1st_pos + 2 ; L6_D_aa = 'D'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_D,L6_D_aa) )

if BeNTF2_obj.ring_obj.connection_type=='BBGB':
	L6_PK = L6_1st_pos + 1; L6_PK_aa = 'PK'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_PK,L6_PK_aa) )
	L6_G = L6_1st_pos + 2 ; L6_G_aa = 'G'
	resfile_handle.write( '%d A PIKAA %s\n'%(L6_G,L6_G_aa) )
	
resfile_handle.close()

