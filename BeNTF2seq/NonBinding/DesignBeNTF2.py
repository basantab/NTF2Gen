#!/sw/intel-python-2017.2.045/intelpython3/bin/python3

from BeNTF2_toolkit import *
import json
import sys
import argparse
import random, string
import pickle
import numpy as np
import random
import os
import copy
import glob

argparser = argparse.ArgumentParser(description='Design BeNTF2 backbone.')
# XMLs
argparser.add_argument('-database', type=str,help='Database folder')
# trials
argparser.add_argument('-nstruct', type=int,help='Max number of structures to output. May output less if some of them fail.')
argparser.add_argument('-prefix', type=str,help='Prefix to add to output names.')
argparser.add_argument('-input_pdb', type=str,help='sheet pdb')
#argparser.add_argument('-sheet_dict', type=str,help='sheet dict fname')
args = argparser.parse_args()
#out_fname = sys.argv[1]
#dict_name = sys.argv[2]

#handle = open(dict_name,'r')
#sheet_dict = json.load(handle)
#handle.close()

def randomword(length):
	return ''.join(random.choice(string.ascii_uppercase) for i in range(length))

def MakeMoverFromXML(xmlfile,replaces,pose):
	print("Will attempt to create mover now")
	parser = protocols.rosetta_scripts.RosettaScriptsParser()
	replaces_container = utility.vector1_string(0)

	for i in replaces:
		replaces_container.append(i)

	modified_pose = False
	tag = parser.create_tag_from_xml( xmlfile, replaces_container)
	in_mover = parser.generate_mover_for_protocol( pose, modified_pose, tag, options )
	print("Propterly created mover")
	return in_mover

def add_SS_labels(POSE,basic_NTF2Obj):
	'''
	info = POSE.pdb_info()
	for i,pos in enumerate(basic_NTF2Obj.NTF2_bp.bp_data):
		resn = i+1
		ss_label = "SS_%s"%(pos[2][0])
		info.add_reslabel(resn,ss_label)
	'''
	SS_list = [i[2][0] for i in basic_NTF2Obj.NTF2_bp.bp_data]
	SS_string = ''.join(SS_list)
	core.pose.add_comment(POSE,"SECSTRUCT",SS_string)

def add_BeNTF2_dict_comment(POSE,dict_line):
	core.pose.add_comment(POSE,"BENTF2DICT",dict_line)

def detect_pocket_residues(POSE,basic_NTF2Obj):
	bp_data = copy.copy(basic_NTF2Obj.NTF2_bp.bp_data)
	bp = Blueprint(data=bp_data)
	bp.reindex_blueprint()
	if basic_NTF2Obj.NTF2_dict["Opening"] == "Classic":
		E2_N_pos = bp.segment_dict['E1'].bp_data[0][0] + 3
		E3_N_pos = bp.segment_dict['E3'].bp_data[0][0]

		E3_N_CA = POSE.residue(E3_N_pos).xyz('CA')
		E2_N_CA = POSE.residue(E2_N_pos).xyz('CA')

		midpoint = numeric.xyzVector_double_t( (E3_N_CA[0]+E2_N_CA[0])/2, (E3_N_CA[1]+E2_N_CA[1])/2, (E3_N_CA[2]+E2_N_CA[2])/2)

	else:
		H2_N_pos = bp.segment_dict['H2'].bp_data[0][0]
		E1_C_pos = bp.segment_dict['E1'].bp_data[0][0] + 3
		
		H2_N_CA = POSE.residue(H2_N_pos).xyz('CA')
		E1_C_CA = POSE.residue(E1_C_pos).xyz('CA')
		
		midpoint = numeric.xyzVector_double_t( (H2_N_CA[0]+E1_C_CA[0])/2, (H2_N_CA[1]+E1_C_CA[1])/2, (H2_N_CA[2]+E1_C_CA[2])/2)

	virt_type = core.pose.virtual_type_for_pose(POSE)
	virt = core.conformation.ResidueFactory.create_residue(virt_type)
	virt.set_xyz('ORIG',midpoint)
	last_res = len(POSE)
	POSE.append_residue_by_jump(virt,last_res)
	protein_selector = core.select.residue_selector.ResidueIndexSelector()
	virt_atm_selector = core.select.residue_selector.ResidueIndexSelector()
	protein_selector.set_index('1-%d'%last_res)
	virt_atm_selector.set_index('%d'%(last_res+1))
	ss_selector = core.select.residue_selector.SecondaryStructureSelector()
	ss_string = ('').join([i[2][0] for i in bp.bp_data])+'L'
	ss_selector.set_pose_secstruct(ss_string)
	ss_selector.set_selected_ss('HE')
	grp1_selector = core.select.residue_selector.AndResidueSelector()
	grp1_selector.add_residue_selector(ss_selector)
	grp1_selector.add_residue_selector(protein_selector)
	interf_selector = core.select.residue_selector.InterGroupInterfaceByVectorSelector()
	interf_selector.group1_selector(grp1_selector)
	interf_selector.group2_selector(virt_atm_selector)
	interf_selector.cb_dist_cut(20.0)
	interf_selector.vector_dist_cut(16.0)
	interf_selector.vector_angle_cut(90.0)
	interf_selector.nearby_atom_cut(0.0)
	positions_bool = interf_selector.apply(POSE)
	position_vector = []
	for i,n in enumerate(positions_bool):
		if n:
			position_vector.append(i+1)
	return position_vector

def add_pocket_labels(POSE,positions):
	info = POSE.pdb_info()
	pocket_label = 'Pckt'
	for i in positions:
		info.add_reslabel(i,pocket_label)

def add_SS_positions(POSE,NTF2_obj):
	info = POSE.pdb_info()
	H1_label = 'H1pos'
	H1_positions = NTF2_obj.get_H1_positions()
	for i in H1_positions:
		info.add_reslabel(i,H1_label)
	
	if NTF2_obj.NTF2_dict["has_cHelix"]:
		H4_label = 'H4pos'
		H4_positions = NTF2_obj.get_H4_positions()
		for i in H4_positions:
			info.add_reslabel(i,H4_label)
	
	long_HP_label = 'longHPpos'
	long_HP_positions = NTF2_obj.get_long_arm_inward_pos()
	for i in long_HP_positions:
		info.add_reslabel(i,long_HP_label)

def add_H3_positions(POSE,NTF2_obj):
	info = POSE.pdb_info()
	H3_label = 'H3'
	H3_positions = NTF2_obj.get_H3_pos()
	for i in H3_positions:
		info.add_reslabel(i,H3_label)

def get_angle(pose,pos):
	pos_ca = pose.residue(pos).xyz('CA')
	print("Got gly position: %d"%pos)
	pos_minus_ca = pose.residue(pos-2).xyz('CA')
	pos_plus_ca = pose.residue(pos+2).xyz('CA')
	vec1 = pos_minus_ca - pos_ca
	vec2 = pos_plus_ca - pos_ca
	vec1.normalize_or_zero()
	vec2.normalize_or_zero()
	angle = math.degrees( math.acos( rosetta.numeric.sin_cos_range(vec2.dot_product(vec1)) ) )
	return angle

def RescueMainBulge(POSE,NTF2_obj):
	# Phe position at the main bulge
	gly_pos = NTF2_obj.get_gly_resc_gly()[0]
	angle_at_gly = get_angle(POSE,gly_pos)
	print("Measuring Main bulge angle:")
	print("Angle: %0.3f"%angle_at_gly)
	if angle_at_gly < 147.5:
		print("returning True")
		return True
	else:
		print("REturning Flase")
		return False

def RescueSecBulge(POSE,NTF2_obj):
	# Phe position at the main bulge
	if NTF2_obj.NTF2_dict['ring_dict']['sheet_dict']['Second_bulge_E3']:
		gly_pos = NTF2_obj.get_gly_resc_gly()[1]
		angle_at_gly = get_angle(POSE,gly_pos)
		if angle_at_gly < 147.5:
			return True
		else:
			return False
	else:
		return False

def RescueCurvedLongArm(POSE,NTF2_obj):
	if NTF2_obj.NTF2_dict['ring_dict']['sheet_dict']['CurvedLongArm']:
		gly_pos = NTF2_obj.get_gly_resc_gly()[1]
		angle_at_gly = get_angle(POSE,gly_pos)
		if angle_at_gly < 147.5:
			return True
		else:
			return False
	else:
		return False

def ExtremeCurve(POSE,NTF2_obj):
	phe_pos = NTF2_obj.get_gly_resc_phe()[1]
	angle_at_phe = get_angle(POSE,phe_pos)
	if angle_at_phe < 140.0:
		return True
	else:
		return False

def KeepBackbone(data,NTF2_dict):
	cacheable_data = data.get(2)
	avdeg_H3 = cacheable_data.map()['avdeg_H3']
	avdeg_H1 = cacheable_data.map()['avdeg_H1']
	is_TP = 1 if NTF2_dict['Opening'] == 'Tropical' else 0
	print("BB features for keeping or not: avdeg_H3: %0.3f avdeg_H1: %0.3f is_TP: %d"%(avdeg_H3, avdeg_H1, is_TP))
	if avdeg_H3 <= 10.08:
		if is_TP:
			return True
		else:
			return False
	else:
		if avdeg_H1 <= 10.13:
			return False
		else:
			return True

def GetSecondPositionInAllBulges(NTF2_obj):
	sheet_obj = NTF2_obj.ring_obj.sheet_obj
	has_second_bulge = sheet_obj.Second_bulge_E3
	bp = Blueprint( data=[ [ j for j in i ] for i in NTF2_obj.NTF2_bp.bp_data ] )
	bp.reindex_blueprint()
	bp.freeze_all()
	E3 = bp.segment_dict['E3']
	positions = {}
	if sheet_obj.ExtendedE4:
		main_bulge_pos_from_C = -1*(sheet_obj.base_width)-3
	else:
		main_bulge_pos_from_C = -1*(sheet_obj.base_width)-1
	positions['main_bulge'] = E3.bp_data[main_bulge_pos_from_C][0]

	E6 = bp.segment_dict['E6']
	E6_bulge_INDEX_from_N = 2*(sheet_obj.short_arm_l)
	positions['E6_bulge'] = E6.bp_data[E6_bulge_INDEX_from_N][0]
	
	if has_second_bulge:
		if sheet_obj.ExtendedE4:
			main_bulge_pos_from_C = -1*(sheet_obj.base_width)-3
		else:
			main_bulge_pos_from_C = -1*(sheet_obj.base_width)-1
		second_bulge_position_from_C = main_bulge_pos_from_C - ( 2*sheet_obj.Second_b_place + 3 )
		positions['SecBulge'] = E3.bp_data[second_bulge_position_from_C][0]
	
	return positions

def GetLoop2Positions(NTF2_obj):
	bp = Blueprint( data=[ [ j for j in i ] for i in NTF2_obj.NTF2_bp.bp_data ] )
	bp.reindex_blueprint()
	bp.freeze_all()
	positions = {}
	L2 = bp.segment_dict['L2']
	H1c = bp.segment_dict['H1'].bp_data[-1][0]
	positions['H1c'] = H1c
	positions['L2_1'] = L2.bp_data[0][0]
	positions['L2_2'] = L2.bp_data[1][0]
	return positions

def GetLoop3Positions(NTF2_obj):
	bp = Blueprint( data=[ [ j for j in i ] for i in NTF2_obj.NTF2_bp.bp_data ] )
	bp.reindex_blueprint()
	bp.freeze_all()
	positions = {}
	L3 = bp.segment_dict['L3']
	positions['L3_2'] = L3.bp_data[1][0]
	positions['L3_3'] = L3.bp_data[2][0]
	return positions

def GetTRPlockPositions(NTF2_obj):
	bp = Blueprint( data=[ [ j for j in i ] for i in NTF2_obj.NTF2_bp.bp_data ] )
	bp.reindex_blueprint()
	bp.freeze_all()
	positions = {}
	E5 = bp.segment_dict['E5']
	E6 = bp.segment_dict['E6']
	E4 = bp.segment_dict['E4']
	positions['trp'] = E6.bp_data[1][0]
	positions['lys'] = E5.bp_data[-2][0]
	positions['small'] = E4.bp_data[1][0]
	return positions

def DesignStep(POSE,NTF2_obj,ssstring,NTF2_dict,step=None):
	replaces = ['SS=%s'%ssstring]

	if step == 1:
		replaces.append('aa_comp=%s/general2.comp'%db)
		xmlfile = '%s/DesignStage1.xml'%db
		#NTF2_obj.write_blueprint('test1.bp')
	elif step == 2:
		replaces.append('aa_comp=%s/general2.comp'%db)
		xmlfile = '%s/DesignStage2_hyak.xml'%db
	elif step == 3:
		replaces.append('aa_comp=%s/general2.comp'%db)
		xmlfile = '%s/DesignStage3_hyak.xml'%db
	else:
		raise SystemExit('Invalid stage selection! Exiting')

	info = POSE.pdb_info()

	if NTF2_obj.ring_obj.sheet_obj.short_arm_l == 2:
		TRP_lock_positions = GetTRPlockPositions(NTF2_obj)
		info.add_reslabel(TRP_lock_positions['lys'],'lys')
		info.add_reslabel(TRP_lock_positions['trp'],'trp')
		info.add_reslabel(TRP_lock_positions['small'],'small')
	
	phe_positions = NTF2_obj.get_gly_resc_phe()
	gly_positions = NTF2_obj.get_gly_resc_gly()

	L3_positions = GetLoop3Positions(NTF2_obj)
	info.add_reslabel(L3_positions['L3_2'],'L3_2')
	info.add_reslabel(L3_positions['L3_3'],'L3_3')

	L2_positions = GetLoop2Positions(NTF2_obj)
	info.add_reslabel(L2_positions['H1c'],'H1c')
	info.add_reslabel(L2_positions['L2_1'],'L2_1')
	info.add_reslabel(L2_positions['L2_2'],'L2_2')	

	bulge_positions = GetSecondPositionInAllBulges(NTF2_obj)
	bulge_label = 'bulge_A'
	info.add_reslabel(bulge_positions['main_bulge'],bulge_label)
	info.add_reslabel(bulge_positions['E6_bulge'],bulge_label)
	if NTF2_obj.ring_obj.sheet_obj.Second_bulge_E3:
		info.add_reslabel(bulge_positions['SecBulge'],bulge_label)
	
	if RescueMainBulge(POSE,NTF2_obj):
		#phe_pos_label = 'GRF'
		gly_pos_label = 'GRG'
		#info.add_reslabel(phe_positions[0],phe_pos_label)
		info.add_reslabel(gly_positions[0],gly_pos_label)

	if RescueSecBulge(POSE,NTF2_obj):
		#phe_pos_label = 'GRF_sec'
		gly_pos_label = 'GRG_sec'
		#info.add_reslabel(phe_positions[1],phe_pos_label)
		info.add_reslabel(gly_positions[1],gly_pos_label)
		#if ExtremeCurve(POSE,NTF2_obj):
		#	replaces.append('GRG_sec=F')
		#	replaces.append('GRF_sec=GDN')
		#else:
		replaces.append('GRG_sec=G')
		replaces.append('GRF_sec=F')
	else:
		replaces.append('GRG_sec=G')
		replaces.append('GRF_sec=F')

	if RescueCurvedLongArm(POSE,NTF2_obj):
		gly_pos_label = 'GRG'
		info.add_reslabel(gly_positions[1],gly_pos_label)
		
	replaces.append('hard_fa=beta_nov16')
	replaces.append('sspredexec=/gscratch/baker/basantab/utils/psipred4.01/runpsipred_csbuild_single')
	replaces.append('longHPcomp=%s/longHP.comp'%db)
	replaces.append('sheet_comp=%s/sheet.comp'%db)
	replaces.append('H1comp=%s/H1.comp'%db)
	replaces.append('pocketcomp=%s/pocket.comp'%db)
	replaces.append('pocket_charge=%s/neutral.charge'%db)
	replaces.append('charge=%s/negative.charge'%db)
	replaces.append('surface_Helix=%s/surface_Helix_filter.comp'%db)
	replaces.append('surface_Strand=%s/surface_Strand_filter.comp'%db)
	replaces.append('boundary_Helix=%s/boundary_Helix_filter.comp'%db)
	replaces.append('boundary_Strand=%s/boundary_Strand_filter.comp'%db)
	replaces.append('core_Helix=%s/core_Helix_filter.comp'%db)
	replaces.append('core_Strand=%s/core_Strand_filter.comp'%db)
	replaces.append('core_comp=%s/core.comp'%db)
	replaces.append('surface_comp=%s/surface.comp'%db)
	mover = MakeMoverFromXML(xmlfile,replaces,POSE)
	try:mover.apply(POSE)
	except RuntimeError:
		print('Something went wrong while running mover')
		return False
	else:
		print('Done applying mover')
		data = POSE.data()
		if mover.get_last_move_status() in failed_mover_statuses :
			print('Pose did not pass filters inside mover')
			handle = open('./step3_done','w')
			handle.close()
			return False
		else:
			if step == 3:
				print('Succesfully applied mover 3!')
				handle = open('./step3_done','w')
				handle.close()
				return True
			if step == 1:
				if KeepBackbone(data,NTF2_dict):
					print('Succesfully applied mover 1 and passed filters')
					POSE.dump_file( '%s_BasicBeNTF2_step1_%04d.pdb'%(prefix,trial) )
					return True
				else:
					POSE.dump_file( 'REJECTED_%s_BasicBeNTF2_%04d.pdb'%(prefix,trial) )
					print('Did not pass criterion for well-formed backbone')
					return False
			if step == 2:
				POSE.dump_file( '%s_BasicBeNTF2_step2_%04d.pdb'%(prefix,trial) )
				print('Succesfully applied mover 2!')
				return True

# Initialize Rosetta:

db = args.database

general_flags = ['-beta',\
		'-out:file:pdb_comments',\
		'-corrections:beta_nov16',\
		'-out:file:pdb_comments',\
		'-mute all']

prefix = ''

init( extra_options=" ".join(general_flags) )
options = basic.options.process()
failed_mover_statuses_names = [ "FAIL_RETRY", "FAIL_DO_NOT_RETRY", "FAIL_BAD_INPUT"]
failed_mover_statuses = [ protocols.moves.mstype_from_name( name ) for name in failed_mover_statuses_names ]

for trial in range(args.nstruct):
	'''
	if not args.prefix:
		prefix = ''
	else:
		prefix = args.prefix
	prefix = randomword(8)+'_'+prefix
	'''
	prefix = args.input_pdb.split('/')[-1].split('.')[0]
	print("Prefix: %s Trial: %d"%(prefix,trial))
	
	POSE = pose_from_file( filename=args.input_pdb)
	pdb_file_handle = open(args.input_pdb,'r')
	pdb_file_lines = pdb_file_handle.readlines()
	ssstring = [ line.split()[1] for line in pdb_file_lines if 'SECSTRUCT' in line ][0]
	print("Secondary structure string: %s"%ssstring)
	NTF2_dict_string = ' '.join([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ][0])
	print([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ])
	NTF2_dict = json.loads(NTF2_dict_string)
	NTF2_obj = CreateBasicNTF2fromDict(NTF2_dict,db=db)
	NTF2_obj.NTF2_bp.reindex_blueprint()
	add_H3_positions(POSE,NTF2_obj)
	#### STEP 1 ####
	step1_files = glob.glob('*_BasicBeNTF2_step1_*.pdb')
	step2_files = glob.glob('*_BasicBeNTF2_step2_*.pdb')
	if (len(step1_files) == 0) and (len(step2_files) == 0 ) :
		success = DesignStep(POSE,NTF2_obj,ssstring,NTF2_dict,step=1)
		if not success: continue
		print('Passed step 1')
	elif (len(step2_files) == 0 ):
		print("Step1 files found: %s"%(' '.join(step1_files)) )
		if not (len(step1_files) == 1):
			raise ValueError('You have more than one step1 file in this dir!')
		prefix = args.input_pdb.split('/')[-1].split('.')[0]
		POSE = pose_from_file( filename=step1_files[0])
		pdb_file_handle = open(step1_files[0],'r')
		pdb_file_lines = pdb_file_handle.readlines()
		ssstring = [ line.split()[1] for line in pdb_file_lines if 'SECSTRUCT' in line ][0]
		print("Secondary structure string: %s"%ssstring)
		NTF2_dict_string = ' '.join([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ][0])
		print([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ])
		NTF2_dict = json.loads(NTF2_dict_string)
		NTF2_obj = CreateBasicNTF2fromDict(NTF2_dict,db=db)
		NTF2_obj.NTF2_bp.reindex_blueprint()
	#### STEP 2 ####
	if len(step2_files) == 0:
		success = DesignStep(POSE,NTF2_obj,ssstring,NTF2_dict,step=2)
		if not success: continue
		print('Passed step 2')
	else:
		print("Step2 files found: %s"%(' '.join(step2_files)) )
		if not (len(step2_files) == 1):
			raise ValueError('You have more than one step2 file in this dir!')
		prefix = args.input_pdb.split('/')[-1].split('.')[0]
		POSE = pose_from_file( filename=step2_files[0])
		pdb_file_handle = open(step2_files[0],'r')
		pdb_file_lines = pdb_file_handle.readlines()
		ssstring = [ line.split()[1] for line in pdb_file_lines if 'SECSTRUCT' in line ][0]
		print("Secondary structure string: %s"%ssstring)
		NTF2_dict_string = ' '.join([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ][0])
		print([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ])
		NTF2_dict = json.loads(NTF2_dict_string)
		NTF2_obj = CreateBasicNTF2fromDict(NTF2_dict,db=db)
		NTF2_obj.NTF2_bp.reindex_blueprint()
	#### STEP 3 ####
	success = DesignStep(POSE,NTF2_obj,ssstring,NTF2_dict,step=3)
	if not success: continue
	print('Passed step 3')

	try: POSE.dump_file( '%s_BasicBeNTF2_designed_%04d.pdb'%(prefix,trial) )
	except RuntimeError:
		print('Unable to create file, skipping')
		continue
