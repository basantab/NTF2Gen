#!/software/conda/envs/pyrosetta/bin/python3.7

from BeNTF2_toolkit import *
import json
import sys


def compare_SheetDict(sheet1,sheet2):
	#print("Comparing sheets...")
	sheet_defining_attributes = [ key for key in NTF2_sheet.sheet_atributes_Dict.keys() ]
	for attr in sheet_defining_attributes:
		if sheet1[attr] != sheet2[attr]:
			#print("They are different")
			return False
	#print("They are the same")
	return True

def compare_RingDict(ring_dict1,ring_dict2):
	#print("Comparing rings...")
	ring_defining_attributes = RingConnection.ring_attributes
	for attr in ring_defining_attributes:
		if ring_dict1[attr] != ring_dict2[attr]:
			#print("They are different")
			return False
	#print("They are the same")
	return True

def compare_BeNTF2FinalPartsDict(dict1,dict2):
	#print("Comparing Final NTF2 parts...")
	BeNTF2_defining_attributes = ["loopTwoABEGO", "Opening","h1_len","h2_len","has_cHelix"]
	for attr in BeNTF2_defining_attributes:
		if dict1[attr] != dict2[attr]:
			#print("They are different")
			return False
	#print("They are the same")
	return True

def compare_BeNTF2Dict(dict1,dict2):
	sheet_dict1 = dict1["ring_dict"]["sheet_dict"]
	sheet_dict2 = dict2["ring_dict"]["sheet_dict"]
	sheet_same = compare_SheetDict(sheet_dict1,sheet_dict2)
	if not sheet_same: return False
	ring_dict1 = dict1["ring_dict"]
	ring_dict2 = dict2["ring_dict"]
	ring_same = compare_RingDict(ring_dict1,ring_dict2)
	if not ring_same: return False
	BeNTF2_same = compare_BeNTF2FinalPartsDict(dict1,dict2)
	if not BeNTF2_same: return False
	return True

def compare_two(fname1, fname2):
	pdb1_handle = open(fname1,'r')
	pdb2_handle = open(fname2,'r')
	dict1_str = ' '.join([line[:-1] for line in pdb1_handle.readlines() if 'BENTF2DICT' in line ][0].split()[1:])
	dict2_str = ' '.join([line[:-1] for line in pdb2_handle.readlines() if 'BENTF2DICT' in line ][0].split()[1:])
	pdb1_handle.close()
	pdb2_handle.close()
	dict1 = json.loads(dict1_str)
	dict2 = json.loads(dict2_str)
	return compare_BeNTF2Dict(dict1,dict2)

def get_dict_fname(fname):
	pdb1_handle = open(fname,'r')
	try: dict1_str = ' '.join([line[:-1] for line in pdb1_handle.readlines() if 'BENTF2DICT' in line ][0].split()[1:])
	except: print('Offending file: %s'%fname)
	pdb1_handle.close()
	dict1 = json.loads(dict1_str)
	return dict1

if __name__ == "__main__":
	fname_list = sys.argv[1]
	fname_list_handle = open(fname_list,'r')
	PDBs = [ i[:-1] for i in fname_list_handle.readlines() ]
	fname_list_handle.close()
	print("Done reading file names")
	full_list = {}
	PDBs_dicts = { fname:get_dict_fname(fname) for fname in PDBs }
	print("Done converting files to dictionaries")
	n_cluster = 0
	for n,i in enumerate(sorted(PDBs_dicts.keys())):
		are_the_same = []
		uniq_one = True
		for target in sorted([ key for key in PDBs_dicts.keys()])[:n]:
			if compare_BeNTF2Dict( PDBs_dicts[i], PDBs_dicts[target]) :
				full_list[i] = full_list[target]
				uniq_one = False
				break
		if uniq_one:
			full_list[i] = n_cluster
			n_cluster += 1
	print("Unique dict NTF2s")
	for name,val in full_list.items():
		print("%s,%d"%(name,val))
