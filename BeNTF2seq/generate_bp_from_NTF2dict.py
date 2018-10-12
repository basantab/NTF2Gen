#!/software/miniconda3/envs/pyrosetta3/bin/python

from BeNTF2_Koga_toolkit import *
import json
import sys
import argparse

argparser = argparse.ArgumentParser(description='Design BeNTF2 backbone.')
# XMLs
argparser.add_argument('-database', type=str,help='Database folder')
argparser.add_argument('-input_pdb', nargs='+', type=str, help='sheet pdb')
argparser.add_argument('-output_bp', nargs='+', type=str, help='output bp name')
#argparser.add_argument('-sheet_dict', type=str,help='sheet dict fname')
args = argparser.parse_args()

db = args.database

for pdb,bp_name in zip(args.input_pdb,args.output_bp):
	pdb_file_handle = open(pdb,'r')
	pdb_file_lines = pdb_file_handle.readlines()
	NTF2_dict_string = ' '.join([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ][0])
	print([ line.split()[1:] for line in pdb_file_lines if 'BENTF2DICT' in line ])
	NTF2_dict = json.loads(NTF2_dict_string)
	NTF2_obj = CreateBasicNTF2fromDict(NTF2_dict,db=db)
	NTF2_obj.write_blueprint(bp_name)
