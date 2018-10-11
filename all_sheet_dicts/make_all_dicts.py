#!/software/miniconda3/envs/pyrosetta3/bin/python3

from BeNTF2toolkit import *
import json

for n,d in CreateAllPossibleSheetDicts():
	fname = "%04d_sheet.dict"%n
	handle = open(fname,'w')
	json.dump(d,handle)
	handle.close()
