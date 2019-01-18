#!/software/miniconda3/envs/pyrosetta3/bin/python

import sys
from glob import glob

# Take dock filename as input, search for the original scaffold. Open both, retrieve two comment lines from original,
# write new file with same name as dock but inseting the two comment lines plus a blank line right before the
# first REMARK line.

dock_fname = sys.argv[1]
prefix = dock_fname.split('/')[-1].split('_')[0]
original_fname = glob('/home/basantab/NTF2_project/20180411_designBBsInHyak/finished_designs_by_length*/???/%s*.pdb'%prefix)[0]

dock_handle = open(dock_fname,'r')
dock_lines = dock_handle.readlines()
dock_handle.close()

ori_handle = open(original_fname,'r')
ori_lines = ori_handle.readlines()
ori_handle.close()

important_comment_lines = []
for n,line in enumerate(ori_lines):
	if "##Begin comments##" in line:
		important_comment_lines.append(ori_lines[n+1])
		important_comment_lines.append(ori_lines[n+2])
		important_comment_lines.append(ori_lines[n+3])
		important_comment_lines.append('\n')
		break

first_remark_i = 0
for n,line in enumerate(dock_lines):
	if "##End comments##" in line:
		first_remark_i = n
		break

for line in reversed(important_comment_lines):
	dock_lines.insert(first_remark_i,line)

output_fname = dock_fname.split('/')[-1]
output_handle = open(output_fname,'w')
for line in dock_lines: output_handle.write(line)
output_handle.close()

