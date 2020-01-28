from BeNTF2_toolkit import *
import json
import argparse
import subprocess
import copy
from Blueprint import *
import sys
import glob
import pandas as pd
from Bio import SeqIO
from optparse import OptionParser

sys.path.append("/home/basantab/NTF2Gen")
sys.path.append("/home/basantab/fragment_tools")
sys.path.append("/home/basantab/scripts")


################################################################################
# Options
################################################################################

parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")
parser.add_option("-p", "--pdb_list", type="string", dest="pdblistPath", metavar="STR", help="list with all paths to relevant pdbs")
parser.add_option("-a", "--fastaPath", type="string", dest="fastaPath", metavar="STR", help="path to fasta files with all sequences")
(opts, args) = parser.parse_args()

#-----------------------------------

fastas = opts.fastaPath
pdblistPath = opts.pdblistPath
alnout = 'aln_seqs.fasta'

with open(pdblistPath, 'r') as f_open:
    pdbs = [line.split()[0] for line in f_open]

#################################

# Get dictionaries from the pdbs
pdb_data = {}
for PDB_fname in pdbs:
    pdb_handle = open(PDB_fname,'r')
    print(PDB_fname)
    dict_line = [line[:-1] for line in pdb_handle.readlines() if 'BENTF2DICT' in line ][0]
    pdb_handle.close()
    dict_string = ' '.join(dict_line.split()[1:])
    BeNTF2dict = json.loads(dict_string)
    for key in BeNTF2dict['ring_dict']:
        if 'sheet_dict' in key:
            continue
        BeNTF2dict[key] = BeNTF2dict['ring_dict'][key]
    for key in BeNTF2dict['ring_dict']['sheet_dict']:
        BeNTF2dict[key] = BeNTF2dict['ring_dict']['sheet_dict'][key]
    for key in BeNTF2dict['c_helix_dict']:
        BeNTF2dict[key+'c_helix_dict'] = BeNTF2dict['c_helix_dict'][key]
    pdb_data[PDB_fname] = BeNTF2dict

df = pd.DataFrame.from_dict(pdb_data, orient='index')

# Get sequences
with open(fastas, "rt") as f_open:
    pdb_seqs = SeqIO.parse(f_open, "fasta")
    pdb_seqs_str = {str(x.id): str(x.seq) for x in pdb_seqs}

# Remove the sequences that are not in the fasta file
df = df[df.index.isin(list(pdb_seqs_str.keys()))]

# The necessary lengths for making the alignment
connection_type2length = {'BA':1, 'GBA': 2, 'GB': 2, 'ClassicDirect': 0, 'BulgeAndB': 4, 'BBGB': 3}
connection_type_order = ['BA', 'GBA', 'GB', 'ClassicDirect', 'BulgeAndB','BBGB']
miscounts_E3 = ['BulgeAndB', 'ClassicDirect', 'GBA','BA', 'BBGB'] # 
H3_types = [15, 14, 10, 11]
total_connection_type_length_at_h3len = [8, 4, 4, 8]
connection_types_length_at_h3len = {15: ['BulgeAndB', 'GBA', 'GB'], 
                                    14: ['BA', 'ClassicDirect', 'BBGB'], 
                                    10: ['ClassicDirect', 'BBGB', 'BA'], 
                                    11: ['GBA', 'GB', 'BulgeAndB']}

E3_insertion_positions_rel2_E3end = [-16, -15, -14, -13, -12, -11, -10, -8, -6, -5]


max_H1 = 23
max_H2 = 11
max_hairpin = 6
max_H3 = 15
max_E3 = 16 + 2
max_E4 = 14 + 2
max_E5 = 16 + 2
max_E6 = 11 + 2
max_term_helix_len = 11
h3_connection_segmentation = False

f_aln_out = open(alnout, 'w')


for pdbidindex,pdbid in enumerate(pdb_seqs_str):
    aln_seq = ''
    gaps_inserted = 0
    
    name = pdbid.split('/')[-1].split('_')[0]
    
    # H1 length variation
    h1_len = df.loc[pdbid]['h1_len']
    gaps_h1 = max_H1-h1_len
    aln_seq += gaps_h1*'-'
    aln_seq += pdb_seqs_str[pdbid][0:h1_len]
    gaps_inserted += gaps_h1
   
    # H2 length variation
    h2_len = df.loc[pdbid]['h2_len']
    gaps_h2 = max_H2-h2_len
    aln_seq += gaps_h2*'-'
    h2_end = h1_len + h2_len
    E1_start = h2_end + 8
    aln_seq += pdb_seqs_str[pdbid][h1_len:E1_start]
    gaps_inserted += gaps_h2
    
    # Hairpin length variation
    hairpin_len = df.loc[pdbid]['hairpin_len']
    hairpin_partA_gaps = max_hairpin - hairpin_len
    E1 = pdb_seqs_str[pdbid][E1_start:E1_start+hairpin_len]
    E1E2_loop = pdb_seqs_str[pdbid][E1_start+hairpin_len:E1_start+hairpin_len+2]
    E2 = pdb_seqs_str[pdbid][E1_start+hairpin_len+2:E1_start+2*hairpin_len+2]
    aln_seq += E1   
    aln_seq += hairpin_partA_gaps*'-'
    aln_seq += E1E2_loop
    aln_seq += hairpin_partA_gaps*'-'
    aln_seq += E2
    E2_end = E1_start+2*hairpin_len+2
    
    # H3 length variation
    h3_len = df.loc[pdbid]['h_len']
    h3_gaps = max_H3 - h3_len
    aln_seq += pdb_seqs_str[pdbid][E2_end:E2_end+1+h3_len]
    aln_seq += '-' * h3_gaps
    h3_end = E2_end+1+h3_len
        
    # Connection type variation
    connection_type = df.loc[pdbid]['connection_type']
    this_len = connection_type2length[connection_type] #+ miscount_E3
    gaps_inserted = 0
    for c_type in connection_type_order:
        if c_type!=connection_type:
            c_gaps = connection_type2length[c_type]
            aln_seq += '-' * c_gaps
        else:
            aln_seq += pdb_seqs_str[pdbid][h3_end:h3_end+this_len]

    ################### E3 ######################
    # I align this based on the alignment of AAB* and resi 65-88, 
    # which alignment E3 and E4. This let's me insure that the 
    # up and down pattern of the beta-stand is intact. The
    # alternative is to figure out whether the strand start
    # with up or down, but that results in unalignable bulges
    # the patterning of which is possibly somewhat environment
    # independent.
    miscount_E3 = 1 if connection_type in miscounts_E3 else 0
    E3_start = h3_end + this_len
    sec_bulge = df.loc[pdbid]['sec_bulge_E3_pos'] + miscount_E3
    if connection_type=='GBA':
        sec_bulge = 1
    if connection_type=='BA':
        sec_bulge = 1
    E3_len = df.loc[pdbid]['E3_l'] + miscount_E3 
    E3_end = E3_start + E3_len 
    has_not_sec_bulge = 0 if sec_bulge>0 else 1
    
    init_E3_gaps = max_E3 - E3_len - has_not_sec_bulge
    main_bulge = df.loc[pdbid]['main_bulge_E3_pos'] 
    E3_w_gaps = '-' * init_E3_gaps + pdb_seqs_str[pdbid][E3_start:E3_end] 

    rel_mbulge_position = main_bulge - 1 - E3_len
    rel_sbulge_position = sec_bulge - 1 - E3_len
    
    E3 = ''
    for i, c in enumerate(E3_w_gaps[::-1]):
        i = -i - 1 # we count from -1 and down
        if i in E3_insertion_positions_rel2_E3end and (i!=rel_mbulge_position and i!=rel_sbulge_position):
            E3 += '-'
        E3 += c
    
    # Some very E3s are not being measured correctly for some reason.
    # The fudge added due to this is assigned to the end of the strand.
    # One day we should figure out why this is necessary.
    if len(E3) == 26:
        E3 += '-'
    aln_seq += E3[::-1]
    
    # E3E4
    E3E4 = pdb_seqs_str[pdbid][E3_end:E3_end+2]
    aln_seq += E3E4
    
    ################### E4 ######################
    # Now we change reference frame to the global
    # alignment. This way at least a few of the
    # H1 facing positions will be structurally 
    # aligned
    long_arm_l = df.loc[pdbid]['long_arm_l'] 
    basewidth = df.loc[pdbid]['base_width']
    E4_c_term_gaps = int((26-(((long_arm_l*2)+1)*2+2+(basewidth-3)*2+2))/2)
    E4_len = df.loc[pdbid]['E4_l']
    E4_gaps_total = max_E4 - E4_len
    E4_n_term_gaps = int(E4_gaps_total - E4_c_term_gaps+2)
    aln_seq += '-' * E4_n_term_gaps

    E4_start = E3_end+2
    E4_end = E4_start + E4_len
    aln_seq += pdb_seqs_str[pdbid][E4_start:E4_start+E4_len]
    aln_seq += '-' * E4_c_term_gaps
    
    ################## E4E5 ######################
    E4E5 = pdb_seqs_str[pdbid][E4_end:E4_end+2]
    aln_seq += E4E5
    
    ################### E5 ######################
    E5_start = E4_end+2
    E5_len = df.loc[pdbid]['E5_l']
    E5_end = E5_start + E5_len
    E5_gaps_cterm = (2 - df.loc[pdbid]['short_arm_l'])*2
    E5_gaps_nterm = E4_c_term_gaps
    E5 = pdb_seqs_str[pdbid][E5_start:E5_start+E5_len]
    aln_seq += '-' *E5_gaps_nterm + E5 + '-' *E5_gaps_cterm
    
    # E5E6
    E5E6 = pdb_seqs_str[pdbid][E5_end:E5_end+2]
    aln_seq += E5E6
    

    ################## E6 ######################
    E6_start = E5_end+2
    E6_n_term_gaps = E5_gaps_cterm
    E6_len = df.loc[pdbid]['E6_l']
    is_hungry_tropical_pitcher = (df.loc[pdbid]['Opening']=='Classic' and df.loc[pdbid]['has_cHelix'])
    #c_helix_connetion_type = df.loc[pdbid]['loopABEGOc_helix_dict']
    has_helix = df.loc[pdbid]['h1_lenc_helix_dict'] > 0

    total_length = len(pdb_seqs_str[pdbid])
    E6_end = E6_start+E6_len
    term_helix_len = df.loc[pdbid]['h1_lenc_helix_dict'] if df.loc[pdbid]['h1_lenc_helix_dict'] > 0 else 0
    E6_2_helix_linker_length = int(total_length - E6_end - term_helix_len - 1)
    #print(pdbid, E6_2_helix_linker_length)
    
    if not has_helix:
        E6 = pdb_seqs_str[pdbid][E6_start:E6_start+E6_len+1]
    else:
        E6 = pdb_seqs_str[pdbid][E6_start:E6_start+E6_len]
    aln_seq += '-'*E6_n_term_gaps + E6

    # Gaps to complete the c-term E6 length variation
    if not has_helix:
        E6_gaps_c_term = max_E6 - E6_len - E6_n_term_gaps + 2
    else:
        E6_gaps_c_term = max_E6 - E6_len - E6_n_term_gaps + 3

    aln_seq += '-' * E6_gaps_c_term


    ################### c-term_helix ######################
    c_helix_type = df.loc[pdbid]['Opening']
    
    ctermHelix = ''

    # Add the E6 to helix linker
    aln_seq += '-' * (2-E6_2_helix_linker_length)
    aln_seq += pdb_seqs_str[pdbid][E6_end:E6_end+E6_2_helix_linker_length]
    
    # Add the helix
    helix_gaps = max_term_helix_len - term_helix_len
    if has_helix:
        aln_seq += pdb_seqs_str[pdbid][E6_end+E6_2_helix_linker_length:]
    else:
        aln_seq += '-'
    aln_seq += '-' * int(max_term_helix_len - term_helix_len)
    
    if len(aln_seq) != 183:
        print("There is a sequence, which is not getting the expected length. Time for debugging...")
        print(aln_seq, pdbid[-40:], len(aln_seq))
           
    f_aln_out.write('>'+pdbid+'\n')
    f_aln_out.write(aln_seq+'\n')


f_aln_out.close()




