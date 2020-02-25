#!/software/conda/envs/pyrosetta/bin/python3

import argparse
import numpy as np
import scipy
import scipy.spatial
import matplotlib as mpl
import string
import pyrosetta
import pandas as pd
from collections import Counter
import csv
import glob
import operator
from Bio import SeqIO
from shutil import copyfile
import os
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
aas = 'FIWLVMYCATHGSQRKPNED'

argparser = argparse.ArgumentParser(description='Generate PSSMs for biasing Rosetta design.')
argparser.add_argument('-alignment_file', type=str,help='File outputb from ../fit_PSSM_model/scripts/pdb2aln.py')
args = argparser.parse_args()

# The alignment used to fit the weights of the PSSM
alnOrigPath= '../fit_PSSM_model/aln_seqs.fasta'

# The alignment of sequences, we should generate new PSSMs for.
newAln = args.alignment_file # It is assumed that the id in this fasta file is the path to the pdb file.

pdbids = SeqIO.to_dict(SeqIO.parse(open(newAln,'r'),"fasta"))
pdb_list = [pdbid for pdbid in pdbids]

# Define needed functions
def make_pose2alnMap(alnOrig):
    posePos2OrigAlnPos = []
    for n in range(alnOrig.shape[0]): 
        pose2aln = {}
        seqOrig = alnOrig[n]
        posePos = 1
        for i,alnPos in enumerate(seqOrig):
            if alnPos!=20:
                pose2aln[posePos] = i
                posePos += 1
        posePos2OrigAlnPos.append(pose2aln)
    return posePos2OrigAlnPos

def parse_a3m(filename):
    seqs = []
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))

    # read file line by line
    for line in open(filename,"r"):
        # skip labels
        if line[0] != '>':
            # remove lowercase letters and right whitespaces
            seqs.append(line.rstrip().translate(table))
    # convert letters into numbers
    aas = 'FIWLVMYCATHGSQRKPNED-'
    alphabet = np.array(list(aas), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in seqs], dtype='|S1').view(np.uint8)
    for i in range(alphabet.shape[0]):
        msa[msa == alphabet[i]] = i

    # treat all unknown characters as gaps
    msa[msa > 20] = 20
    return msa

# Load the PSSM from the precalculated weights
alnOrig = parse_a3m(alnOrigPath)
oneBodyTerms = np.loadtxt('../fit_PSSM_model/TF_weights.txt')
X1_pssm_map = np.loadtxt('../fit_PSSM_model/TF_weights_annotation.txt', dtype='U9')
posePos2OrigAlnPos = make_pose2alnMap(alnOrig)
use_gaps = False
if use_gaps:
    n_aa_types = 21
else:
    n_aa_types = 20
pssm_length = int(len(oneBodyTerms)/n_aa_types)
pssm_height = n_aa_types
pssm = np.reshape(oneBodyTerms, [pssm_length,pssm_height])

# Figure out how many times we have seen each amino acid at each position
aa_occurance_count = np.zeros(pssm.shape)
for pos in range(0, alnOrig.shape[1]):
    aas_t, counts = np.unique(alnOrig[:,pos], return_counts=True)
    for aa,c in zip(aas_t,counts):
        if not use_gaps and aa==20:
            continue
        aa_occurance_count[pos, aa] += c

#### Make pdbs with PSSM for pymol plotting
rosetta_pssm_order = 'ARNDCQEGHILKMFPSTWYV'
pssm_aas = 'FIWLVMYCATHGSQRKPNED'
aa2no_pssm = {aa:i for i,aa in enumerate(pssm_aas)}
no2aa_pssm = {i:aa for i,aa in enumerate(pssm_aas)}
aa2no_ros = {aa:i for i,aa in enumerate(rosetta_pssm_order)}

# First we set the weight of non-observed amino acids to -100
#pssm[np.where(aa_occurance_count==0)] = -100.0
pssm[np.where(aa_occurance_count==0)] = 0.0

# There will be cases, where for all amino acids with counts, gets
# the same score (0 for instance) [CASE A]. We would probably allow design of those 
# positions (they are equally good and we don't have evidence
# that they are bad). For other positions we might have a distribution
# of amino acid scores [CASE B]. For CASE B we might want to 
# disallow design of amino acids with a score of 0. To do this,
# I set the pssm score of CASE A amino acid to 1. Then we can use
# SeqProfCon to descriminate CASE A and B, assuming that no
# weight from the PSSM fitting was equal or larger than 1.
for i, pssm_vec in enumerate(pssm):
    pssm_vec_for_present_aas = pssm_vec[np.where(aa_occurance_count[i]>0)]
    if len(pssm_vec_for_present_aas) and np.all(pssm_vec_for_present_aas == pssm_vec_for_present_aas[0]) and pssm_vec_for_present_aas[0]!=-1.:
        #pssm_vec[np.where(aa_occurance_count[i]>0)] = 100.0
        #This was giving previously unobserved regions, i.e., unaligned, very strong biases towards a specific amino-acids
        pssm[i] = pssm[i]

# Here I rescale the PSSM by the max and min of all values, such that max=10, min=-10, and leave all 0.0 as 0.0
# similar to the BLOSSUM matrix, as discussed here: https://www.rosettacommons.org/node/3615

flat_pssm2 = np.interp(oneBodyTerms, [oneBodyTerms.min(), oneBodyTerms.max()], [-10,10])
pssm2 = np.reshape(flat_pssm2, [pssm_length,pssm_height])
pssm2[np.where(pssm==0)] = 0.0

for i, pssm_vec in enumerate(pssm2):
    pssm_vec_for_present_aas = pssm_vec[np.where(aa_occurance_count[i]>0)]
    if len(pssm_vec_for_present_aas) and np.all(pssm_vec_for_present_aas == pssm_vec_for_present_aas[0]) and pssm_vec_for_present_aas[0]!=-1.:
        # Now I give a slight advantage to identities that are the only one observed for certain positions
        pssm2[i][np.where(pssm2[i]!=0)] = 1.0

# Then we iterate over all pdbs and write the pssm scores following the rules
outdir = './designs_w_pssm/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

alnNew = parse_a3m(newAln)
posePos2NewAlnPos = make_pose2alnMap(alnNew)
    
for pdbidx,pdb in enumerate(pdb_list):
    outpdb = outdir+pdb.split('/')[-1]
    copyfile(pdb, outpdb)
    pdbres2OrigAlnMap = posePos2NewAlnPos[pdbidx]
    with open(outpdb,'a') as f_open:
        f_open.write('AAorder ' + pssm_aas + '\n')
        for pdbres in pdbres2OrigAlnMap:
            origAlnResNum = pdbres2OrigAlnMap[pdbres]
            resn = no2aa_pssm[alnNew[pdbidx][origAlnResNum]]
            f_open.write('PSSM ' + str(pdbres) + ' ' + ' '.join([str(x) for x in pssm2[origAlnResNum]]) + '\n')

# Make PSSMs for rosetta
end_str = '\n\n' + '                      K         Lambda\n' + 'PSI Ungapped         0.1334     0.3157\n' + 'PSI Gapped           0.0408     0.2670\n'

for pdbidx,pdb in enumerate(pdb_list):    
    outpssmPath = outdir+pdb.split('/')[-1].replace('.pdb','.pssm')
    pdbres2OrigAlnMap = posePos2NewAlnPos[pdbidx]
    with open(outpssmPath,'w') as outpssm:
        outpssm.write('\n')
        outpssm.write('Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n')
        outpssm.write('           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n')
        for pdbres in pdbres2OrigAlnMap:
            origAlnResNum = pdbres2OrigAlnMap[pdbres]
            resn = no2aa_pssm[alnNew[pdbidx][origAlnResNum]]
            outpssm.write(' ' + str(pdbres) + ' ' + resn + ' ' + ' '.join([str(pssm2[origAlnResNum][aa2no_pssm[aa]]) for aa in rosetta_pssm_order]) + '\n')
    # Also write out the occurance counts
    outAACountsPath = outdir+pdb.split('/')[-1].replace('.pdb','.counts')
    pdbres2OrigAlnMap = posePos2NewAlnPos[pdbidx]
    with open(outAACountsPath,'w') as outCounts:
        outCounts.write('\n')
        outCounts.write('Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n')
        outCounts.write('           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n')
        for pdbres in pdbres2OrigAlnMap:
            origAlnResNum = pdbres2OrigAlnMap[pdbres]
            resn = no2aa_pssm[alnNew[pdbidx][origAlnResNum]]
            outCounts.write(' ' + str(pdbres) + ' ' + resn + ' ' + ' '.join([str(aa_occurance_count[origAlnResNum][aa2no_pssm[aa]]) for aa in rosetta_pssm_order]) + '\n')
    
