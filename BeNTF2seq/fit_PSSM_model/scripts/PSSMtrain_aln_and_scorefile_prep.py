import sys
import numpy as np
import scipy
import scipy.spatial
import string
import pandas as pd
from collections import Counter
import csv
import glob
import os
import seaborn as sns
from sklearn import linear_model
import pickle
import time
from optparse import OptionParser
from sklearn.model_selection import train_test_split
from scipy.cluster.hierarchy import ward, fcluster
import scipy.cluster.hierarchy as sch
from sklearn.model_selection import KFold
from scipy import stats
import matplotlib.pyplot as plt
from Bio import SeqIO
import pyrosetta
pyrosetta.init('-mute all')


################################################################################
# Options
################################################################################

parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")

# parser.add_option("-g", "--max_gap_freq", type="float", dest="maxGapFreq", metavar="FLOAT", help="Maximum frequency of gaps at position in alignment")
opts_maxGapFreq = 0

# parser.add_option("-o", "--exclude_unstable_cutoff", type="float", dest="exclude_cutoff", metavar="FLOAT", help="Do not train on design with a stability lower than cutoff")
opts_exclude_cutoff = 0.0 #None

# parser.add_option("--clustering_threshold", type="float", dest="clustering_threshold", metavar="FLOAT", help="Clustering threshold used to decide what is a different environment")
opts_clustering_threshold = 1000.0 # 0.5

# parser.add_option("--squashing_strength", type="float", dest="squashing_strength", metavar="FLOAT", help="Controls the sigmoid used to squash the distance vector. At 10 the distance is linear almost all the way. At 3, it is distance increase linearly, until around 10 AA where it flattens out.")
opts_squashing_strength = 1.25

parser.add_option("--pdbsPathFile", type="string", dest="pdbsPathFile", metavar="STR", help="...")
#opts_pdbsPathFile = 'all_pdbs.txt'

# parser.add_option("--maintain_united_rep",
#                   action="store_true", dest="maintain_united_rep", default=False,
#                   help="If environment splitting is done at a position, also retain an all-aa representation for the position")
opts_maintain_united_rep = True

# parser.add_option("--environment_clustering",
#                   action="store_true", dest="cluster_by_environments", default=False,
#                   help="If environment clustering enabled, the amino acids at one number-position of the protein, may be split into different positions, if their environment differ")
opts_cluster_by_environments = False

parser.add_option("-x", "--inputAlnFile", type="string", dest="inputAlnFile", metavar="STR", help="...")
#opts_inputAlnFile = 'beta_barrels_struc_protData.fasta'

parser.add_option("-s", "--ec50valuesPath", type="string", dest="score_and_EC50_vals", metavar="STR", help="Path to dataframe containing stability scores and ec50 values")
#opts_score_and_EC50_vals = '/home/norn/beta_barrels/190708_find_PSSM_chymo/bb_r1_chip_stability.csv'

parser.add_option("--objectiveScore", type="string", dest="objectiveScore", metavar="STR", help="Which column is score df do you want to train against (uusually stabilityscore)")

parser.add_option("--calibrate_ec50_limits",
                  action="store_true", dest="calibrate_ec50_limits", default=False,
                  help="If calibrate_ec50_limits, calibration figures will be made")
#opts_calibrate_ec50_limits = False

opts_tryp_lower_limit  = 0.5 
opts_tryp_upper_limit  = 4.5
opts_chymo_lower_limit = 0.5
opts_chymo_upper_limit = 4.5

parser.add_option("-i", "--runID", type="string", dest="runID", metavar="STR", help="Unique run ID")
    
(opts, args) = parser.parse_args()

if opts_cluster_by_environments == False and opts.inputAlnFile == None:
    raise ValueError("You need to either provide an input alignment or \
                     or let the script by the alignment for you using \
                     the cluster_by_enviroment flag") 


################################################################################
# Settings
################################################################################
outAlnPath = opts.inputAlnFile + '_' + opts.runID + '.cln'
model_pickle_filename = 'model_' + str(opts_clustering_threshold) + '_' + opts.pdbsPathFile + '_' + str(opts_squashing_strength) + '.pickle'

with open(opts.pdbsPathFile,'r') as f_open:
    pdb_list_all = [line.split()[0] for line in f_open]


################################################################################
# Util
################################################################################

aas = 'FIWLVMYCATHGSQRKPNED-'
no2aa = {i:aa for i,aa in enumerate(aas)}

################################################################################
# Functions
################################################################################

# a function to extract cb distances between all pairs of residues PyRosetta's pose
def get_pair_distances(pose, pdb):
    # number of residues in the structure
    nres = pyrosetta.rosetta.core.pose.nres_protein(pose)

    # extract coordinates of all non-hydrogen atoms
    # saving residue indices they belong to
    xyz = []
    aidx = []
    
    # three anchor atoms
    N = np.stack([np.array(pose.residue(i).atom('N').xyz()) for i in range(1,nres+1)])
    Ca = np.stack([np.array(pose.residue(i).atom('CA').xyz()) for i in range(1,nres+1)])
    C = np.stack([np.array(pose.residue(i).atom('C').xyz()) for i in range(1,nres+1)])
    
    # recreate Cb given N,Ca,C
    b = Ca - N
    c = C - Ca
    a = np.cross(b, c)
    Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca # (xyz)s over all cbeta
    
    # Calculate distance matrix
    condensed_dist_M = scipy.spatial.distance.pdist(Cb)
    dist_M = scipy.spatial.distance.squareform(condensed_dist_M)
    return dist_M

def make_dist_M_array(pdb_list):
    dist_Ms = []
    for pdb_idx,pdb in enumerate(pdb_list):
        pose = pyrosetta.pose_from_file(pdb)
        dist_M = get_pair_distances(pose, pdb)
        dist_Ms.append(dist_M)
    dist_Ms_arr = np.stack(dist_Ms)
    return dist_Ms_arr

def get_sequences(pdb_list, outAlnPath=None):
    sequences = []
    for pdb_idx,pdb in enumerate(pdb_list):
        pose = pyrosetta.pose_from_file(pdb)
        sequences.append(pose.sequence())
    
    if outAlnPath is not None:
        with open(outAlnPath,'w') as f_out:
            for pdb,seq in zip(pdb_list, sequences):
                f_out.write('>'+pdb+'\n')
                f_out.write(seq + '\n')
    
    return sequences

def filter_alignment(alnPathIn, keep_ids, alnPathOut):
    missing_in_aln = []
    
    align = SeqIO.to_dict(SeqIO.parse(alnPathIn, "fasta"))
    alignShortIds = [x.split('/')[-1].replace('.pdb','') for x in align.keys()]
    alignSeqs = [str(align[i].seq) for i in align.keys()]
    alignRenamed = {i:j for i,j in zip(alignShortIds, alignSeqs)}

    with open(alnPathOut,'w') as f_out:
        for shortId in keep_ids:
            if shortId in alignRenamed:
                f_out.write('>' + str(shortId) + '\n')
                f_out.write(str(alignRenamed[shortId]) + '\n')
            else:
                print("missing ", shortId)
                missing_in_aln.append(shortId)

def get_open_top_n_bot_idxes(df, opts_tryp_lower_limit, opts_tryp_upper_limit, opts_chymo_lower_limit, opts_chymo_upper_limit, calibrate_ec50_limits):
    # Get indices for where the designs max out the assay
    if 'stabilityscore' in opts.objectiveScore:
        is_tryp_scored = (df['stabilityscore_t_MGM'] < df['stabilityscore_c_MGM'])
        is_chymo_scored = (df['stabilityscore_t_MGM'] >= df['stabilityscore_c_MGM'])
    elif 'ec50_t' in opts.objectiveScore:
        is_tryp_scored = 1
        is_chymo_scored = 0
    elif 'ec50_c' in opts.objectiveScore:
        is_tryp_scored = 0
        is_chymo_scored = 1
    else:
        sys.exit("expected different input for objectiveScore flag")
       
    tryp_lower_limit = opts_tryp_lower_limit
    tryp_upper_limit = opts_tryp_upper_limit
    chymo_lower_limit = opts_chymo_lower_limit
    chymo_upper_limit = opts_chymo_upper_limit
    
    has_open_top    = (is_tryp_scored & (df['ec50_t'] > tryp_upper_limit)) | (is_chymo_scored & (df['ec50_c'] > chymo_upper_limit))
    has_open_bottom = (is_tryp_scored & (df['ec50_t'] < tryp_lower_limit)) | (is_chymo_scored & (df['ec50_c'] < chymo_lower_limit))
       
    if calibrate_ec50_limits:
        df['open_top'] = has_open_top
        df['open_bottom'] = has_open_bottom
        df['middle'] = (has_open_top + has_open_bottom == 0)

        g = sns.distplot(df['ec50_t'], label='ec50_t')
        g = sns.distplot(df['ec50_c'], label='ec50_c')
        fig = g.get_figure()
        fig.savefig("cal_ec50dists.png")
        plt.legend()
        plt.figure()
        
        subdftop = df[df['ec50_t'] > tryp_upper_limit]
        subdfbot = df[df['ec50_t'] < tryp_lower_limit]
        subdfmid = df[(df['ec50_t'] < tryp_upper_limit) & (df['ec50_t'] > tryp_lower_limit)]
        
        stabilityscore_postscript = ''
        if 'stabilityscore_t_MGM' in subdftop:
            stabilityscore_postscript += '_MGM'

        g = sns.scatterplot(subdftop['ec50_pred_t' + stabilityscore_postscript],subdftop['stabilityscore_t' + stabilityscore_postscript], label='openTopPoints')
        g = sns.scatterplot(subdfbot['ec50_pred_t' + stabilityscore_postscript],subdfbot['stabilityscore_t' + stabilityscore_postscript], label='openBottomPoints')
        g = sns.scatterplot(subdfmid['ec50_pred_t' + stabilityscore_postscript],subdfmid['stabilityscore_t' + stabilityscore_postscript], label='midPoints')
        fig = g.get_figure()
        fig.savefig("cal_stabScore_vs_predEC50_t.png")
    
        plt.figure()
        subdftop = df[df['ec50_c'] > chymo_upper_limit]
        subdfbot = df[df['ec50_c'] < chymo_lower_limit]
        subdfmid = df[(df['ec50_c'] < chymo_upper_limit) & (df['ec50_c'] > chymo_lower_limit)]
        g = sns.scatterplot(subdftop['ec50_pred_c' + stabilityscore_postscript],subdftop['stabilityscore_c' + stabilityscore_postscript], label='openTopPoints')
        g = sns.scatterplot(subdfbot['ec50_pred_c' + stabilityscore_postscript],subdfbot['stabilityscore_c' + stabilityscore_postscript], label='openBottomPoints')
        g = sns.scatterplot(subdfmid['ec50_pred_c' + stabilityscore_postscript],subdfmid['stabilityscore_c' + stabilityscore_postscript], label='midPoints')
        fig = g.get_figure()
        fig.savefig("cal_stabScore_vs_predEC50_c.png")
        
        plt.figure()
        subdftop = df[df['open_top']]
        subdfbot = df[df['open_bottom']]
        subdfmid = df[df['middle']]
        g = sns.distplot(subdftop['stabilityscore' + stabilityscore_postscript], label='openTopPointsc')
        g = sns.distplot(subdfbot['stabilityscore' + stabilityscore_postscript], label='openBottomPointsc')
        g = sns.distplot(subdfmid['stabilityscore' + stabilityscore_postscript], label='midPointsc')
        plt.legend()
        fig = g.get_figure()
        fig.savefig("cal_stabilityDists.png")
        
        sys.exit("Exiting as you asked for opts.calibrate_ec50_limits")
    
    return has_open_top, has_open_bottom


################################################################################
# Setup for training
################################################################################
# load the protease stabilities for the list of pdbs provided
pdblistShortNamesInput = [n.split('/')[-1].replace('.pdb','') for n in pdb_list_all] # split('_')[0] this might be needed for chip2...

# Load the score file
df = pd.read_csv(opts.score_and_EC50_vals) # skiprows=range(0,9)
df = df.dropna(subset=['ec50_c', 'ec50_t'])
df = df[df['name'].isin(pdblistShortNamesInput)]

# Do not use data with high credibility intervals
df = df[(df['ec50_95ci_c']*np.log10(3) < 2) & (df['ec50_95ci_t']*np.log10(3) < 2)]

# Annotate which sequences that max out the assay
has_open_top, has_open_bottom = get_open_top_n_bot_idxes(df, opts_tryp_lower_limit, opts_tryp_upper_limit, opts_chymo_lower_limit, opts_chymo_upper_limit, opts.calibrate_ec50_limits)

# Update the pdblistShortNames and let user know if we dropped any pdbs
pdblistShortNames = [n for n in pdblistShortNamesInput if n in df['name'].tolist()]
input_len = len(pdblistShortNamesInput)
filtered_len = len(pdblistShortNames)
if input_len != filtered_len:
    print("Dropped", input_len-filtered_len,"sequences as they had missing stability values or high 95ci intervals")

scores = np.array(df[opts.objectiveScore].tolist()) # 6882 
keep_ids = [x.split('/')[-1].replace('.pdb','') for x in df['name'].tolist()]

# Filter the alignment
filter_alignment(opts.inputAlnFile, keep_ids, outAlnPath)

# Dump the filtered data
d = {'scores': scores, 'name': keep_ids, 'has_open_top': has_open_top, 'has_open_bottom': has_open_bottom}
df_out = pd.DataFrame.from_dict(d)
df_out.to_csv('annotated_stability_data.csv')



print(len(scores))
print(len(df))
print(len(keep_ids))

print("Training files now pickled")
