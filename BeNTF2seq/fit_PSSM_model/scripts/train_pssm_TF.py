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

start = time.time()
print("Started modeling at time", start)

# import keras
import tensorflow as tf
#import keras.backend as K
# from keras.models import Sequential
# from keras.layers import Dense
# from keras import regularizers




################################################################################
# Options
################################################################################

parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")

parser.add_option("-r", "--l1ratio", type="float", dest="l1ratio", metavar="STR", help="...")
#opts.l1ratio = 0.0001

parser.add_option("-a", "--alpha", type="float", dest="alpha", metavar="STR", help="...")


parser.add_option("--cv_mode",
                  action="store_true", dest="cv_mode", default=False,
                  help="Run cross-validation on training split")


parser.add_option("-x", "--inputAlnFile", type="string", dest="inputAlnFile", metavar="STR", help="...")


parser.add_option("-i", "--runID", type="string", dest="runID", metavar="STR", help="Unique run ID")

parser.add_option("-v", "--learningRate", type="float", dest="learningRate", metavar="STR", help="Unique run ID")



(opts, args) = parser.parse_args()

# opts.runID = 'test'
# opts.inputAlnFile = 'beta_barrels_struc_protData.fasta_test.cln'
# opts.cv_mode = True
# opts.alpha = 0.0005
# opts.l1ratio = 0.0001
# opts.learningRate = 0.001


################################################################################
# Util
################################################################################

aas = 'FIWLVMYCATHGSQRKPNED-'
no2aa = {i:aa for i,aa in enumerate(aas)}

################################################################################
# Train function
################################################################################
def tf_fit(alpha, l1_ratio, aln1hot_train, normalized_scores_train, open_top_idxes, open_bottom_idxes, learningRate, epochs=500, couplings=False, test_aln1hot = None):
    features = aln1hot_train.shape[1]
    nseqs = aln1hot_train.shape[0]
    tf.reset_default_graph()

    x = tf.placeholder(tf.float32, shape=(None, features)) 
    y = tf.placeholder(tf.float32, shape=(None,))
    top_idxes = tf.constant(open_top_idxes, dtype=tf.float32)
    bottom_idxes = tf.constant(open_bottom_idxes, dtype=tf.float32)
    y_pred = tf.layers.dense(x, 1, activation=tf.identity, name="dense")
    d = y_pred[:, 0] - y
    
    # We use a loss function which does not penalize with stability scores
    # that are above the maximum that the assay can give
    err_top = tf.multiply(tf.multiply(d, top_idxes), tf.dtypes.cast(d<0, dtype=tf.float32))
    err_bot = tf.multiply(tf.multiply(d, bottom_idxes), tf.dtypes.cast(d>0, dtype=tf.float32))
    err_middle = d * tf.constant((open_top_idxes + open_bottom_idxes) == 0, dtype=tf.float32)
    mse_mod_loss = tf.reduce_mean(tf.pow(err_top + err_bot + err_middle, 2))

    l1_weight = alpha*l1_ratio*2
    l2_weight = (1-l1_ratio) * alpha
    
    wts = tf.trainable_variables()
    loss_L2 = tf.reduce_sum([ tf.nn.l2_loss(w) for w in wts if 'bias' not in w.name ]) * l2_weight
    loss_L1 = tf.reduce_sum([ tf.math.abs(w) for w in wts if 'bias' not in w.name ]) * l1_weight
    cost = mse_mod_loss + loss_L2 + loss_L1
    
    # Optimize
    train = tf.train.AdamOptimizer(learningRate).minimize(cost) 

    # instantiating Session
    config = tf.ConfigProto(gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.9))
    sess = tf.Session(config=config)
    sess.run(tf.global_variables_initializer())

    # Train
    losses = []
    for i in range(1, epochs+1):
        if i % int(epochs/10) == 0:
            print(i)

        # Use all the data in every epoch
        sess.run(train, feed_dict={x:aln1hot_train, y:normalized_scores_train})
        losses.append(sess.run(cost, feed_dict={x:aln1hot_train, y:normalized_scores_train}))
        
        if (i > 200):
            if losses[i-101]-losses[i-1] < 0.001:
                break

    # Return
    weights = sess.run(tf.trainable_variables())
    if test_aln1hot is not None:
        predictions = sess.run(y_pred, feed_dict={x:test_aln1hot, y:np.zeros(test_aln1hot.shape[0])}) 
        return sess, predictions, weights, losses
    else:
        return weights

# read A3M and convert letters into
# integers in the 0..20 range
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

################################################################################
# Load data and train
################################################################################

# Load the annotated protease stability
df = pd.read_csv('annotated_stability_data.csv')
names = df['name']
scores = np.array(df['scores'])
has_open_top = np.array(df['has_open_top'])
has_open_bottom = np.array(df['has_open_bottom'])

# Load the alignment
use_gaps = False 
aln = parse_a3m(opts.inputAlnFile)
M=aln.shape[0]
L=aln.shape[1]
ns = 21 if use_gaps else 20
print("Cleaned alignment has", M, "entries with ", L, "positions.")

# Setup a 1hot matrix for primary sequence and a map for that WITH GAPS
if use_gaps:
    X1 = np.reshape(np.eye(21)[aln], [M,L*ns])
else:
    X1 = np.reshape(np.eye(21)[aln][:,:,:ns], [M,L*ns])
X1_pssm_map = np.empty(X1.shape[1], dtype='U9')
pos = 0
for i in range(L):
    for j in range(0, ns):
        aa = no2aa[j]
        X1_pssm_map[pos] = str(i) + aa
        pos += 1
print("We will be fitting ", X1.shape, " weights for 1 body terms")

# Check that everything is alright
m = len(scores)
print("Scorefile has " + str(m) + " scores")
if m != M:
    raise ValueError("The number of protease stabilities " + str(m) +" does NOT match with the number of sequences "+str(M)) 

print(X1_pssm_map.shape)
print(X1.shape)
print(scores.shape)
print(has_open_top.shape)
print(has_open_bottom.shape)

# Split into training and test
aln1hot_train, aln1hot_test                 =  train_test_split(X1, test_size=0.10, random_state=42)
scores_train_raw, scores_test_raw           =  train_test_split(scores, test_size=0.10, random_state=42)
has_open_top_train, has_open_top_test       =  train_test_split(has_open_top, test_size=0.10, random_state=42)
has_open_bottom_train, has_open_bottom_test =  train_test_split(has_open_bottom, test_size=0.10, random_state=42)
indices_train, indices_test                 =  train_test_split(range(0, len(scores)), test_size=0.10, random_state=42)

# Normalize the data
ave = np.mean(scores_train_raw)
std = np.std(scores_train_raw)
scores_train = (scores_train_raw-ave)/std
scores_test = (scores_test_raw-ave)/std

if opts.cv_mode:
    kf = KFold(n_splits=5, random_state=42)
    
    scores = []
    for split_idx, (train_index, test_index) in enumerate(kf.split(scores_train)):
        kfold_aln1hot_train = aln1hot_train[train_index,:]
        kfold_scores_train = scores_train[train_index]
        kfold_aln1hot_test = aln1hot_train[test_index,:]
        kfold_scores_test = scores_train[test_index]

        sess, predictions, weights, losses = tf_fit(opts.alpha, opts.l1ratio, kfold_aln1hot_train, 
                                            kfold_scores_train, has_open_top_train[train_index], 
                                            has_open_bottom_train[train_index], 
                                            opts.learningRate, epochs=5000, couplings=False, test_aln1hot=kfold_aln1hot_test)
        
        score_cv_test = stats.pearsonr(predictions.flatten(), kfold_scores_test)[0] ** 2
        scores.append(str(score_cv_test))
    print("TF CV", opts.l1ratio, opts.alpha, opts.learningRate, opts.runID, ','.join(scores))
else:
    sess, predictions, weights, losses = tf_fit(opts.alpha, opts.l1ratio, aln1hot_train, 
                                        scores_train, has_open_top_train, 
                                        has_open_bottom_train, 
                                        opts.learningRate, epochs=5000, test_aln1hot=aln1hot_test)
    score_test = stats.pearsonr(predictions.flatten(), scores_test)[0] ** 2
    print("R sqaured on holdout set", score_test)

    # Make some diagnostics plots
    plt.plot(losses)
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.savefig("convergence_against_holdout.png")

    # Then train on the full dataset
    weights = tf_fit(opts.alpha, opts.l1ratio, X1,
                                        (scores-ave)/std, has_open_top,
                                        has_open_bottom,
                                        opts.learningRate, epochs=5000)
    
    # And save the weights
    weights[0].shape
    oneBodyTerms = weights[0]*std + ave/len(weights[0])/21
    np.savetxt("TF_weights.txt", oneBodyTerms)
    np.savetxt("TF_weights_annotation.txt", X1_pssm_map, fmt='%s')


