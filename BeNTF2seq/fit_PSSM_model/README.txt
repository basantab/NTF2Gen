# Make the alignments of designs with well-estimated EC50 values (contained in pdbPaths.txt)
/software/conda/envs/pyrosetta/bin/python scripts/pdb2aln.py -p pdbPaths.txt -a chip12_all_seqs.fasta

# Make the inputfiles for training using lower and upper bounds for the EC50 values
/software/conda/envs/pyrosetta/bin/python scripts/PSSMtrain_aln_and_scorefile_prep.py --objectiveScore stabilityscore_MGM --runID MGM --inputAlnFile aln_seqs.fasta --ec50valuesPath ../molten_globule_model/chip12_comp_USM_n_stability.csv --pdbsPathFile pdbPaths.txt

# Run hyperparameter optimization while fitting weights of the PSSM
python scripts/submit_hyperparameter_search.py > commands
# Run those commands (most likely on a HPC-cluster)

# Then analyse whether the hyperparameter search found a solution on the
# edge of the tested grid. If it did, extend the grid search changing the
# ranges in submit_hyperparameter_search.py
grep "TF CV" cv_logs/* | sed 's/,/ /g' > cv_summary.txt
# analysis_cv.ipynb will calculate mean scores, and identify the best hyperparameter. 
# Then train again, but this time on the entire dataset, not cross-validating: 
source activate /net/software/conda/envs/tensorflow; python scripts/train_pssm_TF.py -r 0.004642 -a 0.077426 --cv_mode --runID MGM --inputAlnFile aln_seqs.fasta_MGM.cln --learningRate 0.002154 | grep "R sqaured on holdout set"  > final_score.log

# To analyse the resulting weights, run analyse_models_TF.ipynb

# To add PSSM to pdbs
BENJA write how here

