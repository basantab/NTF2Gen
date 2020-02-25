1) Generate a text file containing the path to all your NTF2 backbones. The PDB files listed must contain the BENTF2DICT line. For an example of this text file see ./pdb.list

2) Generate a FASTA file with the sequences of the PDB files listed in the above file. The FASTA name must coincide with the path listed above. See ./fasta.list as an example.

3) Generate the NTF2 structural alignment using ../fit_PSSM_model/scripts/pdb2aln.py (e.g. "python3 ../fit_PSSM_model/scripts/pdb2aln.py -p pdb.list -a fasta.list" )

4) Provide the alignment file to create the PSSMs, e.g.: "python3 ./generate_PSSMs.py -alignment_file aln_seqs.fasta" PSSM files will be created in ./designs_w_pssm

5) Design the backbones using each specific PSSM file for each backbone, using ./DesignBeNTF2_test1.py. E.g.: python3 ./DesignBeNTF2_test1.py -database ~/NTF2Gen -nstruct 1 -input_pdb ./AAAGRQOG_BasicBeNTF2.pdb -threads 1 -sspredexec ~/psipred -pssm_fname ./designs_w_pssm/AAAGRQOG_BasicBeNTF2.pssm 
