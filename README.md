# NTF2Gen: An enumerative algorithm for full-atom models of proteins belonging to the NTF2-like superfamily<br/>
<br/>

If you find the scripts in this repository useful, please cite https://doi.org/10.1101/2020.03.23.003913.<br/>
<br/>

This algorithm samples a wide diversity of protein structures by carrying out backbone
sampling at two levels. At the top level, sampling is carried out in the space of high-level
parameters that define the overall properties of the NTF2 fold: for example, the overall sheet length and
curvature, the lengths of the helices that complement the sheet, the placement of the pocket opening and
the presence or absence of C-terminal elements. We then convert each choice of high-level
parameters into structure blueprint/constraints pairs, which guide
backbone structure sampling at successive stages of fold assembly. In a final sequence design step, for each generated
backbone, low energy sequences are identified through combinatorial sequence optimization using
RosettaDesign.<br/>
<br/>

The NTF2Gen repository contains all the tools for de novo design of NTF2-like proteins. The main script
is CreateBeNTF2_backbone.py, which manages the construction of NTF2 backbones, followed by
DesignBeNTF2.py (BeNTF2seq/Nonbinding, or DesignBeNTF2_test1.py at BeNTF2seq/design_with_PSSM to design using PSSMs), which designs sequence on a given backbone generated by the previous script. To
generate backbones from a specific set of parameters, use CreateBeNTF2PDBFromDict.py.<br/>
<br/>

The fundamental building blocks of the backbone generation protocol are Rosetta XML protocols
(included in the repository) that are specialized instances of the BlueprintBDRMover Rosetta fragment
assembly mover. All backbone quality checks and filters previous to design are
implemented either in the XML files or the python scripts. The design script is also based on a set of XML
protocols, one for each design stage. The glycine placement in highly curved strand positions
and the selection of pocket positions are managed by DesignBeNTF2.py (or DesignBeNTF2_test1.py at BeNTF2seq/design_with_PSSM to design using PSSMs). Pocket positions are selected by placing a virtual atom in the midpoint between the H3-S3 connection and
the S6 bulge, and choosing all positions whose C<sub>α</sub>-C<sub>β</sub> vector is pointing towards the virtual atom (the V<sub>atom</sub>-C<sub>α</sub>-C<sub>β</sub> angle is smaller than 90º), and their C<sub>α</sub> is closer than 8Å.

# Dependencies
pyrosetta\*<br/>
pandas<br/>
h5py<br/>
pickle<br/>
random<br/>
math<br/>
numpy<br/>
glob<br/>
json<br/>
sys<br/>
argparse<br/>
string<br/>
os<br/>
copy<br/>
<br/>
\*pyrosetta is free (with a subscription) for academic use: http://www.pyrosetta.org/dow


