This directory contains all the output files after executing the first step of the QSproteome protocol.

It contains:

#### Q8WV44_V1_5_nodiso1.pdb: 

This is the example input pdb file (Q8WV44_V1_5.pdb) where residues with a pLDDT score below 40 are filtered out.

#### Q8WV44_V1_5_nodiso2.pdb: 

Starting from Q8WV44_V1_5_nodiso1.pdb, a median pLDDT score was computed, and residues with a pLDDT score
below 75 and below the median value were discarded.

#### Q8WV44_V1_5_nodiso3.pdb: 

Starting from Q8WV44_V1_5_nodiso2.pdb, a single linkage clustering was applied on the contact matrix of the remaining residues and the largest cluster was retained, thus eliminating disconnected structural parts.

#### Q8WV44_V1_5_FULL.txt: 

table of intra and interchain contacts. It is organized as follow:
- code = code of the pdb file in input
- chain1 = chain id of the first residue
- chain2 = chain id of the second residue
- res1 = residue number of the first residue
- res2 = residue number of the second residue
- rescode1 = one-letter code of the first residue
- rescode2 = one-letter code of the second residue
- dmin = distance 1 (?)
- dmax = distance 2 (?)
- davg = distance 3 (?)
  
#### Q8WV44_V1_5_diso_info.csv:

csv file that contain five columns:
- chain = chain of the model
- resnum = residue number
- nodiso1 = TRUE if the residue is present in the structure nodiso1, FALSE if filtered out
- nodiso2 = TRUE if the residue is present in the structure nodiso2, FALSE if filtered out
- nodiso3 = TRUE if the residue is present in the structure nodiso3, FALSE if filtered out

#### Q8WV44_V1_5_probability_scores.csv:

Contains five probability scores.
- PAE1 = 
- PAE2 = 
- PAE3 = 
- PAE_interface = probability of the model to be a physiological homomer based on a logistic regression model using the PAE4 metric (mean PAE values of interchain contacts in the core structures).
- dimer_proba = probability that the dimer is a physiological one.
