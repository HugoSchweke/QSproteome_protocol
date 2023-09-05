This directory contains all the output files after executing the first step of the QSproteome protocol.

It contains:
- Q8WV44_V1_5_nodiso1.pdb: example input pdb file (Q8WV44_V1_5.pdb) where residues with a pLDDT score below 40 are filtered out.
- Q8WV44_V1_5_nodiso2.pdb: Starting from Q8WV44_V1_5_nodiso1.pdb, a median pLDDT score was computed, and residues with a pLDDT score
below 75 and below the median value were discarded.
- Q8WV44_V1_5_nodiso3.pdb: Starting from Q8WV44_V1_5_nodiso2.pdb, we applied single linkage clustering on the contact matrix of the remaining residues and retained the largest cluster, thus eliminating disconnected
structural parts.
- Q8WV44_V1_5_FULL.txt:
- Q8WV44_V1_5_diso_info.csv:
- Q8WV44_V1_5_probability_scores.csv:
