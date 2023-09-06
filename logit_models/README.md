This directory contains the logistic regression models stored as Rdata files. 

Here is a brief description of the files (for more detailed information, please refer to the manuscript (bioRxiv preprint doi: https://doi.org/10.1101/2023.06.09.544317)):

#### logit_model_FULL_pae3.RDS 

logistic regression model using the *PAE3* metric (mean inter-chain PAE values of interface residues in the core structure).

#### logit_model_FULL_pae4.RDS

logistic regression model using the *PAE4* metric (mean PAE values of interchain contacts in the core structures).

#### logit_model_FULL_con3.RDS

logistic regression model using the *CON3* metric (number of interchain residue-residue contacts in the core structures).

#### logit_model_FULL_repre.RDS

logistic regression model using the *repre* metric (the consistency of the five models as determined from the structural superposition of their core structure before the single linkage clustering was applied. Each model was assigned a score ranging from 0 (where no other model was structurally similar) to 8 when a model reciprocally matched all the other four models).

#### logit_model_FULL_pae4.con3.RDS

logistic regression model using a combination of the *PAE4* and *CON3* metrics.
