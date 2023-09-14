- protocol_QSproteome_single_uniprot.pl
Calls functions from general.pm and functions_get_contacts.pm. Invokes process_and_analyze_AF_model.R for processing and analyzing AF models. Invokes select_best_rmsd_clashes_byfile.R to select the best RMSD clashes


- general.pm
Contains utility functions for contacts calculation.


- functions_get_contacts.pm
Contains functions to get contacts, which are used in protocol_QSproteome_single_uniprot.pl


- launch_reconstruct_ananas_cN.sh
Script that launches a reconstruction process using AnAnaS and Molprobity.


- process_and_analyze_AF_model.R
An R script that is called by protocol_QSproteome_single_uniprot.pl to process and analyze AF models. Calculate and remove the disordered regions, and all the PAE/probability scores.


- select_best_rmsd_clashes_byfile.R
An R script that is called by protocol_QSproteome_single_uniprot.pl to select the best symmetry based on RMSD and clash scores.
