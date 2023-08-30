#FIND_DISO:
#	cd bin/002_update_residue/ ;\
#	perl 0000_compute_disorder_EXE.pl $(SPECIES) 20 ;\
#	cd .. ;\
#

my $PDBFILE = $ARGV[0];

require("/media/elusers/users/hugo/15_alphafold/37_revision_Cell/functions_get_contacts.pm");

### STEP 1
## Calculate the contacts in the model
my @PDBFILES = ($PDBFILE);             

process_contacts(@PDBFILES, "test", 1);

### STEP 2
## Calculate nodiso1, 2, 3 
## Calculate clashes
my $CONTACTFILE = 
system("Rscript 01_remove_low_plddt.R $PDBFILE $CONTACTFILE")

### STEP 3
## Calculate PAE1, 2, 3
##

##CONTACT1:
##	cd bin/005_contacts ;\
##	perl 001_process_contacts_EXE.pl $(SPECIES) 40 ;\
##	cd .. ;\
#
#
#
##CONTACT2:
##	cd bin/005_contacts ;\
##	perl 02_concat_contact.pl $(SPECIES) ;\
##	cd .. ;\
#
#
#
#####
##### ==> this was supposed to be done from a file, so check the perl script and maybe change it.
#####
###   That is not working for permission reasons ==> I've got to copy the file and then run it in mysqlQueryBrowser logged in as root
###
###   1. sudo cp ../../results/05_contacts/all/contacts-ct-2022-11-10_all.txt /var/lib/mysql-files/
###   2. ALTER TABLE res_contact DISABLE KEYS
###   3. LOAD DATA INFILE '/var/lib/mysql-files/contacts-ct-2022-11-10_all.txt' INTO TABLE res_contact FIELDS TERMINATED BY ' '
###   4. ALTER TABLE res_contact ENABLE KEYS
###
### CONTACT3:
###	cd bin/005_contacts ;\
###	that's not working for permission issues -- perl 02_add_res_inter_FAST.pl $(SPECIES) ;\
###	cd .. ;\
#
#
####
#### Update protein length, clashes, and pLDDT of interface residues
####
#### ==> UPDATE complex SET length='$protSize', clash_ct='$Nclash_ct', clash_ct_int='$Nclash_ct_int', clash_res='$Nclash_res', clash_res_int='$Nclash_res_int', pLDDTint='$pLDDTint' WHERE code ='$CODE_LONG'
####
#### ---> CONVERT TO R AND READ THE CONTACT FILE, PRINT ERROR MESSAGE IF SOMETHING IS WRONG
####
#CLASH_AND_CO:
#	cd bin/002_update_residue/ ;\
#	perl 003_add_complex_info_EXE.pl $(SPECIES) 10 ;\
#	cd .. ;\
#
####
#### Calculates the pae_cplx 1/2/3/4
####
#### UPDATE complex SET pae_mono =",MONO.mean,", pae_mono3 =",MONO.3.mean,", pae_cplx=",
###           DIM.mean,", ppae_cplx1=",PAE1,", pae_cplx2=",PAE2,", pae_cplx3=",PAE3,", pae_cplx4=",ct.score2,
###          ", n_res_intf=",n_res_in_contact,
###          ", n_con_intf_diso1=",sum(res_in_contact_diso1),
###          ", n_con_intf_diso2=",sum(res_in_contact_diso2),
###          ", n_con_intf_diso3=",sum(res_in_contact_diso3),
###          ", n_con_intf=",nrow(res_in_contact)," WHERE code = '",CODE,"
### YOU CAN LOAD EVERYTHING IN R AGAIN AND GENERATE THE DIFFERENT PAEs
#GET_PAE:
#	cd bin/002_update_residue/ ;\
#	perl 01_compute_pLDTT_EXE.pl $(SPECIES) 40 ;\
#	cd .. ;\
#
### rewrites PDB files into the nodiso75 and nodiso80 folders
#FILTER_DISO:
#	cd bin/002_update_residue/ ;\
#	perl 02_filter_diso_EXE.pl $(SPECIES) 20 ;\
#	cd .. ;\
#
####
#### Benchmark stuff must be done first, dimerProbabilities are used in the selection of representatives, so this should be done before
####
#DIMER_PROBA:
#	cd bin/003_extract_scores/ ;\
#	perl 004_add_infoDB2DB_DimerProba_EXE.pl $(SPECIES) 40 ;\
#	cd .. ;\
#
### Runs representatives
#AF_REPRE_SUPERPOSE1:
#	cd bin/010_align_af2pdb/ ;\
#	perl 01_align_af2self_EDL_EXE.pl 1 $(SPECIES) 40 ;\
#	cd .. ;\
#
#AF_REPRE_SUPERPOSE2:
#	cd bin/010_align_af2pdb/ ;\
#	perl 01_align_af2self_EDL_EXE.pl 2 $(SPECIES) 40 ;\
#	cd .. ;\
#
#AF_REPRE_SUPERPOSE3:
#	cd bin/010_align_af2pdb/ ;\
#	perl 01_align_af2self_EDL_EXE.pl 3 $(SPECIES) 40 ;\
#	cd .. ;\
#
####
#### AF repre count 
#### Analyses superpositions to count representatives
#### //!\\ here I should use the latest trick implemented in the ref2ref clustering to treat Cn(n>=3) correctly
#AF_REPRE_N:
#	cd bin/003_extract_scores/ ;\
#	perl 002_add_infoDB2DB_nAfRepre_EXE.pl  $(SPECIES) 40 ;\
#	cd .. ;\
#
####
#### AF REPRE -- ASSIGNS REPRESENTATIVE STRUCTURES
#### I need the probas to select the best representative
#### So I do it later
#AF_REPRE_SET:
#	cd bin/010_align_af2pdb/ ;\
#	perl 02_get_af_repre_structure_EDL_EXE.pl 1 $(SPECIES) 20 ;\
#	perl 02_get_af_repre_structure_EDL_EXE.pl 2 $(SPECIES) 20 ;\
#	perl 02_get_af_repre_structure_EDL_EXE.pl 3 $(SPECIES) 20 ;\
#	cd .. ;\
#
#AF2PDB:
#	cd bin/03_13_struct_ali ;\
#	perl 001_struct_ali_EXE.pl $(SPECIES) 10 ;\
#	cd .. ;\
#
#### QUERY FIRST, TARGET NEXT
#AF2AF:
#	cd bin/03_13_struct_ali ;\
#	perl 002_struct_ali_af2af_EXE.pl pf ec 10 ;\
#	cd .. ;\
#
#
#RUN_BENCHMARK:
#	cd bin/03_13_struct_ali ;\
#	perl 003_struct_ali_bench_EXE.pl $(SPECIES) 10 ;\
#	cd .. ;\
#
#REFSET_SUPP:
#	cd bin/03_13_struct_ali ;\
#	perl 005_struct_ali_bench_inspectREFSET_EXE.pl ../010_analyze/hs_REF_v3.csv 10 ;\
#
#
#
###
### Goes over the structures from AF2 prediction for yeast
### --> copies PDB files to folder /data5/elevy/01_3dcomplexV0/data/BU_all_renum
### --> Add entries to complex table
#start_update_sc:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/YEAST_QS_2/ 2 V1 sc AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_ec:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/EC_QS_2/ 2 V1 ec AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_hs:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/HS_QS_2/ 2 V1 hs AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_pd:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/PD_QS_2/ 2 V1 pd AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_hs_mono:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/HS_QS_1/ 1 V0 hm AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#### TBD
#start_update_bs:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/BS_QS_2/ 2 V1 bs AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_ca:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/CA_QS_2/ 2 V1 ca AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_mm:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data4/01_AFDB/MM_QS_2/ 2 V1 mm AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#start_update_h5:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data5/01_AFDB/HS50_QS_2/ 2 V2 h5 AF_round1_cov70_id50_mmseq2_hugo ;\
#	cd .. ;\
#
#start_update_pf:
#	cd bin/001_UPDATE ;\
#	perl 01_ADD_SPECIES.pl /data5/01_AFDB/PF_QS_2/ 2 V1 pf AF_round1_cov70_id50_mmseq2_collabfold ;\
#	cd .. ;\
#
#
#### This will process the FASTA file and populate the proteome table
#### 
#### Unlike the chain table, each sequence is present once only and all sequences will be there even those witohut AF prediction.
#fill_proteome:
#	cd bin/002_update_chain/ ;\
#	perl 02_add_proteome_EXE.pl sc hs ec tt ca at mm bs ;\
#	cd .. ;\
#
#
#
#### Will run FASTA on all proteomes as well
#### -> concatenate all proteomes into ALL_FASTA and runs FASTA on 40 bunches
####
#
#### Will run STRUCT SUPERPOS on all proteomes as well
#### Since the complexes are fairly small (dimers it OK to run 10 jobs in parallel as memory shouldnt be an issue)
####
#
#### ASA --> we should use freesasa
####
####
#
#### Contacts
####
####
#CONTACTS_ec:
#	cd bin/005_contacts ;\
#	perl 001_process_contacts_EXE.pl ec 20 ;\
#	cd .. ;\
#
#CONTACTS_sc:
#	cd bin/005_contacts ;\
#	perl 001_process_contacts_EXE.pl sc 20 ;\
#	cd .. ;\
#
#CONTACTS_hs:
#	cd bin/005_contacts ;\
#	perl 001_process_contacts_EXE.pl hs 20 ;\
#	cd .. ;\
#
#CONTACTS_pd:
#	cd bin/005_contacts ;\
#	perl 001_process_contacts_EXE.pl pd 20 ;\
#	cd .. ;\
#
#CONTACT2_ec:
#	cd bin/005_contacts ;\
#	perl 02_concat_contact.pl ec ;\
#	perl 02_add_res_inter.pl ec ;\
#	cd .. ;\
#
#CONTACT2_sc:
#	cd bin/005_contacts ;\
#	perl 02_concat_contact.pl sc ;\
#	perl 02_add_res_inter.pl sc ;\
#	cd .. ;\
#
#CONTACT2_hs:
#	cd bin/005_contacts ;\
#	perl 02_concat_contact.pl hs ;\
#	perl 02_add_res_inter.pl hs ;\
#	cd .. ;\
#
#CONTACT2_pd:
#	cd bin/005_contacts ;\
#	perl 02_concat_contact.pl pd ;\
#	perl 02_add_res_inter.pl pd ;\
#	cd .. ;\
#
#CONTACT3_sc:
#	cd bin/005_contacts ;\
#	perl 002_add_complex_info_EXE.pl sc 10 ;\
#	cd .. ;\
#
#CONTACT3_ec:
#	cd bin/005_contacts ;\
#	perl 002_add_complex_info_EXE.pl ec 10 ;\
#	cd .. ;\
#
#CONTACT3_hs:
#	cd bin/005_contacts ;\
#	perl 002_add_complex_info_EXE.pl hs 10 ;\
#	cd .. ;\
#
#CONTACT3_pd:
#	cd bin/005_contacts ;\
#	perl 002_add_complex_info_EXE.pl pd 10 ;\
#	cd .. ;\
#
#RUN_BENCHMARK:
#	cd bin/03_13_struct_ali ;\
#	perl 003_struct_ali_bench_EXE.pl pd 10 ;\
#	cd .. ;\
#
