rm(list=ls())

suppressMessages(library(RMySQL))
suppressMessages(library(igraph))
suppressMessages(library(stringr))
suppressMessages(library(bio3d))
suppressMessages(library(rjson))

# args = commandArgs(trailingOnly = TRUE)

# if (length(args)!=2) {
#   stop("One arguments expected --- USAGE: Rscript 0000_remove_disorder.R CODE CONTACTS", call.=FALSE)
# }

# CODE = args[1] 
# JSONFILE = args[2]
# CONTACTS = args[3]
# OUTPATH = args[4]

# CODE = "/media/elusers/users/hugo/15_alphafold/37_revision_Cell/Q8WV44_V1_5.pdb"
CODE = "Q8WV44_V1_5"
JSONFILE = "data/Q8WV44_rank_2_model_5_ptm_seed_0_pae.json.bz2"
CONTACTS = read.table("data/Q8WV44_V1_5_FULL_CONTACTS.txt")
OUTPATH = "outpath/"

colnames(CONTACTS) = c("code", "chain1", "chain2", "res1", "res2", 
                       "resname1", "resname2", "n_contacts", "dmin", "dmax", "davg")
CONTACTS$resid1 = paste(CONTACTS$chain1, CONTACTS$res1, sep = "_")
CONTACTS$resid2 = paste(CONTACTS$chain2, CONTACTS$res2, sep = "_")
CONTACTSALL = CONTACTS # intra and inter-chains contacts
CONTACTS = CONTACTS[CONTACTS$chain1 != CONTACTS$chain2,] # we only check inter-chains contacts

Nclash_ct = sum(CONTACTS$dmin < 2)
Nclash_ct_int = sum(CONTACTS$dmin < 2 & (CONTACTS$chain1 != CONTACTS$chain2))
Nclash_res = length(unique(CONTACTS$res1[CONTACTS$dmin < 2]))
Nclash_res_int = length(unique(CONTACTS$res1[CONTACTS$dmin < 2 & (CONTACTS$chain1 != CONTACTS$chain2)]))


#################### READ PDB FILE ########################

pdb_file <- "/media/elusers/users/hugo/15_alphafold/37_revision_Cell/Q8WV44_V1_5.pdb"
pdb_data <- read.pdb(pdb_file)
pdb_dataframe <- as.data.frame(pdb_data$atom)
pdb_dataframe_plddt = pdb_dataframe[!duplicated(pdb_dataframe[,c("chain","resno","b")]), c("chain","resno","b")]
colnames(pdb_dataframe_plddt) = c("chain", "resnum", "bfact")


#################### Part getting contact matrix ########################
indices = paste(pdb_dataframe_plddt$chain, pdb_dataframe_plddt$resnum, sep = "_")
N = nrow(pdb_dataframe_plddt)
CTmat = matrix(ncol = N, nrow = N, data = 0)
colnames(CTmat) = indices
rownames(CTmat) = indices

# populating the matrix with contacts
CTmat[as.matrix(CONTACTSALL[,c("resid1", "resid2")])] = 1

####################### Part subsetting pLDDT ###########################
index.nodiso = pdb_dataframe_plddt$bfact > 40

if(sum(index.nodiso)>5){
  nodiso.median = median(pdb_dataframe_plddt$bfact[index.nodiso])
} else {
  nodiso.median = 40
}

index.nodiso75 = (pdb_dataframe_plddt$bfact > 75) | (pdb_dataframe_plddt$bfact > nodiso.median)
index.removed75 = which(! index.nodiso75)
CTmat[index.removed75, ] = 0
CTmat[, index.removed75] = 0

# ### Now transforms the CT matrix into a graph (see my R course from last year)
CTgraph = graph_from_adjacency_matrix(CTmat, mode = "undirected")

### Get the residue indexes of the largest component
all.comp = clusters(CTgraph)
largest.comp = all.comp$membership[all.comp$membership == which.max(all.comp$csize)]
res.to.keep = names(largest.comp)
res.to.update= as.data.frame(t(data.frame(str_split(res.to.keep, "_"))))
res.to.update$nodiso3 = TRUE
row.names(res.to.update) = NULL
colnames(res.to.update) = c("chain", "resnum", "nodiso3")

data.diso = data.frame(chain = pdb_dataframe_plddt$chain,
                       resnum = pdb_dataframe_plddt$resnum,
                       nodiso1 = index.nodiso,
                       nodiso2 = index.nodiso75)
data.all.diso = merge(data.diso, res.to.update, by = c("chain", "resnum"), all.x = TRUE)
data.all.diso[is.na(data.all.diso$nodiso3), "nodiso3"] = FALSE

################################## PART PAE #####################################

res_in_contact = CONTACTS[,c("res1", "res2", "chain1", "chain2")]

res_diso1 = data.all.diso[data.all.diso$nodiso1,c("resnum", "chain")]
res_diso2 = data.all.diso[data.all.diso$nodiso2,c("resnum", "chain")]
res_diso3 = data.all.diso[data.all.diso$nodiso3,c("resnum", "chain")]

isA = res_diso1[,2]=="A"
res.diso1.A = paste0(res_diso1[isA,1],res_diso1[isA,2])
res.diso1.B = paste0(res_diso1[!isA,1],res_diso1[!isA,2])

isA = res_diso2[,2]=="A"
res.diso2.A = paste0(res_diso2[isA,1],res_diso2[isA,2])
res.diso2.B = paste0(res_diso2[!isA,1],res_diso2[!isA,2])

isA = res_diso3[,2]=="A"
res.diso3.A = paste0(res_diso3[isA,1],res_diso3[isA,2])
res.diso3.B = paste0(res_diso3[!isA,1],res_diso3[!isA,2])

n_res_in_contact = length( unique(c(res_in_contact[,1],res_in_contact[,2])))

res_in_contact_diso1 = 
  CONTACTS$res1 %in% res_diso1$resnum[res_diso1$chain==CONTACTS$chain1[1]] &
  CONTACTS$res2 %in% res_diso1$resnum[res_diso1$chain==CONTACTS$chain2[1]]

res_in_contact_diso2 = 
  CONTACTS$res1 %in% res_diso2$resnum[res_diso2$chain==CONTACTS$chain1[1]] &
  CONTACTS$res2 %in% res_diso2$resnum[res_diso2$chain==CONTACTS$chain2[1]]

res_in_contact_diso3 = 
  CONTACTS$res1 %in% res_diso3$resnum[res_diso3$chain==CONTACTS$chain1[1]] &
  CONTACTS$res2 %in% res_diso3$resnum[res_diso3$chain==CONTACTS$chain2[1]]

n_contact = nrow(CONTACTS)
if(is.na(n_contact)) {
  n_contact=0
}



## Read the json file
if(file.exists(JSONFILE)){
  
  tmp2 = fromJSON(file=JSONFILE)
  mat = matrix(ncol=length(tmp2[[1]]$residue1)^0.5,nrow=length(tmp2[[1]]$residue1)^0.5,data=tmp2[[1]]$distance, byrow=TRUE)
  L = dim(mat)[1]    
  PROT.L = L/2
  chA = paste0(1:PROT.L,"A")
  chB = paste0(1:PROT.L,"B")
  colnames(mat)=c(chA,chB)
  rownames(mat)=c(chA,chB)
  
  ###
  ### This is to define pae_cplx
  ###
  idx = floor(seq(from=1,to=L,len=20))
  x1.1 = idx[2]
  x1.2 = x1.1 + (idx[8]-idx[2])
  x2.1 = idx[12]
  x2.2 = x2.1 + (idx[8]-idx[2])
  
  ### Monomer  score
  MONO1 = mean( mat[ c(x1.1:x1.2) , c(x1.1:x1.2) ] )
  MONO2 = mean( mat[ c(x2.1:x2.2) , c(x2.1:x2.2) ] )
  MONO.mean = round(mean(c(MONO1,MONO2)),2)
  
  MONO1.3 = mean( mat[ res.diso3.A , res.diso3.A  ] )
  MONO2.3 = mean( mat[ res.diso3.B , res.diso3.B  ] )
  MONO.3.mean = round(mean(c(MONO1.3,MONO2.3), na.rm=TRUE),2)
  if(is.na(MONO.3.mean)){
    MONO.3.mean = 100
  }
  ### Dimer  score
  DIM1 = mean( mat[ c(x1.1:x1.2) , c(x2.1:x2.2) ] )
  DIM2 = mean( mat[ c(x2.1:x2.2) , c(x1.1:x1.2) ] )
  DIM.mean = round(mean(c(DIM1, DIM2)),2)
  
  #################################
  ### pae_cplx1
  if(length(res_diso1)>0 && length(unique(res_diso1[,2]))>1 ){
    PAE1 = round(
      ((mean( mat[ res.diso1.A , res.diso1.B ]) + mean( mat[ res.diso1.B , res.diso1.A ]))/2)
      ,1)
  } else {
    PAE1 = 50
  }
  
  #################################
  ### pae_cplx2
  ##
  if(length(res_diso2)>0 && length(unique(res_diso2[,2]))>1){
    PAE2 = round(
      ((mean( mat[ res.diso2.A , res.diso2.B ]) + mean( mat[ res.diso2.B , res.diso2.A ]))/2)
      ,1 )
  } else {
    PAE2 = 50
  }
  #################################
  ### pae_cplx3
  if(length(res_diso3)>0 && length(unique(res_diso3[,2]))>1){
    PAE3 = round(
      ((mean( mat[ res.diso3.A , res.diso3.B ]) + mean( mat[ res.diso3.B , res.diso3.A ]))/2)
      ,1 )
  } else {
    PAE3 = 50
  }
  #################################
  ### pae_cplx4
  ct.score2 = 50
  if(nrow(CONTACTS)>5){

    contacts.A = paste0(CONTACTS[,"res1"],CONTACTS[,"chain1"])
    contacts.B = paste0(CONTACTS[,"res2"],CONTACTS[,"chain2"])
    
    ct.score2 = round(
      (mean(mat[ contacts.A, contacts.B ])+mean(mat[ contacts.B, contacts.A ]))/2
      ,2)
  }    
}


########################## CALCULATE DIMER PROBA ###############################

# data = dbGetQuery(mydb, paste0("SELECT pae_cplx2, pae_cplx3, pae_cplx4, af_repre1_N, af_repre2_N, n_con_intf_diso3 FROM complex WHERE code = '",CODE,"'"))
data = data.frame("pae_cplx2" = PAE2,
                  "pae_cplx3" = PAE3, 
                  "pae_cplx4" = ct.score2, 
                  "af_repre1_N" = 0, 
                  "af_repre2_N" = 0, 
                  "ndiso3" = sum(res_in_contact_diso3))
  
# repre = dbGetQuery(mydb, paste0("SELECT af_repre1_N FROM complex WHERE code = '",CODE,"'"))

if(is.na(data[5])){
  data[5] = 0
}

# colnames(data) = c("pae_cplx2","pae_cplx3", "pae_cplx4", "af_repre1_N", "af_repre2_N", "ndiso3")

my.logit.pae3      = readRDS( "/data5/elevy/01_3dcomplexV0/bin/010_analyze/logit_model_17_2_2022_FULL_pae3.RDS") 
my.logit.pae4      = readRDS("/data5/elevy/01_3dcomplexV0/bin/010_analyze/logit_model_17_2_2022_FULL_pae4.RDS")
my.logit.con3      = readRDS("/data5/elevy/01_3dcomplexV0/bin/010_analyze/logit_model_17_2_2022_FULL_con3.RDS")
my.logit.repre     = readRDS("/data5/elevy/01_3dcomplexV0/bin/010_analyze/logit_model_17_2_2022_FULL_repre.RDS")
my.logit.pae4.con3 = readRDS("/data5/elevy/01_3dcomplexV0/bin/010_analyze/logit_model_17_2_2022_FULL_pae4.con3.RDS")

#logodds      = predict(my.logit, newdata=data)
proba.pae3      = predict(my.logit.pae3, newdata=data, type="response")
proba.pae4      = predict(my.logit.pae4, newdata=data, type="response")
proba.con3      = predict(my.logit.con3, newdata=data, type="response")
proba.repre     = predict(my.logit.repre, newdata=data, type="response")
proba.pae4.con3 = predict(my.logit.pae4.con3, newdata=data, type="response")
proba.all   = max(c(proba.pae3, proba.pae4, proba.con3, proba.pae4.con3))

cat("dimer_proba_pae3` =",round(proba.pae3,5),", 
    `dimer_proba_pae4`=",round(proba.pae4,5), ", 
    `dimer_proba_con3`=",round(proba.con3,5),", 
    `dimer_proba_repre`=",round(proba.repre,5), ", 
    `dimer_proba_pae4_con3`=",round(proba.pae4.con3,5), ", 
    `dimer_proba_max`=",round(proba.all,5), "\n")

#### THEN we write everything:
## 1- pdb file model nodiso1/2/3


## 2- list info per residue -> nodiso mainly


## 3- Contacts file -> ok already done at previous step


## 4- A file with X column for all the PAE score etc
df_towrite = data.frame(PAE1 = PAE1,
                        PAE2 = PAE2,
                        PAE3 = PAE3,
                        PAE_interface = ct.score2,
                        dimer_proba = round(proba.all,5))

pdb_nodiso1 =
pdb_nodiso2 =
pdb_nodiso3 =

write.csv(df_towrite, paste0(OUTPATH, "/", CODE, "_probability_scores.csv"))