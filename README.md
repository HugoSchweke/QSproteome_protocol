# QSproteome


This repository contains the scripts and protocol from  ”An atlas of protein homo-oligomerization across domains of life.”, Schweke et al ([1](#ref-1)).

It contains all the necessary scripts to predict homomerization, detect symmetry and reconstruct full size complex from an AlphaFold homodimeric model.


# Table of contents

- [Aims](#Aims)
- [Installation](#Installation)
- [How it works](#How-it-works)
- [Usage of QSPROTEOME](#Usage-of-QSPROTEOME)
- [Supporting the project](#Supporting-the-project)
- [Contacts](#Contacts)
- [Licence](#Licence)
- [How to cite QSPROTEOME](#How-to-cite-QSPROTEOME)
- [References](#References)

# Aims
[Go to the top](#Table-of-contents)

QSPROTEOME is a protocol that detect homo-oligomerization from AlphaFold homodimeric predictions. The pipeline proposed here can:

1- Calculate a probability score that the homodimer given in input is a physiological one. 

2- Detect cyclic or trans symmetry and, if needed, reconstruct the full size complex with the AnAnas software ([2](#ref-2)).

Optionnally, if the user has installed the AlphaFold Big Bang (link to ColabFold git?) , it is possible to use this reconstructed complex as a template to model the full size complex.

</div>


# Installation
[Go to the top](#Table-of-contents)

### Requirements

QSPROTEOME is a tool that requires a UNIX-based OS system. It is written in perl (version XX), R (version 4.3.1) and bash. It may optionally require AnAnas ([2](#ref-2)) and molprobity ([3](#ref-3)) if the user wants to reconstruct full homomeric complexes.


## How to install QSPROTEOME
[Go to the top](#Table-of-contents)

- Install the following packages in R: igraph, stringr, bio3d, rjson, sys.

- Clone the git repository to your machine
```bash
# clone QSPROTEONE on your machine
git clone [https://github.com/HugoSchweke/QSproteome_protocol]
```

To enable symmetry detection and full size complexes reconstruction, AnAnaS and the Phenix software suite must be installed as follow: 

- Install the Phenix software suite. You can download it here: [Link](https://phenix-online.org/download/)

- Export the path to phenix in your bashrc:
```bash
cat "export PHENIX_CLASHSCORE=PATHTOPHENIX/build/bin/phenix.clashscore" >> ~/.bashrc
source ~/.bashrc
```

- Install the AnAnaS software. You can download it here: [Link](https://team.inria.fr/nano-d/software/ananas/)

- export the path to AnAnaS in your bashrc:
```bash
cat "export ANANAS=PATHTOANANAS/bin/ananas" >> ~/.bashrc
source ~/.bashrc
```


# How it works
[Go to the top](#Table-of-contents)


### QSPROTEOME workflow: inputs/outputs

QSPROTEOME needs in input an AlphaFold model of a homodimer in pdb format, as well as the associated json file provided by AF.
<br>

Six outputs are generated, plus two additional output if the user requested the reconstruction of the full size complex: 

**Main outputs**
- A pdb file where residues are filtered out according of to the *nodiso1* definition (residues with a pLDDT score below 40 are filtered out)
- A pdb file where residues are filtered out according of to the *nodiso2* definition (starting from the nodiso1 file, a median pLDDT score is computed, and residues with a pLDDT score
below 75 and below the median value are discarded.)
- A pdb file where residues are filtered out according of to the *nodiso3* definition (Starting from the nodiso2 file, a single linkage clustering is applied on the contact matrix of the remaining residues and the largest cluster is retained, thus eliminating disconnected structural parts.)
- a csv file that indicating which residues are filtered out following the *nodiso* definitions
  
<details>
<summary>Example of a table of disorder format (.csv)</summary>
 
<pre> 
chain,resnum,nodiso1,nodiso2,nodiso3
A,1,FALSE,FALSE,FALSE
A,2,FALSE,FALSE,FALSE
A,3,FALSE,FALSE,FALSE
A,4,FALSE,FALSE,FALSE
A,14,TRUE,FALSE,FALSE
A,15,TRUE,FALSE,FALSE
A,16,TRUE,FALSE,FALSE
A,17,TRUE,TRUE,FALSE
A,18,TRUE,TRUE,FALSE
A,19,TRUE,TRUE,FALSE

- chain = chain of the model
- resnum = residue number
- nodiso1 = TRUE if the residue is present in the structure nodiso1, FALSE if filtered out
- nodiso2 = TRUE if the residue is present in the structure nodiso2, FALSE if filtered out
- nodiso3 = TRUE if the residue is present in the structure nodiso3, FALSE if filtered out
 </pre>
</details>
 
- a contact file containing information regarding all the residues in contact in the input pdb file.

<details>
<summary>Example of a table of contacts format (.txt)</summary>

<pre> 
code chain1 chain2 res1 res2 rescode1 rescode2 d1 d2 d3
Q8WV44_V1_5 B B 9 13 N T 2 2.806 3.287 3.046
Q8WV44_V1_5 B B 9 12 N Q 1 3.172 3.172 3.172
Q8WV44_V1_5 B B 10 13 P T 3 3.015 3.243 3.125
Q8WV44_V1_5 B B 10 14 P L 5 3.147 4.014 3.480
Q8WV44_V1_5 B B 11 15 V Q 2 3.128 3.553 3.341
Q8WV44_V1_5 B B 11 14 V L 4 3.182 3.670 3.502
Q8WV44_V1_5 B B 12 15 Q Q 2 3.354 3.796 3.575
Q8WV44_V1_5 B B 12 9 Q N 1 3.172 3.172 3.172
Q8WV44_V1_5 B B 12 16 Q E 3 3.177 3.760 3.451
Q8WV44_V1_5 B B 13 16 T E 1 3.502 3.502 3.502

- code = code of the pdb file in input
- chain1 = chain id of the first residue
- chain2 = chain id of the second residue
- res1 = residue number of the first residue
- res2 = residue number of the second residue
- rescode1 = one-letter code of the first residue
- rescode2 = one-letter code of the second residue
- dmin = distance 1
- dmax = distance 2
- davg = distance 3
 </pre>
</details>


- a probability scores file containing scores such as interaction probability (i.e. the probability that the interaction modeled by AlphaFold is a physiological one).

<details>
<summary>Example of a probability scores file (.csv)</summary>

<pre> 
PAE1,PAE2,PAE3,PAE_interface,dimer_proba
10.8,6.7,4,4.68,0.97528


- PAE1 = PAE score of the nodiso1 residues
- PAE2 = PAE score of the nodiso2 residues
- PAE3 = PAE score of the nodiso3 residues
- PAE_interface = PAE score of the interface residues
- dimer_proba = Probability of the interaction to be a physiological one.
 </pre>
</details>

**Optional outputs**


# Usage of QSPROTEOME
[Go to the top](#Table-of-contents)

Once you have [installed the QSPROTEOME package](#how-to-install-qsproteome), you should be ready to use it. 

#### QSPROTEOME options

<details>
<summary>List of all QSPROTEOME arguments</summary>

<pre>usage: QSPROTEOME [-h] (-pdb PDB | -json JSON)

options:
  -h, --help       show this help message and exit
  --pdb            Path to an AlphaFold model PDB file. *required*
  --json           Path to the json file produced by AlphaFold along with the model. *required*
  --outpath        Path where all output files will be written. *required*
  --reconstruct    Boolean. If specified, the script will invoke AnAnas to reconstruct the full size complex (if needed), starting from the cropped models (nodiso3)
</pre>
</details>


#### The example directory
To guide the user in the usage of QSPROTEOME, we will make use of files that you can find in the `example/` directory. To execute the script using the pdb file 

```bash
cd scripts
perl protocol_QSproteome_single_uniprot.pl --pdb ../example/P32907_V1_1.pdb --json ../example/P32907_rank_1_model_1_ptm_seed_0_pae.json.bz2 --outpath ../../test
```
With this command the script protocol_QSproteome_single_uniprot.pl will calculate the interaction probability and detect the best possible symmetry of the pdb file P32907_V1_1.pdb using structural information from the pdb file and information contained in the json file P32907_rank_1_model_1_ptm_seed_0_pae.json.bz2. It will write all the results in the directory ../../test. 

Now, let's have a look at the results. In the probability file:
<details>
<summary>P32907_V1_1_probability_scores.csv (.csv)</summary>
<pre> 
PAE1,PAE2,PAE3,PAE_interface,dimer_proba
10.8,6.7,4,4.68,0.97528
 </pre>
</details>

We can see that the dimer probability (column dimer_proba) is 0.97528. The interaction predicted by AlphaFold is thus most likely a physiological one, thus *P32907_V1_1* forms a homomer. 
Here we did not specified the --reconstruct option, so we have no information about the symmetry of this homomer, and the full size complex is not reconstructed. 

We can get these informations using the --reconstruct option:

```bash
perl protocol_QSproteome_single_uniprot.pl --pdb ../example/P32907_V1_1.pdb --json ../example/P32907_rank_1_model_1_ptm_seed_0_pae.json.bz2 --outpath ../../test --reconstruct
```
This command is similar to the previous one, but the script will also detect the best matching cyclic symmetry, and, if needed (i.e., if a symmetry superior to C2 is detected), reconstruct the full cyclic complex based on the nodiso3 pdb file.

<br>

Let's look at *P32907_V1_1_nodiso3_all_csym.dat*, the file containing the result of the symmetry detection with AnAnaS:

<details>
<summary>P32907_V1_1_nodiso3_all_csym.dat (.dat)</summary>
<pre> 
symmetry av.rmsd clashscore
c2 18.620165 NA
c3 11.111322 NA
c4 6.065964 128.54
c5 2.808240 85.27
c6 0.654587 44.23
c7 1.110083 47.91
c8 2.284813 61.57
c9 3.210861 68.65
c10 2.808240 NA
c11 1.610034 607.33
c12 0.654587 NA
 </pre>
</details>


We can see that, according to AnAnaS, c6 symmetry has the lowest rmsd, as well as the lowest clashscore (please take note that there are no clashscore values for c2 and c3 because the rmsd is too high (>7A) and for c10 and c12 because those have exactly the same rmsd than c5 and c6, meaning that in that case they are just two c5 or c6 superposed, respectively). 
<br>
P32907_V1_1 forms a homohexamer. The complex reconstructed with AnAnaS using the trimmed file P32907_V1_1_nodiso3.pdb can be found in the output directory under the name *P32907_V1_1_nodiso3_c6.pdb*.

Please take not that using such a trimmed pdb file  is important, as a reconstructed complex using the full length model can lead to a lot of clashes, due to the low plddt flexible regions. 

<br>



If we take the model P25298_V1_5.pdb 

```bash
perl protocol_QSproteome_single_uniprot.pl --pdb ../example/P25298_V1_5.pdb --json ../example/P25298_rank_1_model_5_ptm_seed_0_pae.json.bz2 --outpath ../../test --reconstruct
```

Now let's have a look at the symmetry file:

<details>
<summary>P32907_V1_1_nodiso3_all_csym.dat (.dat)</summary>
<pre> 
symmetry av.rmsd clashscore
c2 0.125649 NA
 </pre>
</details>

Here the AlphaFold model is homodimer of C2 symmetry. In that case, no reconstruction using AnAnaS is necessary, as the full complex involves only two subunits.

# Supporting the project
[Go to the top](#Table-of-contents)

- If you find a bug or have a suggestion for a new feature, please report it via an [issue](https://github.com/HugoSchweke/QSproteome_protocol/issues)
- If you find QSPROTEOME useful, consider starring the repository


# Contacts
[Go to the top](#Table-of-contents)

If you have any question regarding QSPROTEOME, you can contact us:
- [@emmanuel.levy@weizmann.ac.il](mailto:@emmanuel.levy@weizmann.ac.il) (project leader and original code author)
- [@hugo.schweke@weizmann.ac.il](mailto:hugo.schweke@weizmann.ac.il) (original code author)


# Licence
[Go to the top](#Table-of-contents)

This project is under the MIT License terms. Please have a look at the LICENSE file for more details.


# How to cite QSPROTEOME
[Go to the top](#Table-of-contents)

If QSPROTEOME has been useful to your research, please cite us:

> An atlas of protein homo-oligomerization across domains of life
Hugo Schweke, Tal Levin, Martin Pacesa, Casper A. Goverde, Prasun Kumar, Yoan Duhoo, Lars J. Dornfeld, Benjamin Dubreuil, Sandrine Georgeon, Sergey Ovchinnikov, Derek N. Woolfson, Bruno E. Correia, Sucharita Dey, Emmanuel D. Levy
bioRxiv 2023.06.09.544317; doi: https://doi.org/10.1101/2023.06.09.544317 [Link](https://www.biorxiv.org/content/10.1101/2023.06.09.544317v1)


Moreover, if you use the pipeline of the homomer structure prediction in your research, please cite the following papers:
<br>
> Pagès, Guillaume, Elvira Kinzina, and Sergei Grudinin. 2018. “Analytical Symmetry Detection in Protein Assemblies. I. Cyclic Symmetries.” Journal of Structural Biology 203 (2): 142–48. [Link](https://doi.org/10.1016/j.jsb.2018.04.004)

> Williams, Christopher J., Jeffrey J. Headd, Nigel W. Moriarty, Michael G. Prisant, Lizbeth L. Videau, Lindsay N. Deis, Vishal Verma, et al. 2018. “MolProbity: More and Better Reference Data for Improved All-Atom Structure Validation.” Protein Science: A Publication of the Protein Society 27 (1): 293–315. [Link](https://doi.org/10.1107/S0907444909042073)

# References
[Go to the top](#Table-of-contents)

<a id="ref-1"></a>

> (1) Hugo Schweke, Tal Levin, Martin Pacesa, Casper A. Goverde, Prasun Kumar, Yoan Duhoo, Lars J. Dornfeld, Benjamin Dubreuil, Sandrine Georgeon, Sergey Ovchinnikov, Derek N. Woolfson, Bruno E. Correia, Sucharita Dey, Emmanuel D. Levy. ”An atlas of protein homo-oligomerization across domains of life.” bioRxiv 2023.06.09.544317.


<a id="ref-2"></a>

> (2) Pagès, Guillaume, Elvira Kinzina, and Sergei Grudinin. 2018. “Analytical Symmetry Detection in Protein Assemblies. I. Cyclic Symmetries.” Journal of Structural Biology 203 (2): 142–48.


<a id="ref-3"></a>

> (3) Christopher J. Williams, Jeffrey J. Headd, Nigel W. Moriarty, Michael G. Prisant, Lizbeth L. Videau, Lindsay N. Deis, Vishal Verma, et al. 2018. “MolProbity: More and Better Reference Data for Improved All-Atom Structure Validation.” Protein Science: A Publication of the Protein Society 27 (1): 293–315.
