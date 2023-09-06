# QSproteome


This repository contains the scripts and protocol from ([1](#ref-1)) 

It contains all the necessary scripts to calculate predict, detect and reconstruct its full size complex from an AlphaFold homodimeric model.

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

3- Optionnally, if the user has installed the AlphaFold Big Bang (?) , he can use this reconstructed complex as a template to model the full size complex.

</div>


# Installation
[Go to the top](#Table-of-contents)

### Requirements

QSPROTEOME is a tool that requires a UNIX-based OS system. It is written in perl (version 3.7), R (version 4.3.1) and bash. It and may optionally require AnAnas ([2](#ref-2)) and molprobity ([3](#ref-3)) if the user wants to reconstruct full homomeric complexes.

All those requirements (including APBS) are met in a [predefined Docker image](https://hub.docker.com/r/lopesi2bc/surfmap/tags) that we recommend the user to use. 

<details open>
<summary><b>For a usage of the docker image</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [Docker](https://docs.docker.com/get-docker/)

</details>

<details>
<summary><b>For a usage on your local OS</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [R >= 3.6](https://cran.r-project.org/)
- [APBS](https://github.com/Electrostatics/apbs/releases) (optional - only if you want to compute electrostatic potential)
 
</details>
<br>

> :bell: Please note that **whether you want to use the Docker image of SURFMAP or not, you will still need to [install the SURFMAP package](#How-to-install-SURFMAP)**. Indeed the package contains internal features that make the use of the Docker image totally transparent for the user who will not have to enter 'complex' commands for the connection of useful mounting points. In fact, the SURFMAP commands are almost exactly the same between the use of the docker image or not (see [here](#cmd_docker_or_not)).




## How to install QSPROTEOME
[Go to the top](#Table-of-contents)

First, make sure you meet the [system requirements](#requirements) outlined earlier and consider the [recommendation](#recommendation). Then, follow instructions described in option 1 or 2 if you're not interested in accessing/modifying the source code, otherwise prefer option 3. 

<a id="install_option1"></a>
<details open>
<summary><h4>Option 1: from the archive (git not required)</h4></summary>

First download an archive of our latest release <a href="https://github.com/i2bc/SURFMAP/releases/latest" target="_blank">here</a>.

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP v2.0.0
python3 -m pip install SURFMAP-2.0.0.zip # (or .tar.gz) 
```
</details>


# How it works
[Go to the top](#Table-of-contents)


### QSPROTEOME workflow: inputs/outputs

<div align="center">
  <img src="./doc/images/surfmap_workflow.png" width="70%"/>

<i>The figure above represents the main steps of the SURFMAP worflow to compute the projection on a 2D map of a protein surface feature. More details about each step can be found in our article: see the [published version](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01269) or its [free version](https://www.biorxiv.org/content/10.1101/2021.10.15.464543v1)</i>
</div>
<br>

QSPROTEOME needs in input an AlphaFold model of a homodimer in pdb format, as well as the associated json file provided by AF.
<br>
<br>

[Using a PDB file as input](#from-a-pdb-structure) is the most classic usage of SURFMAP. In this case, two outputs are generated: 
- the 2D map projection in a PDF format (PNG is also available)
- a matrix text file written in a SURFMAP-specific format

The matrix text file contains all information about each projected surface residue and their associated feature value. As the above figure shows, this text file is the direct input for the last step of the SURFMAP workflow as it is read to generate the 2D map projection.
<br>
<br>

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

# Usage of QSPROTEOME
[Go to the top](#Table-of-contents)

Once you have [installed the QSPROTEOME package](#how-to-install-qsproteome), you should be ready to use it. 

#### The example directory
To guide the user in the usage of QSPROTEOME, we will make use of files that you can find in the `example/` directory. You can see where this directory is located on your machine with the following command:

```bash
python3 -c "import surfmap; print(surfmap.PATH_TO_EXAMPLES)"
```

Please note that for all command examples illustrated below, we will make [use of the Docker image of SURFMAP](#use-surfmap-with-docker-or-not).


#### QSPROTEOME options

<details>
<summary>List of all QSPROTEOME arguments</summary>

<pre>usage: QSPROTEOME [-h] (-pdb PDB | -json JSON)

options:
  -h, --help        show this help message and exit
  -pdb PDB          Path to an AlphaFold model PDB file. *required*
  -json JSON        Path to the json file produced by AlphaFold along with the model. *required*
</pre>
</details>


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
