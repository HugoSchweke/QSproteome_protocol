# QSproteome_protocol
This repository contains all the necessary script and data to calculate everything from input pdb model + json file to dimer probability.
The pipeline is based on perl and R.

Directories:
data_test: contains files for testing


Requirements:
- perl (version XX)
- R (version > 4.0.0)

 general.pm : set of perl functions necessary for the contact calculations part.
  functions_get_contacts.pm  : perl functions to calculate

# SURFMAP

<div align="center">
  <img src="./doc/images/toc_Schweke_SURFMAP_cmyk.png" width="80%"/>  
</div>

# Table of contents

- [Aims](#Aims)
- [Installation](#Installation)
- [How it works](#How-it-works)
- [Usage of SURFMAP](#Usage-of-SURFMAP)
- [Supporting the project](#Supporting-the-project)
- [Contacts](#Contacts)
- [Licence](#Licence)
- [How to cite SURFMAP](#How-to-cite-SURFMAP)
- [References](#References)

# Aims
[Go to the top](#Table-of-contents)

<div>
<img src="./doc/images/TOC_Schweke_manuscript_revisions_forGitHub.png" width="60%" align="right"/>

QSPROTEOME is a free standalone and easy-to-use command-line interface (CLI) software that enables the fast and automated 2D projection of either predefined features of protein surface (electrostatic potential, Kyte-Doolittle hydrophobicity, Wimley-White hydrophobicity, stickiness and surface relief) or any descriptor encoded in the temperature factor column of a PDB file. The 2D maps computed by QSPROTEOME can be used to analyze and/or compare protein surface properties.
</div>


# Installation
[Go to the top](#Table-of-contents)

### Requirements

QSPROTEOME is a CLI tool that requires a UNIX-based OS system. It is written in perl (version 3.7), R (version 4.3.1) and bash. It relies on the already included MSMS software ([1](#ref-1)) and may optionally require APBS ([2](#ref-2)) if the user wants to perform electrostatics calculations.

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

SURFMAP accepts as input either a *PDB file* or a *text file in a SURFMAP-specific matrix format*.
<br>
<br>

[Using a PDB file as input](#from-a-pdb-structure) is the most classic usage of SURFMAP. In this case, two outputs are generated: 
- the 2D map projection in a PDF format (PNG is also available)
- a matrix text file written in a SURFMAP-specific format

The matrix text file contains all information about each projected surface residue and their associated feature value. As the above figure shows, this text file is the direct input for the last step of the SURFMAP workflow as it is read to generate the 2D map projection.
<br>
<br>

[Using a text file in a SURFMAP-specific matrix format as input](#from-a-surfmap-matrix-file) represents a special case that could be useful if the user wants to generate a 2D map from an internally pre-processed matrix, such as to normalize or average with other matrices.

<details>
<summary>Example of a table of contacts format (.txt)</summary>

<pre> code chain1 chain2 res1 res2 rescode1 rescode2
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
</pre>
</details>



# Usage of QSPROTEOME
[Go to the top](#Table-of-contents)

Once you have [installed the QSPROTEOME package](#how-to-install-surfmap), you should be ready to use SURFMAP. 

#### The example directory
To guide the user in the usage of QSPROTEOME, we will make use of files that you can find in the `example/` directory of SURFMAP. You can see where this directory is located on your machine with the following command:

```bash
python3 -c "import surfmap; print(surfmap.PATH_TO_EXAMPLES)"
```

Please note that for all command examples illustrated below, we will make [use of the Docker image of SURFMAP](#use-surfmap-with-docker-or-not).


#### QSPROTEOME options

<details>
<summary>List of all SURFMAP options</summary>

<pre>usage: surfmap [-h] (-pdb PDB | -mat MAT | -v) -tomap TOMAP [-proj PROJ] [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth] [--png] [--keep]
               [--docker] [--pqr PQR] [-ff FF] [-verbose VERBOSE]

options:
  -h, --help        show this help message and exit
  -pdb PDB          Path to a PDB file
  -mat MAT          Input matrix. If the user gives an imput matrix, SURFMAP will directly compute a map from it.
  -v, --version     Print the current version of SURFMAP.
  -tomap TOMAP      Specific key of the feature to map. One of the following: stickiness, kyte_doolittle, wimley_white, electrostatics,
                    circular_variance, bfactor, binding_sites, all.
  -proj PROJ        Choice of the projection. Argument must be one of the following: flamsteed, mollweide, lambert. Defaults to flamsteed.
  -res RES          File containing a list of residues to map on the projection. Expected format has the following space/tab separated column
                    values: chainid resid resname
  -rad RAD          Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the
                    smoother the surface. Defaults to 3.0
  -d D              Output directory where all files will be written. Defaults to &apos;./output_SURFMAP_$pdb_$tomap&apos; with $pdb and $tomap based on
                    -pdb and -tomap given values
  -s S              Value defining the size of a grid cell. The value must be a multiple of 180. Defaults to 5.0
  --nosmooth        If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)
  --png             If chosen, a map in png format is computed (default: only pdf format is generated)
  --keep            If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)
  --docker          If chosen, SURFMAP will be run on a docker container (requires docker installed).
  --pqr PQR         Path to a PQR file used for electrostatics calculation. Option only available if &apos;-tomap electrosatics&apos; is requested.
                    Defaults to None.
  -ff FF            Force-field used by pdb2pqr for electrostatics calculation. One of the following: AMBER, CHARMM, PARSE, TYL06, PEOEPB,
                    SWANSON. Defaults to CHARMM.
  -verbose VERBOSE  Verbose level of the console log. 0 for silence, 1 for info level, 2 for debug level. Defaults to 1.
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
