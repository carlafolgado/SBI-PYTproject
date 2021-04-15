**BioComplex Builder Manual**
===================================
*by* Arnau Llinàs, Carla Folgado i Oriol Canal
*MSc Bioinformatics for Health Sciences, Pompeu Fabra University, Barcelona, 2021.*


<br />

# Introduction

BioComplex Builder is a bioinformatic tool written in Python program whose objective is to model macrocomplexes of proteins from PDB files of interacting chains and based on this information, construct the most likely macro-complex in a PDB file. BioComplex Builder is not only compatible with pairs of protein chains, in addition it can deal with DNA-proteins interactions to reconstruct the final complex.

You can source the code at the [github repository](https://github.com/carlafolgado/SBIPYTproject).
<br /><br />

# Biological background


Proteins are supramolecular structures that are comprised of 20 natural amino acids joined by peptide bonds, which constructs the primary structure. They play a central role in biological organisms by transferring information across the cell (i.e cell signalling), catalysing thousands of complex chemical reactions and transformingchemical energy into biological work.[1](
https://pubs.rsc.org/en/content/getauthorversionpdf/C3CS60474H)

The function of a protein is directly dependent on its 3D structure which is determined by the sequence of amino acids. However, proteins rarely act alone. It has been revealed that over 80% of proteins do not operate alone, but in complexes.[2](https://www.hindawi.com/journals/ijpro/2014/147648/)

Proteins interact with  other proteins or other biological molecules (DNA, RNA, lipids among others) through different interactions as electrostatic forces, hydrogen bonding or hydrophobic effect. The final result of these interactions is the formation of macrocomplexes allowing to acquire new functional capabilities to proteins. So, proteins have to be studied in the context of their interacting partners to understand their function.

So, protein-protein interactions (PPi) are crucial for the formation of macromolecules that are the basis of most cellular processes as signal transduction, activation or inhibition of proteins. To reveal the function mechanisms in cells, it is important to identify PPIs that take place in the organism.

In order to identify and study the PPIs, it has been used high-throughput experiments and computational methods. Small-scale experimental approaches can detect high quality PPI but with high time cost. On the other hand, computational methods and high throughput methods can detect a large amount of novel PPIs but with low quality and consequently with high number of false positives.[3](https://academic.oup.com/bib/article/18/5/798/2562794)

<br /><br />

# BioComplex algorithm

The main objective of the program is to reconstruct protein complexes (which may include DNA) from  PDB files of interacting chains.
<br />

### Parsing PDB files and stochiometry file

From the input PDB files (which can be compressed or not) containing information about two interacting chains that can be whether protein-protein interaction or DNA-protein interactions, the algorithm parse these files and transform it to PDB objects. In addition, if the stochiometry file is given in the arguments, it parses the file to return a dictionary where the key is the uniprot_id and the values the number of appearances.

<br />

### Constructing model without DNA chain as template

  **First, the common chains in different PDB files are identified and are considered the same PDB object**. The following step is to look for the most common chain of the PDB files of interacting chains and take it as core which will be added to the complex. Then by superimposing  the complex structure with a common chain of a PDB file, the new chain is placed on the interacting site. The superimposition RMSD is checked and if the value is lower than a threshold and the new coordinates of the new chain have no clashes, the new chain is added to the macrocomplex. In this way the complex now have 3 chains. By using an iterative algorithm the different PDB interactions are superposed to the complex and the new chains that have not clashes and good RMSD value are added to the final complex. The algorithm of BioComplex Builder is recursive until it has exhausted all the PDB objects. Finally it is obtained a PDB file for the macrocomplex.

**The solution given is not heuristic (add some information HEREEE!!!**
<br />

### Constructing model using DNA template

In order to construct the complex with a DNA sequence, the PDB file of the whole DNA sequence (template DNA) have to be given using the argument --dna. The DNA template sequence will be used as the initial complex in order to reconstruct the complete complex.

Then it starts an iterative process for each PDB object given as input:

* First, for the PDB objects that we have parsed, we extract the DNA sequence (interacting DNA) and it is performed a local alignment between the template DNA and the interacting DNA in order to find the position where the interacting chains have to be plced. For the interacting DNA that have been alignmed with a %identity higher than **99,9**. 
* From the local alignment the initial and final position of the alignment are stored in order to identify the region of the template DNA where the alignment has occurred.  

* In order to perform the superimposition, the superimposer method needs 2 objects with the same number of atoms. However, sometimes the final and the initial residues of the interacting DNA are not complete and consequently it can't be applied the superimposer method.  To avoid this problem, the algorithm detects if the number of atoms between the interacting chain and the nucleotides aligned in the template chain are not the same. If it is the case, the initial and final residues of the interacting chain are removed allowing to have a correct superimposition between the DNA sequences. 

* Once solved this problem, we can superimpose the DNA sequence that contain the interaction pair to the template DNA sequence io order to obtain the correct coordinates of the chain that interacts with the DNA.

* Finally, the correct coordinates of the chain are added in the final complex if no clashes occurs between the chain added and the complex.

* Performing this process iterating for all the PDB objects that has been given as input, the final complex is built.


<br /><br />


# Installation
### Prerequisites:
- Python 3.0: `https://www.python.org/download/releases/3.0/`
- Biopython package: `https://biopython.org/wiki/Download`

### Installation from GitHub repositori

```shell
$ git clone https://github.com/carlafolgado/SBI-PYTproject.git
$ cd SBI-PYT_Project
$ sudo python3 setup.py install
```
<br /><br />

# Organization of the program scripts

The program have been organized in different scripts:

* **builder.py:** <p>Where workflow of the program can be found. </p>
* **arguments.py:** <p>In this script are defined what arguments requires the program, and which ones are optional. It also automatically generates help and usage messages and issues errors when users give the program invalid arguments. </p>
* **utilities.py:**  <p>Where the functions for the executation of the program are stored.</p>
* **DNAbased_utilities.py:** <p>Where functions for the executation of the program when a DNA chain is given for the construction of the complex are stored.</p>

<br /><br />

# Argument options

BioComplex Builder have different arguments which some of them are mandatory to run the program and other are optional which help to deal with other information as how to give the total chain of DNA where the complex have to be built. 
<br />

### Command line options:

#### Required arguments:

* **-i, --input**
<p>Directory where input PDB files of interacting chains are located</p>


* **-o, --output**
<p>Directory where the PDB file of the macrocomplex constructed by the program will be found.</p>

#### Non required arguments:

* **-s, --stoichiometry**
<p>Path where the stoichiometry file of the complex is located</p>

* **-f, --force**
<p> When used this argument the created files of the program overwrite existing files in specified output path</p>

* **-v, --verbose**
<p> Display verbose progression log of the program</p>

* **--rmsd**
<p> Set rmsd threshold for chain superposition. Default value at 0.5 </p>

* **--dna**
<p> Path to the file with the template DNA when the macrocomplex have to be build using the DNA template as core. </p>







