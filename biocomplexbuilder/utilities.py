import argparse, os, sys, re, copy, gzip
from Bio.PDB import PDBParser, Superimposer, NeighborSearch, PDBIO, Structure
from Bio.PDB.Selection import unfold_entities
from arguments import *
from Bio.PDB.Polypeptide import PPBuilder
from modeller import *
from modeller.scripts import complete_pdb

options=argparser()

def ParsePDB(file, compressed=False):
    """
    Takes a PDB file and parses all the chains returning an instance of PDB Structure() object.
    Accepts compressed files with argument compressed=True.
    """
    if options.total_DNA_path:
        name = file.split('_')[0]
    else:
        name = file.split('.')[0] ###!!!

    # If the files are compressed, open it with gzip and pass the filehandle to Bio.PDB.PDBParser
    if compressed:
        file = gzip.open(file, "rt")
    return PDBParser(QUIET=True).get_structure(name, file)

def FindCoreChain(object_list):
    """
    Takes a list of PDB objects and returns a string with the name of the chain appearing more times in the binary interactions.
    It is assumed that the objects belong to the same PDB structure.
    """
    my_dict = {}
    # Finds the times each chain appears in the list of objects
    for pdb in object_list:
        for chain in pdb.get_chains():
            my_dict[chain.get_id()] = my_dict.setdefault(chain.get_id(), 0) + 1

    # find key with the most appearances
    return max(my_dict, key=my_dict.get)

def CheckClashes(structure, chain):
    """
    Checks for clashes at a radius = 2 between a PDB structure and all the atoms on the provided chain.
    Returns True or False depending if number of different residues clashing exceeds 25.
    Arguments:
        -structure: PDB structure.
        -chain: PDB chain structure to check if it has clashes with the main structure.
    """
    # declare NeighborSearch() object instance with all the atoms from the structure (model 0), that includes all chains of that structure.
    ns = NeighborSearch(unfold_entities(structure[0], 'A'))

    # iterate over atoms in input chain, search for close residues
    clashing_residues = set([])
    for atom in chain.get_atoms():
        close_res = ns.search(atom.get_coord(), radius=2, level="R")

        try:
            close_res.remove(atom.get_parent())
        except ValueError:
            pass

        for res in close_res:
            neighbor_res = (atom.get_parent(), res)
            clashing_residues.add(neighbor_res)
            if len(clashing_residues) > 25:
                return True

    return False

def WritePDB(structure, name):
    """
    Takes a PDB structure and an output name and generates a .pdb file with the provided name in the output directory.
    Function returns None.
    """
    # declare PDBIO object
    io = PDBIO()
    io.set_structure(structure)
    io.save(os.path.join(options.outdir, name))
    return

def SuperimposeStructures(object_list, complex, RMSD_threshold):
    """
    Superimposes chains from objects in object_list to chains in complex. Adds the non-clashing chains to the complex and removes the structure from the object_list.
    Returns the complex with the new added chains, and the updated object_list with

    Arguments:
        -object_list : list of PDB objects that have to be superimposed and added to the complex.
        -complex: main structure to which individual chains from the object_list have to be added after superimposition.
        -RMSD_threshold: threshold for the RMSD value of the superposition between a chain of an object an the same chain on the complex.
        Default value for the program is 0.5.
    """

    # Get core chain to start reconstruction
    core = FindCoreChain(object_list)
    if options.verbose:
        sys.stderr.write("Chain defined as core to superimpose: %s\n" %(core))
        sys.stderr.write("Added to the final complex:\n")

    # Declare Superimpose object
    sup = Superimposer()
    ref_struct = None

    for structure in list(object_list):

        # select the first structure with the core chain to be the reference
        try:
            if core in structure[0] and not ref_struct:
                ref_struct = copy.deepcopy(structure)
                complex.add(ref_struct[0])

        except:
            pass

        # if the structure contains the core chain, superimpose that to the chain with same name in ref structure set before
        if core in structure[0] and (structure is not ref_struct):

            sup.set_atoms(unfold_entities(ref_struct[0][core], 'A'), unfold_entities(structure[0][core], 'A'))
            sup.apply(structure[0])

            RMSD = float(sup.rms)
            print(RMSD)
            # check for clashes before adding new chain to complex
            if RMSD < RMSD_threshold:
                for chain in structure[0]:

                    if chain.get_id() != core:
                        if not CheckClashes(complex, chain):
                            chain_copy = copy.deepcopy(chain)

                            N = 65
                            while chain_copy.get_id() in [a.get_id() for a in complex.get_chains()]:
                                try:
                                    chain_copy.id = chr(N)
                                except ValueError:
                                    pass
                                N += 1

                            complex[0].add(chain_copy)

                            if options.verbose:
                                sys.stderr.write("\tChain %s\n" %(chain.id))


                object_list.remove(structure)

    return (complex, object_list)

def complex_builder(object_list, RMSD_threshold, complex=complex ):
    """
    Takes an object_list and a complex and calls SuperimposeStructures to add the superimposable chains from objects in object_list to the complex.
    The function calls itself recursively until no more objects remain in object_list (which will be removed upon addition to the complex by the SuperimposeStructures function).
    Generates a a file in the output directory with the complex once all possible chains are added, and returns None.
    Arguments:
        -object_list: list of PDB objects. It is assumed that the objects all belong to the same PDB structure.
        -RMSD_threshold: threshold for the RMSD value of the superposition between a chain of an object an the same chain on the complex.
        Default value for the program is 0.5.
        -complex: main PDB structure where non-clashing chains will be added.
    """

    complex, object_list = SuperimposeStructures(object_list, complex, RMSD_threshold)

    if len(object_list):
        return complex_builder(object_list,RMSD_threshold,complex)
    else:
        if options.verbose:

            sys.stderr.write("All chains added\n")

            sys.stderr.write("\n#### Complex finished ####\n")

            sys.stderr.write("\nNew 'complex.pdb' file created in folder called %s \n\n"%(options.outdir))
        WritePDB(complex, "complex.pdb")

def DOPEscoring(complex):
    """
    Takes a PDB object and calculates the DOPE energy score of the model using MODELLER package.
    To do so generates a temp file which is then removed.
    Returns an int with the DOPE energy score value.
    """
    # Generate a PDB file from the complex as a temp file
    WritePDB(complex, "temp_model.pdb")
    log.none()
    env = Environ()
    env.libs.topology.read('${LIB}/top_heav.lib')
    env.libs.parameters.read('${LIB}/par.lib')
    mdl = Model(env)
    mdl.read(os.path.join(options.outdir, "temp_model.pdb"))
    os.remove(os.path.join(options.outdir, "temp_model.pdb"))
    return int(Selection(mdl).assess_dope())
