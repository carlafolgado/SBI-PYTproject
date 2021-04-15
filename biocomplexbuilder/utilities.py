import argparse, os, sys, re, copy, gzip
from Bio.PDB import PDBParser, Superimposer, NeighborSearch, PDBIO, Structure
from Bio.PDB.Selection import unfold_entities
from biocomplexbuilder.arguments import *
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

            RMSD = sup.rms

            # check for clashes before adding new chain to complex
            if RMSD < RMSD_threshold:
                for chain in structure[0]:

                    if chain.get_id() != core:

                        for renamed_chain in rename(object, super_complex).get_chains():
                            if not CheckClashes(complex, chain):
                                complex[0].add(chain_copy)

                                if options.verbose:
                                    sys.stderr.write("\tChain %s\n" %(chain.id))


                object_list.remove(structure)

    return (complex, object_list)

def complex_builder(object_list,RMSD_threshold, complex=complex ):
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

def data_extraction(object_list, threshold = 0.90):
    """
    Takes as input a list of pdb files and a fasta file
    Returns a dictionary of dictionaries. The primary key is the model, the secondary
    key is the chain_id in said model and the value is the fasta_id of said chain.
    This function finds sequences by pairwise sequence alignment. If a pdb chain
    aligns with 0.95 (by default) identity with a fasta sequence, it is given the fasta identifier.
    If an aminoacid sequence is fragmented because of discontinuity or some other reason, the
    function puts it together.
    """
    big_dictionary = {}

    for object in object_list:

        for chain in object.get_chains():
            pdb_seqs = PPBuilder().build_peptides(chain)

            # if it's a DNA sequence, list will be empty
            if len(pdb_seqs)>1:
                pp_seq = "".join(list([str(pp.get_sequence()) for pp in pdb_seqs]))

            elif len(pdb_seqs) == 1:
                pp_seq = pdb_seqs[0].get_sequence()

            else:
                #use our own function to extract the DNA sequence
                pp_seq = extract_DNA_sequence(object)[chain.get_id()]
            print(pp_seq)

            big_dictionary[(object,chain.get_id())] = pp_seq

    i = 0
    for item, seq in big_dictionary.items():

        if type(seq) is int:
            continue

        print(item[1])
        for compare_item, compare_seq in big_dictionary.items():
            print(compare_item[1])
            if type(compare_seq) is not int:
                alignment = seq_comparison(seq, compare_seq, threshold=0.99)
                if alignment:
                    big_dictionary[compare_item] = i
                    print(item, compare_item)
                    print(alignment)
                else:
                    print("not aligned")

        big_dictionary[item] = i
        i += 1
        print(i)


    for item, seq in big_dictionary.items():

        #item[0] is object, [0] for the model, item[1] for the chain
        try:
            item[0][0][item[1]].id = chr(seq+65)
        except ValueError:
            # in the cases two chains are equal in the same object, you can't rename them the same
            # instead we go for the lowercase version of the later (if it's homotrimer in the same PDB file the program breaks)
            item[0][0][item[1]].id = chr(seq+97)
    print(big_dictionary)
    return
