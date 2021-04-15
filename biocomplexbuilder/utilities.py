import argparse, os, sys, re, copy, gzip
from Bio.PDB import PDBParser, Superimposer, NeighborSearch, PDBIO, Structure
from Bio.PDB.Selection import unfold_entities
from .arguments import *
from Bio.PDB.Polypeptide import PPBuilder
from .DNAbased_utilities import *

options=argparser()

def ParsePDB(file, compressed):
    """
    Takes a PDB file and parses all the structures returning instances of PDB Structure() objects.
    """

    name = file.split('.')[0]

    # If the files are compressed, open it with gzip and pass the filehandle to Bio.PDB.PDBParser
    if compressed:
        file = gzip.open(file, "rt")
    return PDBParser(QUIET=True).get_structure(name, file)


def FindCoreChain(object_list):
    """
    Takes a list of PDB objects and returns the name of the chain appearing more times in the binary interactions.
    It is assumed the objects belong to the same PDB structure.
    """
    dict = {}
    # Finds the times each chain appears in the list of objects
    for pdb in object_list:
        for chain in pdb.get_chains():
            dict[chain.get_id()] = dict.setdefault(chain.get_id(), 0) + 1

    # find key with the most appearances
    return sorted(dict, key=dict.get, reverse=True)

def CheckClashes(structure, chain):
    """
    Checks for clashes at a radius = 2 between a PDB structure and all the items on the provided chain.
    Arguments:
        -structure: PDB structure
        -chain: PDB chain structure to check if it has clashes with the structure.
    """
    # declare NeighborSearch() object instance with all the atoms from the structure (model 0), that includes all chains of that structure.
    ns = NeighborSearch(unfold_entities(structure[0], 'A'))

    # iterate over atoms in chain, search for close residues
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
    Takes a PDB structure and an output name and generates a .pdb file with the provided name in the output directory
    """
    # declare PDBIO object
    io = PDBIO()
    io.set_structure(structure)
    io.save(options.outdir+name)
    return

def SuperimposeStructures(object_list, complex, RMSD_threshold):
    """
    Takes a list containing PDB structures and returns a structure object with the correct structures superimposed.
    """
    # Get core chain to start reconstruction
    possible_cores = FindCoreChain(object_list)
    core = possible_cores[0]
    print(core)

    if len(list(complex.get_chains())):
        print(list(complex.get_chains()))
        while core not in [chain.get_id() for chain in list(complex.get_chains())]:
            possible_cores.remove(core)
            try:
                core = possible_cores[0]
            except IndexError:
                sys.stderr.write("!!!! Chain names not compatible\n")
                exit(1)
            print(possible_cores[0])

    if options.verbose:
        sys.stderr.write("Chain defined as core to superimpose: %s\n" %(core))


    # Declare Superimpose object
    sup = Superimposer()
    ref_struct = None

    for structure in list(object_list):
        #print(structure)
        # select the first structure with the core chain to be the reference
        try:
            if core in structure[0] and not ref_struct:
                #print("ref_struct")
                ref_struct = copy.deepcopy(structure)
                complex.add(ref_struct[0])
                print(f"Removing {structure}")
                object_list.remove(structure) ##
                continue ##
                #print(complex)

        except:
            #print("Hola")
            pass
        # if the structure contains the core chain, superimpose that to the ref structure set above
        if core in structure[0] and (structure is not ref_struct):

            sup.set_atoms(list(ref_struct[0][core].get_atoms()),list(structure[0][core].get_atoms()))
            sup.apply(structure[0])

            RMSD = sup.rms

            if RMSD < RMSD_threshold:
                for chain in structure[0]:
                    #print(chain)
                    if chain.id != core:
                        # check for clashes before adding new chain to complex
                        if not CheckClashes(complex, chain):
                            #print("clashes")
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
                                sys.stderr.write("Added to the final complex:\n")
                                sys.stderr.write("\tChain %s\n" %(chain.id))

                print(f"Removing {structure}")
                object_list.remove(structure)

    return (object_list, complex, RMSD_threshold) ##



def complex_builder(object_list,RMSD_threshold, complex=complex ):

    object_list, complex, RMSD_threshold = SuperimposeStructures(object_list, complex, RMSD_threshold)

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
    Takes a possible model and returns the energy DOPE scoring obtained using MODELLER
    """
    # Generate a PDB file from the complex as a temp file
    WritePDB(complex, "temp_model.pdb")
    log.none()
    env = Environ()
    env.libs.topology.read('${LIB}/top_heav.lib')
    env.libs.parameters.read('${LIB}/par.lib')
    mdl = Model(env)
    mdl.read(options.outdir + "temp_model.pdb")
    # os.remove("temp_model.pdb")
    return Selection(mdl).assess_dope()


def EnergyScoring(complex):

    WritePDB(complex, "temp_model.pdb")
    log.none()
    env = Environ()
    env.libs.topology.read('${LIB}/top_heav.lib')
    env.libs.parameters.read('${LIB}/par.lib')
    mdl = complete_pdb(env, "test/complex.pdb", transfer_res_num=True)
    (molpdf, terms) = Selection(mdl).energy(edat=EnergyData(dynamic_sphere=True))
    return molpdf, terms


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
#        model = PDBParser(QUIET = True).get_structure(pdb_file.split(".")[0], pdb_file)[0]

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
                alignment = seq_comparison(seq, compare_seq, threshold=0.95)
                if alignment:
                    big_dictionary[compare_item] = i
                    print(alignment)
                else:
                    print("not aligned")

        big_dictionary[item] = i
        i += 1
        print(i)
    print(big_dictionary)

    for item, seq in big_dictionary.items():

        #item[0] is object, [0] for the model, item[1] for the chain
        try:
            item[0][0][item[1]].id = chr(seq+65)
        except ValueError:
            # in the cases two chains are equal in the same object, you can't rename them the same
            # instead we go for the lowercase version of the later (if it's homotrimer in the same PDB file the program breaks)
            item[0][0][item[1]].id = chr(seq+97)

    return

def information_extraction(object_list):
    """
    Get object list, construct dictonary with dict[Uniprot_id][PDB_id] = [list of objects]
    """
    my_dict = {}
    for object in object_list:
        print(object.id)
        namefile = object.id.split("/")[-1]
        uniprot_id = namefile.split(".")[0]
        PDB_id = namefile.split(".")[2]
        object.id = uniprot_id + "." + PDB_id
        my_dict.setdefault(uniprot_id, {})
        my_dict[uniprot_id].setdefault(PDB_id, []).append(object)

    return my_dict

def stoichiometry_extraction(stoichiometry_path):
    """
    Obtain the path to stechiometry file. Open it and parse. Return dictionary with dict[uniprot_id] = #appearances
    """
    my_dict = {}
    try:
        fh = open(stoichiometry_path, 'r')
        for line in fh:
            parts = line.strip().split(":")
            my_dict.setdefault(parts[0], parts[1])

    except OSError:
        sys.stderr.write("Could not open stoichiometry file %s" %(stoichiometry_path))
        sys.exit(1)
    finally:
        if fh:
            fh.close

    return my_dic
