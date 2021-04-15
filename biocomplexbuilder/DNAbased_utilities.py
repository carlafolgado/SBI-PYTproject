from Bio import pairwise2
import copy
from utilities import *
import random

def rename(object, complex):
    """
    Takes a PDB object and a complex.
    Returns a copy of the PDB object so that no chain in the object has the same id as a chain in the complex.
    Renaming is done checking non-used chain names in the complex, starting form ASCII character A.
    """

    # ASCII A-Za-z encoded on decimal 65 to 122
    object_copy = copy.deepcopy(object)
    for chain in object_copy.get_chains():
        N = 65
        while chain.get_id() in [a.get_id() for a in complex.get_chains()]:
            try:
                chain.id = chr(N)
            except ValueError:
                pass
            N += 1
    return object_copy

def is_DNA(PDB_entity):
    """
    Input a PDB_entity of any level.
    Returns True if the residue sequence is DNA or RNA. If not, return False.
    """
    # if it's atom get residue (parent)
    if PDB_entity.get_level() == "A":
        PDB_entity = PDB_entity.get_parent().get_parent() # get chain
    elif PDB_entity.get_level() == "R":
        PDB_entity = PDB_entity.get_parent() # get chain

    for res in PDB_entity.get_residues():
        if len(res.get_resname().strip()) <= 2:
            #is DNA
            return True
        else:
            return False

def extract_DNA_sequence(PDB_object):
    """
    Takes a PDB structure or model and checks if the contained chains are DNA or RNA.
    Return a dictionary with the names of the chains that are DNA/RNA as keys and the sequence as value.
    If there are no DNA/RNA chains, returns an empty dictoinary.
    """
    sequence = ""
    my_dict = {}

    for chain in PDB_object.get_chains():
        for res in chain:
            if not is_DNA(res):
                # is NOT DNA nor RNA
                break
            else:
                sequence += res.get_resname().strip()[-1]
        if sequence:
            my_dict[chain.id] = sequence
            # restart sequence
            sequence = ""

    return my_dict

def total_DNA_extraction(total_DNA_path, compressed=False):
    """
    Takes the path to the total DNA file and parses it.
    The file name must contain an underscore "_".
    Returns a PDB Structure() instance.
    Accepts compressed files with compressed=True.
    """
    return ParsePDB(total_DNA_path, compressed)

def construct_by_PDB_id(info_dict):
    """
    Takes a dictionary of dictionaries with the form:
    dict{ Uniprot id : {PDB id: [PDB objects],
                        PDB id: [PDB objects]},
                       {PDB id: [PDB objects]}
        }
    The dictionary contains a list of objects that come from the same PDB structure containing the same protein from Uniprot id.
    Overlaps all the objects under the same PDB id and Uniprot id to form the corresponding PDB complex.
    Modifies the same dictionary in place so that the values of the inner dictoinary are now the PDB complex formed.
    Final dictionary is of the form:
    dict{ Uniprot id: { PDB id: PDB object,
                        PDB id: PDB object,},
                      { PDB id: PDB object}
        }
    Returns the modified dictionary.
    """

    for uniprot_id, value in info_dict.items():
        # item is uniprot id, value is dict with PDB_id = [objects]


        for pdb_id, objects in value.items():

            # re-start complex
            complex = Structure.Structure(id = uniprot_id + "." + pdb_id)
            ref_struct = None

            for object in objects:
                if not ref_struct:
                    ref_struct = object[0]
                    complex.add(ref_struct)
                    continue
                for chain in object.get_chains():
                    if chain not in complex.get_chains():
                        complex[0].add(chain)


            info_dict[uniprot_id][pdb_id] = complex

    return info_dict

#Compare chain sequences
def seq_comparison(seq1, seq2, threshold=0.999):
    """
    Takes two strings with DNA sequences and a identity threshold (0-1).
    Uses pariwise alignment to align the sequence, giving a -0.5 penalty for gap opening and -0.1 for gap extension.
    The default identity threshold is set to 0.999 to only identify identic regions on both sequences.
    Returns True if identity is greater than the provided threshold, or False if not.
    Arguments:
        - seq1: string with a DNA sequence (can be previously obtained with extract_DNA_sequence() function).
        - seq2: string with a DNA sequence (can be previously obtained with extract_DNA_sequence() function).
        - threshold: int with identity threshold (0-1). Default is set to 0.999.
    pair of sequences and 	returns True if percentage of identity is greater than the threshold, and False if it is not.
    """
    alignment = pairwise2.align.localxs(seq1, seq2, -0.5, -0.1)

    if alignment:
        score = alignment[0][2]
        ident_perc = score / len(seq1)

        if ident_perc > threshold:
            return alignment
    else:
        return False

def information_extraction(object_list):
    """
    Takes a PDB object list with ids of the form path/Uniprot_id.DNA.PDB_id (coming from the Parsed filenames).
    Extracts the Uniprot id and PDB id information from the object id.
    Returns a dictionary of dictionaries where outer keys are Uniprot ids, inner keys are PDB ids and the values are lists of all the objects that have the given Uniprot id and the PDB id.
    The dictionary will have the form:
    dict{ Uniprot id : {PDB id: [PDB objects],
                        PDB id: [PDB objects]},
                       {PDB id: [PDB objects]}
        }
    """
    my_dict = {}
    for object in object_list:
        namefile = object.id.split("/")[-1]
        uniprot_id = namefile.split(".")[0]
        PDB_id = namefile.split(".")[2]
        object.id = uniprot_id + "." + PDB_id
        my_dict.setdefault(uniprot_id, {})
        my_dict[uniprot_id].setdefault(PDB_id, []).append(object)

    return my_dict

def stoichiometry_extraction(stoichiometry_path):
    """
    Takes the path to the stoichiometry file and parses it.
    Return a dictionary with the protein Uniprot id as keys and the number of appearances in the final complex as values.
    """
    my_dict = {}
    try:
        fh = open(stoichiometry_path, 'r')
        for line in fh:
            parts = line.strip().split(":")
            my_dict.setdefault(parts[0].strip(), int(parts[1].strip()))

    except OSError:
        sys.stderr.write("Could not open stoichiometry file %s" %(stoichiometry_path))
        sys.exit(1)
    finally:
        if fh:
            fh.close

    return my_dict

def Superimpose_on_DNA(total_DNA, complex_dict, RMSD_threshold, out_stoich_dict={}):
    """
    Superimposes different chains onto a DNA template based on the interaction of said chains with a given DNA region.
    The function follows the stoichiometry provided, not adding the same chain more times than dicated by the input.
    Returns a PDB Structure() object with the complex formed.
    Arguments:
        - total_DNA: PDB Structure() object with 2 chains of DNA to use as a template.
        - complex dict: source of chains to be superimposed on to the template DNA. dictionary of the form
            dict{Uniprot id: { PDB id: PDB object,
                               PDB id: PDB object,},
                             { PDB id: PDB object}
                }
        - RMSD_threshold: int with the RMSD threshold for the superimpositions. Those chain-DNA interaction where the DNA superimposes on the template DNA with a higher RMSD than the set threshold will not be added to the final complex.
        - out_stoich_dict: dictionary containing the stoichiometry relations. Keys are Uniprot id and values the number of appearances in the final complex.
        Default value is an empty dictionary, that by the behaviour of the function will allow 1 appearance of each kind of protein on the complex_dict.

    From the dictionary, the PDB objects are accessed randomly and the chains with DNA are extracted and superimposed to the best region on the template DNA. Then the proteic chain is added to the complex if no clashes appear.
    The process is repeated with each chain adding it if there are no clashes and the stoichiometry allows it, until all the stoichiometric conditions are fulfilled.
    If the clashes do not allow to complete all the stoichiometric requirements, the function calls itself and restarts the complex until a complex that fulfills the conditions is obtained.
    """

    # fresh copy of total DNA to re-start the supre_complex
    new_total_DNA_model = total_DNA[0].copy()

    stoich_dict = copy.copy(out_stoich_dict)

    # start super_complex with total DNA

    super_complex = Structure.Structure(id = "super_complex")
    super_complex.add(new_total_DNA_model)

    keys = list(complex_dict.keys())
    random.shuffle(keys)

    for key in keys:
        uniprot_id = key
        value = complex_dict[uniprot_id]

        stoich_dict.setdefault(uniprot_id, 1) # if the uniprot_id was not in stoich, add it with value 1 (to only add it once to the complex)

        for PDB_id, object in value.items():
            if not stoich_dict[uniprot_id] > 0:
                continue


            # extract DNA chains in the object
            for chain, seq in extract_DNA_sequence(object).items():
                # iterate each DNA sequence on our objects with each DNA chain on the template
                for DNA_chain, DNA_seq in extract_DNA_sequence(total_DNA).items():
                    curr_chain = DNA_chain

                    results = (seq_comparison(seq, DNA_seq, 0.999))

                    if results:
                        for result in results:

                            # get alignment positions of sequence
                            initial_pos = result[3]
                            final_pos = result[4]
                            residues_DNA = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos:final_pos]
                            atoms_DNA = unfold_entities(residues_DNA, "A")

                            atoms_compare = unfold_entities(object[0][chain], "A")

                            if not len(atoms_DNA) == len(atoms_compare):
                                # if not same length of atoms, it's because one of the flanking nucleotides is incomplete (proof is that if we slice the 2 extreme nucleotides out, then the superimposer works)
                                residues_DNA = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos+1:final_pos-1]
                                atoms_DNA = unfold_entities(residues_DNA, "A")
                                residues_compare = unfold_entities(object[0][chain], "R")[1:-1]
                                atoms_compare = unfold_entities(residues_compare, "A")

                            sup = Superimposer()
                            sup.set_atoms(atoms_DNA, atoms_compare)
                            sup.apply(object[0])

                            RMSD = float(sup.rms)

                            if RMSD < RMSD_threshold:

                                for renamed_chain in rename(object, super_complex).get_chains():
                                    if not stoich_dict[uniprot_id] > 0:
                                        continue

                                    if is_DNA(renamed_chain):
                                        continue

                                    if not CheckClashes(super_complex, renamed_chain):
                                        orig_struct = renamed_chain.get_parent().get_parent() # obtain structure
                                        super_complex[0].add(renamed_chain)
                                        stoich_dict[uniprot_id] = int(stoich_dict[uniprot_id]) - 1 # convert on int !!


    if sum(stoich_dict.values()) == 0:
        return super_complex
    else:
        super_complex = Superimpose_on_DNA(total_DNA, complex_dict, RMSD_threshold, out_stoich_dict)
        return super_complex
