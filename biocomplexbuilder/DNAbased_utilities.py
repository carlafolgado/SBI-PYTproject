from Bio import pairwise2
import copy
from utilities import *
import random

def DNA_filter(object_list):
    """
    Gets object list, returns list with only those PDB objects containing DNA
    """
    filtered_list = []
    for struct in object_list:
        for chain in struct[0]:
            for res in chain:
                if len(res.get_resname().strip()) == 2:
                    # is DNA
                    filtered_list.append(chain)
                break

    return

def rename(object, complex):
    """
    Input PDB object and complex. Returns object renamed so that it doesn't coincidence with chain names in complex.
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
    Input a PDB_entity and returns True or False if the residue sequence is DNA or not
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
    Input PDB_object, checks if chain is DNA, return None or the sequence
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

def total_DNA_extraction(total_DNA_path, compressed):
    """
    Reads total DNA path, returns PDB object with total DNA
    """
    return ParsePDB(total_DNA_path, compressed)

def construct_by_PDB_id(info_dict):
    """
    Takes dictionary with info_dict[uniprot_id][PDB_id] = [Objects]
    Modifies dictionary in place to return formed PDB complex.
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
            WritePDB(complex, complex.id + ".pdb")

    return info_dict

#Compare chain sequences
def seq_comparison(seq1, seq2, threshold=0.999):
    """
    Uses pariwise alignment to align a pair of sequences and 	returns True if percentage of identity is greater than the threshold, and False if it is not.
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
            my_dict.setdefault(parts[0], int(parts[1]))

    except OSError:
        sys.stderr.write("Could not open stoichiometry file %s" %(stoichiometry_path))
        sys.exit(1)
    finally:
        if fh:
            fh.close

    return my_dict

def Superimpose_on_DNA(total_DNA, complex_dict, RMSD_threshold, out_stoich_dict={}):
    """
    Input total_DNA object, complex_dict with chains etc.
    complex_dict[uniprot_id][PDB_id] = [object_list]
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

                            WritePDB(object, "provant.pdb")
                            RMSD = sup.rms

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
