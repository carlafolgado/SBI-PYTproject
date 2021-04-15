from Bio import pairwise2
import copy

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
    for chain in object_copy[0]:
        N = 65
        while chain.get_id() in [a.get_id() for a in complex.get_chains()]:
            try:
                chain.id = chr(N)
            except ValueError:
                pass
            N += 1
    return object_copy

def extract_DNA_sequence(PDB_object):
    """
    Input PDB_object, checks if chain is DNA, return None or the sequence
    """
    sequence = ""
    my_dict = {}

    for chain in PDB_object.get_chains():
        for res in chain:
            if len(res.get_resname().strip()) != 2:
                # is NOT DNA
                break
            else:
                sequence += res.get_resname().strip()[-1]

        my_dict[chain.id] = sequence
        # restart sequence
        sequence = ""

    return my_dict

def total_DNA_extraction(total_DNA_path, compressed):
    """
    Reads total DNA path, returns PDB object with total DNA
    """
    return ParsePDB(total_DNA_path, compressed)


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

def Superimpose_on_DNA(total_DNA, complex_dict):
    """
    Input total_DNA object, complex_dict with chains etc.
    complex_dict[uniprot_id][PDB_id] = [object_list]
    """
    print(complex_dict)
    # start super_complex with total DNA
    super_complex = Structure.Structure(id = "SuperComplex")
    super_complex.add(total_DNA[0])
    print([a.id for a in super_complex.get_chains()])

    for uniprot_id, value in complex_dict.items():
        for PDB_id, object in value.items():

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
                                print("Length Error")
                                residues_DNA = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos+1:final_pos-1]
                                atoms_DNA = unfold_entities(residues_DNA, "A")
                                residues_compare = unfold_entities(object[0][chain], "R")[1:-1]
                                atoms_compare = unfold_entities(residues_compare, "A")

                            sup = Superimposer()
                            sup.set_atoms(atoms_DNA, atoms_compare)
                            sup.apply(object[0])

                            WritePDB(object, "provant.pdb")

                            for renamed_chain in rename(object, super_complex).get_chains():
                                if not CheckClashes(super_complex, renamed_chain):
                                    print(renamed_chain.get_parent().get_parent())

                                    super_complex[0].add(renamed_chain)


    WritePDB(super_complex, "super_complex.pdb")
    return
