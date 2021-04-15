from Bio import pairwise2

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
