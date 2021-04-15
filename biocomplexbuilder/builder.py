from .arguments import *
from .utilities import *
from .DNAbased_utilities import *

### PARSING AND CHECKING ARGUMENTS ###
options = argparser()

# checking if indir is a directory
if not os.path.isdir(options.indir):
    sys.stderr.write(options.indir + " is not a directory")

# adjusting indir name to add "/" if input was only the folder name
if not options.indir[-1] == "/":
    # indir not ends with /
    options.indir += "/"

# checking if force is false and output directory already exists. If so, stop program
if not options.outdir[-1] == "/":
    options.outdir += "/"

if os.path.isdir(options.outdir):
    if not options.force:
        sys.stderr.write("Output directory already exists\n")
        sys.exit(1)
else:
    # create output directory if it doesn't exist
    try:
        os.mkdir(options.outdir)
        print(options.outdir)
    except OSError:
        sys.stderr.write("Creation of the directory %s failed\n" % path)

# finding if files in indir are compressed or not
compressed = False
for file in os.listdir(options.indir):
    if re.match(r".*(.gz)$", file):
        compressed = True # if one file is compressed, assume all are compressed

RMSD_threshold = options.rmsd

if options.verbose:
    num_pdb=0
    for file in os.listdir(options.indir):
        if os.path.isfile(options.indir + file):
            num_pdb +=1

    sys.stderr.write("\n%s PDB found in %s\n\n" %(num_pdb, options.indir))

### ###

### PARSING PDB OBJECTS ###
PDB_objects = []

for file in os.listdir(options.indir):
    filename = str(options.indir) + str(file) # "/" is added before in this script

    # skip folders located on the input directory
    if os.path.isdir(filename):
        continue

    PDB_objects.append(ParsePDB(filename, compressed))
    if options.verbose:
        sys.stderr.write("\t%s transformed to PDB object\n" %(file))
### ###

if options.verbose:
    sys.stderr.write("\n\t\t####################\n")
    sys.stderr.write("\n\t\tMACROCOMPLEX BUILDER\n\n")
    sys.stderr.write("\t\t####################\n\n")
    sys.stderr.write("## Creating structure object from PDB files...\n")



# Create empty complex
complex = Structure.Structure(id = "complex")
if options.verbose:
    sys.stderr.write("\n## BUILDING MACROCOMPLEX\n")

#data_extraction(PDB_objects)
#complex_builder(PDB_objects,RMSD_threshold, complex)

#print(DOPEscoring(complex))

def random():
    ### 2O61 DNA example
    # total_DNA_path = "examples/2O61/infb_dna.pdb"
    total_DNA = ParsePDB(options.total_DNA_path, compressed)
    super_complex = Structure.Structure(id = "SuperComplex")
    #make a copy to rename chains before adding to super_complex
    total_DNA_copy = copy.deepcopy(total_DNA)

    for chain in total_DNA_copy[0]:

        N = 65
        while chain.get_id() in [a.get_id() for a in super_complex.get_chains()]:
            try:
                chain.id = chr(N)
            except ValueError:
                pass
            N += 1

    super_complex.add(total_DNA_copy[0])
    del total_DNA_copy

    print(extract_DNA_sequence(total_DNA))

    for object in PDB_objects:

        # extract DNA chains in the object
        for item, value in extract_DNA_sequence(object).items():

            # iterate each DNA sequence on our objects with each DNA chain on the template
            for DNA_chain, DNA_seq in extract_DNA_sequence(total_DNA).items():
                curr_chain = DNA_chain
                print(curr_chain)
                result = (seq_comparison(value, DNA_seq, 0.999))

                if result:
                    #print(result)
                    initial_pos = result[0][3]
                    final_pos = result[0][4]
                    # Dice.extract(total_DNA, 'X', initial_pos, final_pos, "test/slicing.pdb")
                    residues = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos:final_pos]
                    #print("!!!!")
                    #print(object)
                    total_DNA_atoms = (unfold_entities(residues, "A"))
                    print([a.get_id() for a in object[0].get_chains()])

                    if len(total_DNA_atoms) == len(unfold_entities(object[0][item], "A")):

                    # print(len(total_DNA_atoms))
                    # print(len(unfold_entities(object[0][item], "A")))
                        sup = Superimposer()
                        sup.set_atoms(total_DNA_atoms, unfold_entities(object[0][item], "A"))
                        sup.apply(object[0])
                        WritePDB(object, "provant.pdb")

                    else:
                        # if not same length of atoms, it's because one of the flanking nucleotides is incomplete (proof is that if we slice the 2 extreme nucleotides out, then the superimposer works)
                        print("Length Error")
                        residues = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos+1:final_pos-1]
                        total_DNA_atoms = unfold_entities(residues, "A")
                        partial_res = unfold_entities(object[0][item], "R")[1:-1]

                        sup = Superimposer()
                        sup.set_atoms(total_DNA_atoms, unfold_entities(partial_res, "A"))
                        sup.apply(object[0])
                        WritePDB(object, "provant2.pdb")

                    # adding all chains of object to super_complex only when alignment is there
                    for chain in object[0]:
                        if not CheckClashes(super_complex, chain):
                            #print(chain.get_id())
                            #print([a.get_id() for a in super_complex.get_chains()])
                            N = 65
                            while chain.get_id() in [a.get_id() for a in super_complex.get_chains()]:
                            #    print([a.get_id() for a in super_complex.get_chains()])
                                try:
                                    chain.id = chr(N)
                                except ValueError:
                                    pass
                                N += 1
                            super_complex[0].add(chain)
                        else:
                            print("clashes")
                            DNA_seq = DNA_seq[:initial_pos] + DNA_seq[final_pos:]
                            threshold = 0.99 * 0.1
                            result2 = (seq_comparison(value, DNA_seq, 0.5))
                            while not result and threshold > 0.6:
                                result2 = seq_comparison(value, DNA_seq, threshold)
                                threshold *= 0.1
                            print(result2)

                            if result2:
                                print(len(result2))
                                for j in range(0, len(result2)):
                                    initial_pos = result2[j][3]
                                    final_pos = result2[j][4]
                                    print(initial_pos)
                                    print(final_pos)

                                    first_seq = result2[0][0][initial_pos-1:final_pos]
                                    split_seq = first_seq.split("-")

                                    print(first_seq)

                                    max_len = max(map(lambda i: len(i), split_seq))
                                    for seq in split_seq:
                                        if len(seq) == max_len:
                                            max_seq = seq


                                    initial_pos = result[0][0].find(max_seq)
                                    final_pos = initial_pos + len(max_seq)

                                    # value is DNA sequence of our current PDB object


                                    residues = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos:final_pos]
                                    total_DNA_atoms = (unfold_entities(residues, "A"))
                                    print(curr_chain)
                                    print(residues)

                                    our_DNA_residues = unfold_entities(object[0][item], "R")[value.find(max_seq):value.find(max_seq)+len(max_seq)]
                                    our_DNA_atoms = unfold_entities(our_DNA_residues, "A")

                                    print(our_DNA_residues)
                                    print(len(total_DNA_atoms))
                                    print(len(our_DNA_atoms))

                                    if len(total_DNA_atoms) == len(our_DNA_atoms) and len(total_DNA_atoms) !=0 :
                                        print("###############################\n\n\n")
                                        print(len(total_DNA_atoms))
                                        print(len(unfold_entities(object[0][item], "A")))
                                        sup = Superimposer()
                                        sup.set_atoms(total_DNA_atoms, our_DNA_atoms)
                                        sup.apply(object[0])
                                        # WritePDB(object, "provant.pdb")

                                    elif len(total_DNA_atoms) !=0 :
                                        # if not same length of atoms, it's because one of the flanking nucleotides is incomplete (proof is that if we slice the 2 extreme nucleotides out, then the superimposer works)
                                        print("Length Error")
                                        residues = unfold_entities(total_DNA[0][curr_chain], "R")[initial_pos+1:final_pos-1]
                                        total_DNA_atoms = unfold_entities(residues, "A")
                                        print(len(total_DNA_atoms))

                                        partial_res = our_DNA_residues[1:-1]
                                        partial_atoms = unfold_entities(partial_res, "A")

                                        sup = Superimposer()
                                        sup.set_atoms(total_DNA_atoms, partial_atoms)
                                        sup.apply(object[0])

                                        # WritePDB(object, "provant2.pdb")


                                    if not CheckClashes(super_complex, chain):
                                            #print(chain.get_id())
                                            #print([a.get_id() for a in super_complex.get_chains()])
                                        copy_chain=copy.deepcopy(chain)
                                        N = 65
                                        while copy_chain.get_id() in [a.get_id() for a in super_complex.get_chains()]:
                                            print([a.get_id() for a in super_complex.get_chains()])
                                            try:
                                                copy_chain.id = chr(N)
                                            except ValueError:
                                                pass
                                            N += 1

                                        super_complex[0].add(copy_chain)

                        WritePDB(super_complex, "super_mega_complex2.pdb")





    WritePDB(super_complex, "super_mega_complex2.pdb")
    return
#random()
