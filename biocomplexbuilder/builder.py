from arguments import *
from utilities import *
from DNAbased_utilities import *

### PARSING AND CHECKING ARGUMENTS ###
options = argparser()

# checking if indir is a directory
if not os.path.isdir(options.indir):
    sys.stderr.write( "> ERROR: "+ options.indir + " is not a directory. Argument -i must be a folder containing PDB files \n")
    exit(1)


# checking if force is false and output directory already exists. If so, stop program
if os.path.isdir(options.outdir):
    if not options.force:
        sys.stderr.write("> ERROR: Output directory already exists, Choose another one or use -f option\n")
        exit(1)
else:
    # create output directory if it doesn't exist
    try:
        os.mkdir(options.outdir)
    except OSError:
        sys.stderr.write("Creation of the directory %s failed\n" % path)

# finding if files in indir are compressed or not
compressed = False
for file in os.listdir(options.indir):
    if re.match(r".*(.gz)$", file):
        compressed = True # if one file is compressed, assume all are compressed

# check if we have stoichiometry as input, and open it if it is a file
if options.stoic:
    if os.path.isfile(options.stoic):
        stoich_dict = stoichiometry_extraction(options.stoic)

if options.total_DNA_path:
    if os.path.isfile(options.total_DNA_path):
        #try:
        total_DNA = total_DNA_extraction(options.total_DNA_path, compressed)
        #except Exception:
        #    sys.stderr.write("Could not parse DNA file %s\n" %(options.total_DNA_path))
        #    sys.exit(1)

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

if options.total_DNA_path:
    if options.verbose:
        sys.stderr.write("\t# DNA set as template for building the complex\n")
    if not options.stoic:
        stoich_dict = {}

    info_dict = information_extraction(PDB_objects)
    stoich_dict = stoich_dict
    complex_dict = construct_by_PDB_id(info_dict)
    scoring_tuples = []
    i = 0
    RMSD_threshold = 15
    while i < int(options.models):

        current_complex = Superimpose_on_DNA(total_DNA, complex_dict, RMSD_threshold, stoich_dict)
        WritePDB(current_complex, "def_complex." + str(i) + ".pdb")
        if options.verbose:
            sys.stderr.write("\nComplex "+ str(i) + " built and saved as"+ "def_complex." + str(i) + ".pdb "+"in %s\n"%(options.outdir))
        score = DOPEscoring(current_complex)
        scoring_tuples.append((str(i), score))

        i += 1

    with open(os.path.join(options.outdir, "scores.dope"), 'w') as fh:
        fh.write("Complex\t\t\t\tDOPE Score\n")
        fh.write("---------------------------------------------------\n")
        scoring_tuples.sort(key=lambda x:x[1])
        for score in scoring_tuples:
            fh.write("\t" + str(score[0]) + "\t\t\t\t%.3f\n" %(score[1]))

    if options.verbose:
        sys.stderr.write("\n New 'scores.dope' created in folder called %s" %(options.outdir))

else:
    # No DNA to superimpose on
    complex_builder(PDB_objects, RMSD_threshold, complex)
