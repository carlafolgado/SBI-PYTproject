import argparse, os, sys, re

def argparser():
    parser = argparse.ArgumentParser(description="This is the program for our SBI-PYT poject")

    parser.add_argument('-i', '--input-directory',
                        dest = "indir",
                        action = "store",
                        required = True,
                        help = "Directory where input files are located")

    parser.add_argument('-o', '--output-directory',
                        dest = "outdir",
                        action = "store",
                        required = True,
                        help = "Directory for output files")

    parser.add_argument('-s', '--stechiometry',
                        dest = "stec",
                        action = "store",
                        required = False,
                        default = None,
                        help = "File with stechiometry of the complex")

    parser.add_argument('-f', '--force',
                        dest = "force",
                        action = "store_true",
                        required = False,
                        default = False,
                        help = "Overwrite existing files in specified output path")

    parser.add_argument('-v', '--verbose',
                        dest = "verbose",
                        action = "store_true",
                        required = False,
                        default = False,
                        help = "Display verbose progression log of the program")

    parser.add_argument('-r','--rmsd',
                        dest="rmsd",
                        action="store",
                        required=False,
                        default=0.5,
                        help="Set rmsd threshold for chain superposition. Default value at 0.5")

    parser.add_argument('-d','--dna',
                        dest="total_DNA_path",
                        action="store",
                        required = False,
                        default = None,
                        help = "File with the template DNA as PDB object")

    return parser.parse_args()
