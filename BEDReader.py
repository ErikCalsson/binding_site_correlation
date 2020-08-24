# packages:  file -> settings -> Project <name> -> Project Interpreter-> +
# imports
import argparse  # parsing arguments from terminal
import pandas as pd  # for data structuring
import numpy as np
import pybedtools
from pybedtools import BedTool

# adding arguments to parser
parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')

# for accessing parsed arguments
args = parser.parse_args()

# opening BED file or file from path
# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
bed_one = pybedtools.BedTool(args.bed1)
bed_two = pybedtools.BedTool(args.bed2)
# intersect both files
interbothBED = bed_one.intersect(bed_two)
# print result to console
print(interbothBED)
# put results in output file
outfile = interbothBED.saveas('intersection-output.bed', trackline='track name="intersection of both files') # TODO here file names
print("saved as: ", outfile.fn)

