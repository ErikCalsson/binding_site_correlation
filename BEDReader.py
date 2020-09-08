# packages:  file -> settings -> Project <name> -> Project Interpreter-> +
# #!!!! https://github.com/numpy/numpy/issues/16494

# imports
import argparse  # parsing arguments from terminal
import pandas as pd  # for data structuring
import numpy as np
import pybedtools

# adding arguments to parser
parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')

# for accessing parsed arguments
args = parser.parse_args()

# opening BED file or file from path
# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'] for access: cromEnd = end, chromStart = start !!!!!!
bed_one = pybedtools.BedTool(args.bed1)
bed_two = pybedtools.BedTool(args.bed2)
# intersect both files
interbothBED = bed_one.intersect(bed_two, s=True, r=True)  # s for overlap on same stand, r for A and B each overlap 90%

# print result to console
# print(interbothBED)

# put results in output file
# TODO here file names
outfile = interbothBED.saveas('intersection-output.bed', trackline='track name="intersection of both files')
print("saved as: ", outfile.fn)


# approach: number ob overlapped found relative to file
# by raw count number or percent
#   use jacquard or make own
# for own: intersection(A,B) / union(A,B) [union = A + B - intersection]

# calculate sequence length
def len_seq(bed_file):
    seq_len = 0
    # calculate length as sum of all chromEnd - chromStart
    for line in bed_file:
        seq_len += int(line[2]) - int(line[1])
    return seq_len

    # bed_file into pandas data_frame
    # data = pd.read_csv(bed_file, sep='\t', comment='t', header=None)
    # calculate length as sum of all chromEnd - chromStart
    # df = pd.DataFrame(data)
    # df[6] = df[2] - df[1]
    # seq_len = ((df[2] - df[1]).sum(axis=0))
    # return seq_len
    # return (df[2] - df[1]).sum(axis=0)


# own overlap quotient:
# TODO make own function, even later maybe in other file
# TODO improve by using log() for each variable
def overlap_quotient(bed_a, bed_b, intersection_ab):
    union_ab = len_seq(bed_a) + len_seq(bed_b) - len_seq(intersection_ab)
    print("normal overlap", len(intersection_ab) / union_ab)
    return np.math.log(len(intersection_ab)) / np.math.log(union_ab)


# degree of overlap for each file:
def overlap_file(bed, intersection_ab):  # log(A) natural, log(A,x) log to basis x
    print("same with log:                ", np.math.log(len(intersection_ab)) / np.math.log(len(bed)))
    return len_seq(intersection_ab) / len_seq(bed)


print("overlap quotient for both:    ", overlap_quotient(bed_one, bed_two, interbothBED))
print("overlap quotient for bed_one :", overlap_file(bed_one, interbothBED))
print("overlap quotient for bed_two :", overlap_file(bed_two, interbothBED))

# example output(normal and log):
# both:  0.8717 and 0.2555
# first: 0.4232 and 0.9151
# second 0.3921 and 0.9083

# own intersect idea: foreach start and stop in A within start and/or stop from B line to new file + count
# use .each from pybedtools
