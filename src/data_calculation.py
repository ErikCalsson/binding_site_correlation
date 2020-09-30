# imports extern
import numpy as np
import scipy.stats as stats

# imports intern
import src.file_reader as file


# intersect both files
inter_both_BED = file.bed_one.intersect(file.bed_two, s=True, r=True)
# s: overlap on same stand, r: both overlap 90% each


# put results in output file
# TODO here file names
# out_file_intersection = inter_both_BED.saveas('intersection-output.bed',
#                                               trackline='track name="intersection of both files')
# print("saved as: ", out_file_intersection.fn)


# calculate sequence length
def len_seq(bed_file):
    seq_len = 0
    # calculate length as sum of all chromEnd - chromStart
    for line in bed_file:
        seq_len += int(line[2]) - int(line[1])
    return seq_len


# save file lengths for reuse
len_one = len_seq(file.bed_one)
len_two = len_seq(file.bed_two)
len_inter = len_seq(inter_both_BED)


# calculate overlap quotient with log:
def overlap_quotient_log(len_first, len_second, len_inter_both):
    union_ab = len_first + len_second - len_inter_both
    # print("normal overlap", len(intersection_ab) / union_ab)
    return np.math.log(len_inter_both) / np.math.log(union_ab)


# calculate overlap quotient lazy:
def overlap_quotient_regular(len_first, len_second, len_inter_both):
    union_ab = len_first + len_second - len_inter_both
    # print("normal overlap", len(intersection_ab) / union_ab)
    return len_inter_both / union_ab


# degree of overlap for each file with log:
def overlap_file_log(len_bed, len_inter_both):
    return np.math.log(len_inter_both) / np.math.log(len_bed)


# degree of overlap for each file lazy:
def overlap_file_regular(len_bed, len_inter_both):
    return len_inter_both / len_bed


# assign output values,
# TODO remove absolute values later ?
both_files_lazy = overlap_quotient_regular(len_one, len_two, len_inter)
first_file_lazy = overlap_file_regular(len_one, len_inter)
second_file_lazy = overlap_file_regular(len_two, len_inter)

# assign output values, relatives with log2
both_files_log = overlap_quotient_log(len_one, len_two, len_inter)
first_file_log = overlap_file_log(len_one, len_inter)
second_file_log = overlap_file_log(len_two, len_inter)


# chi_square test x: a=covered and b=not covered from 1 vs y: c=covered und d=not covered from 2
a = len_one * first_file_log  # = overlap_file_log(bed_one, inter_both_BED)
b = len_one - a  # a + b = len_one
c = len_two * second_file_log  # overlap_file_log(bed_two, inter_both_BED)
d = len_two - c

# TODO adapt for usage of specific parameters for p-value
# reference: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html
sci_out = stats.chi2_contingency([[a, b], [c, d]])
