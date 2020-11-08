# imports extern
import numpy as np
import scipy.stats as stats


# imports intern
import src.file_reader as file
import src.argument_parser as pars


# intersect both files
inter_both_BED = file.bed_one.intersect(file.bed_two, s=True, r=True, wb=True)
# s: overlap on same stand, r: both overlap 90% each

print(inter_both_BED)


# calculate sequence length
def len_seq(bed_file):
    seq_len = 0
    # calculate length as sum of all chromEnd - chromStart
    for line in bed_file:
        seq_len += int(line[2]) - int(line[1])
    return seq_len


# save file lengths for reuse
len_one = len_seq(set(file.bed_one))
len_two = len_seq(set(file.bed_two))
len_inter = len_seq(set(inter_both_BED))


# calculate overlap quotient with log:
def overlap_quotient_log(len_first, len_second, len_inter_both):
    union_ab = len_first + len_second - len_inter_both
    return np.math.log(len_inter_both) / np.math.log(union_ab)


# calculate overlap quotient lazy:
def overlap_quotient_regular(len_first, len_second, len_inter_both):
    union_ab = len_first + len_second - len_inter_both
    return len_inter_both / union_ab


# degree of overlap for each file with log:
def overlap_file_log(len_bed, len_inter_both):
    return np.math.log(len_inter_both) / np.math.log(len_bed)


# degree of overlap for each file lazy:
def overlap_file_regular(len_bed, len_inter_both):
    return len_inter_both / len_bed


# assign % output values
both_files_lazy = overlap_quotient_regular(len_one, len_two, len_inter)
first_file_lazy = overlap_file_regular(len_one, len_inter)
second_file_lazy = overlap_file_regular(len_two, len_inter)

# assign % output values, relatives with log2
both_files_log = overlap_quotient_log(len_one, len_two, len_inter)
first_file_log = overlap_file_log(len_one, len_inter)
second_file_log = overlap_file_log(len_two, len_inter)


# chi_square test x: a=covered and b=not covered from 1 vs y: c=covered und d=not covered from 2
#a = len_one * first_file_log  # = overlap_file_log(bed_one, inter_both_BED)
#b = len_one - a  # a + b = len_one
#c = len_two * second_file_log  # overlap_file_log(bed_two, inter_both_BED)
#d = len_two - c
alpha = 0.05
print("iterln", len_inter)
print("len1:", len_one)
print("len2:", len_two)
a = (len_one - len_inter)
b = (len_one - (len_one - len_inter))
c = (len_two - len_inter)
d = (len_two - (len_two - len_inter))
#print("log1", first_file_log)
#print("1reg", first_file_lazy)
print("a", a)
print("b", b)
#print("log2", second_file_log)
print("c", c)
print("d", d)

# TODO use Fisher Test for better/advanced comparison of p-values?
res = stats.fisher_exact([[a, b], [c, d]])
print('fisher test result: ', res)


# break if both files are the same, when a = c = 0
if a == 0 and c == 0:
    print("same file can't be tested with themselves")
    exit()

# level of significance
if pars.args.alpha is not None:
    alpha = float(pars.args.alpha)


# reference: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html
sci_out = stats.chi2_contingency([[a, b], [c, d]])  # output: x², p-value and degree_of_freedom and expected values


# degree of freedom
if pars.args.freedom is not None:
    freedom = int(pars.args.freedom)
else:
    freedom = sci_out[2]


# calculate chi-square test
def calc_chi(free, alp):
    val_str = 'Values'
    chi_value = stats.chisquare([[a, b], [c, d]], axis=None, f_exp=sci_out[3], ddof=free)
    if chi_value[1] > alp:
        val_str += ' are not '
    else:
        val_str += ' may be '
    val_str += 'statistical significant different with:'
    return chi_value, val_str


# TODO use Fisher Test for better/advanced comparison of p-values?
res = stats.fisher_exact([[a, b], [c, d]])
print('fisher test result: ', res)

# perform chi-square test for program start
chi_results, chi_text = calc_chi(freedom, alpha)
