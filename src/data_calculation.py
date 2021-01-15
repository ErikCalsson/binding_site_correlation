# imports extern
import numpy as np
import scipy.stats as stats
import pybedtools as pt
from itertools import product
import csv

# imports intern
import src.file_reader as file
import src.argument_parser as pars

# intersect both files
inter_both = file.bed_one.intersect(file.bed_two, s=True, r=True, bed=True, sorted=True)
# s: overlap on same stand, r: both overlap 90% each
# inter_both_BED = merge_one.intersect(merge_two, s=True, r=True, bed=True)


# sorting intersection for merge operation later
inter_both_BED = pt.BedTool.sort(inter_both)

# merging files with themselves to remove redundant overlaps
merge_one = pt.BedTool.merge(file.bed_one, s=True)
merge_two = pt.BedTool.merge(file.bed_two, s=True)
merge_inter = pt.BedTool.merge(inter_both_BED, s=True)


# calculate sequence length
def len_seq(bed_file):
    seq_len = 0
    # calculate length as sum of all chromEnd - chromStart
    for line in bed_file:
        seq_len += int(line[2]) - int(line[1])
    return seq_len


# save file lengths for reuse, applying merge to combine overlapping features on same strand into a single one
len_one = len_seq(merge_one)  # pt.BedTool.merge(file.bed_one, s=True))
len_two = len_seq(merge_two)  # pt.BedTool.merge(file.bed_two, s=True))
len_inter = len_seq(merge_inter)  # pt.BedTool.merge(inter_both_BED, s=True))


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
alpha = 0.05
a = np.math.log(len_one - len_inter)
b = np.math.log(len_one - (len_one - len_inter))
c = np.math.log(len_two - len_inter)
d = np.math.log(len_two - (len_two - len_inter))

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
    # reference https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chisquare.html
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
# reference https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
chi_results, chi_text = calc_chi(freedom, alpha)


# convert sequence length to appropriate bp unit
def conv_seq_len(seq_len):
    if seq_len > 1000:
        return str(round(seq_len / 1000, 1)) + "kbp"
    elif seq_len > 1000000:
        return str(round(seq_len / 1000000, 1)) + "mbp"
    else:
        return str(seq_len) + "bp"


# length of file and overlapping features in bp, kbp or mbp
bp_file_one = conv_seq_len(len_one)
bp_over_one = conv_seq_len(a)
bp_file_two = conv_seq_len(len_two)
bp_over_two = conv_seq_len(c)

# BedTools .getfasta for extracting sequence
# https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html


kMin = 1
if pars.args.kmin is not None:
    kMin = int(pars.args.kmin)

seq_one = merge_one.getfasta(fi=file.ref_fasta, s=True)  # bedOUT=True
seq_two = merge_two.getfasta(fi=file.ref_fasta, s=True)
seq_inter = merge_inter.getfasta(fi=file.ref_fasta, s=True)

kMax = kMin + 1  # 2
if pars.args.kmax is not None:
    kMax = int(pars.args.kmax)
else:
    if len_one < len_two and len_one:  # < len_inter:
        kMax = len_one
    else:  # elif len_two < len_one and len_two < len_inter:
        kMax = len_two
# else:
#     kMax = len_inter

# check whether kMin < kMax
if kMax < kMin:
    print("max kMer length < min kMer length")
    exit()
elif kMax == kMin:
    print("max kMer length = min kMer length")
    exit()


# generating dict's from kMin to kMax length for all possible k-mer's
def dic_mer_gen(k_min, k_max):
    # TODO ACGT not the only letters in fasta files!
    seq = ['A', 'C', 'G', 'T']
    tmp_dict = dict()
    for i in range(k_min + 1, k_max + 1):
        tmp_list = [''.join(k_mer) for k_mer in product(seq, repeat=i)]  # k-mer's length i
        for j in tmp_list:
            tmp_dict.update({j: 0})
    return tmp_dict


# template_dict = dic_mer_gen(kMin, kMin + 3)
# dict_one = template_dict  # dic_mer_gen(kMin, kMin + 3)
# dict_two = template_dict  #dict_one

dict_one = dic_mer_gen(kMin, kMax)
dict_two = dic_mer_gen(kMin, kMax)
# dict_one = dic_mer_gen(kMin, kMin + 3)
# dict_two = dic_mer_gen(kMin, kMin + 3)


# https://bioinformatics.stackexchange.com/questions/561/how-to-use-python-to-count-k-mers
# https://stackoverflow.com/questions/49188432/counting-maximum-k-mer-repetition-frequency


# iterating over sequences and counting k-mer appearances
def count_k_mer(dict_k_mer, seq, k_min, k_max):
    for mer_len in range(k_min, k_max + 1):
        for i in range(0, len(seq) - (mer_len - 1)):
            k = i+mer_len
            if seq[i:k] in dict_k_mer:
                dict_k_mer[seq[i:k]] += 1
    return dict_k_mer


# TODO remove following: test of results from small test
test_seq = '>Chr1CCCTAAACCCTAAACCCTAAACCC' \
           '>Chr2CCCTAAACCCTAAACCCTAAACCC' \
           '>Chr3CCCTAAACCCTAAACCCTAAACCC' \
           '>Chr4CCCTAAACCCTAAACCCTAAACCC' \
           '>Chr5CCCTAAACCCTAAACCCTAAACCC'

test_seq2 = '>Chr1ATGTCGTTCCTTTTTCATAAACCC' \
            '>Chr2ATGTCGTTCCTTTTTCAGAAACCC' \
            '>Chr3ATGTCGTTCCTTTTTCATCATCTT' \
            '>Chr4ATGTCGTTCCTTTTTCATCATCTT' \
            '>Chr5ATGTCGTTCCTTTTTCATCATCTT'


dict_one = count_k_mer(dict_one, seq_one, kMin, kMax)
#dict_one = count_k_mer(dict_one, test_seq, kMin, kMin+3)
#f = open("dictOne.txt", "w")
#f.write(str(dict_one))
#f.close()
w = csv.writer(open("dictOne.csv", "w"))
for key, val in dict_one.items():
    w.writerow([key, val])


dict_two = count_k_mer(dict_one, seq_two, kMin, kMax)
#dict_two = count_k_mer(dict_two, test_seq2, kMin, kMin+3)
#f = open("dictTwo.txt", "w")
#f.write(str(dict_two))
#f.close()
w = csv.writer(open("dictTwo.csv", "w"))
for key, val in dict_two.items():
    w.writerow([key, val])

# TODO: results saved as: key_in_one, key_in_both, count
# k_mer_finding = pd.DataFrame({
#    "in one": [],
#    "in both": [],
#    "count": []
# })

in_one = 0
in_both = 0
total_count = 0


# scan both dictionary's and get counts of k_mer's  appearing in both or in one with > 0
for key in dict_one:
    if (dict_one[key] > 0) and (dict_two[key] > 0):
        in_both += dict_one[key] + dict_two[key]
    elif (dict_one[key] > 0) or (dict_two[key] > 0):
        in_one += dict_one[key] + dict_two[key]
    total_count += dict_one[key] + dict_two[key]


#for k in dict_one:
#    if dict_one[k] != 0 and dict_two[k] == 0:
#        in_one += dict_one[k]
#        total_count += dict_one[k]
#    elif dict_one[k] != 0 and dict_two[k] != 0:
#        in_both += dict_one[k] + dict_two[k]
#        total_count += dict_one[k] + dict_two[k]
#    elif dict_one[k] == 0 and dict_two[k] != 0:
#        in_one += dict_two[k]
#        total_count += dict_two[k]

#for k in dict_two:
#    if dict_two[k] != 0:
#        if dict_one[k] != 0:
#            in_both += dict_one[k] + dict_two[k]
#            total_count += dict_one[k] + dict_two[k]
#        else:
#            in_one += dict_two[k]
#            total_count += dict_two[k]
#    else:
#        if dict_one[k] != 0:
#            in_one += dict_one[k]
#            total_count += dict_one[k]


# total count smaller than in other(output is:toal=274, one=225, both=49 instead of  toal=540, one=491, both=49)
#for key_out in dict_one:
#    if dict_one[key_out] != 0:
#        if dict_two[key_out] != 0:
#            in_both += dict_one[key_out] + dict_two[key_out]
#            total_count += dict_one[key_out] + dict_two[key_out]
#        else:
#            in_one += dict_one[key_out]
#            total_count += dict_one[key_out]
#    if dict_two[key_out] != 0:
#        if dict_one == 0:
#            in_one += dict_two[key_out]
#            total_count += dict_two[key_out]


# TODO for chi-square:
# occurrence against non occurrence in first and second

# TODO compare result and calculate conclusion
print("total, one, both")
print(total_count)
print(in_one)
print(in_both)
if in_both != 0:
    kMer_val = in_both / total_count
    print(kMer_val, "% same")
    print(in_one / total_count, "% different")

