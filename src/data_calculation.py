# imports extern
import numpy as np
import scipy.stats as stats
import pybedtools as pt


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

# ref seq: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.3/
# RefSeq; Genomic FASTA (.fna)

# BETTER SOURCE FOR FASTA !!!
#

kMin = 1
if pars.args.kmin is not None:
    kMin = pars.args.kmin

#seq_one = pt.BedTool.getfasta(file.ref_fasta, merge_one, s=True, bedOUT=True)
seq_one = merge_one.getfasta(fi=file.ref_fasta)  # pt.BedTool.getfasta(fi=file.ref_fasta, bed=merge_one, s=True, bedOUT=True)
#seq_two = pt.BedTool.getfasta(file.ref_fasta, merge_two, s=True, bedOUT=True)
seq_two = merge_two.getfasta(fi=file.ref_fasta)
#seq_inter = pt.BedTool.getfasta(file.ref_fasta, merge_inter, s=True, bedOUT=True)
seq_inter = merge_inter.getfasta(fi=file.ref_fasta)

# TODO k-mer analysis

