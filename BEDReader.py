# packages:  file -> settings -> Project <name> -> Project Interpreter-> +
# #!!!! https://github.com/numpy/numpy/issues/16494

# imports
import argparse  # parsing arguments from terminal
import pandas as pd  # for data structuring
import numpy as np
import pybedtools
import dash  # visualisation of data via dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
from dash.dependencies import Input, Output
import scipy.stats as stats

# --------------------------------------------

# start dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# --------------------------------------------
# adding arguments to parser
parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')  # inactive
parser.add_argument("display_type", help=" may be dash for graphics or console for plane text", nargs='?')  # inactive
# last one later for deciding if results be shown in dash or terminal

# for accessing parsed arguments
args = parser.parse_args()

# --------------------------------------------
# opening BED file or file from path
# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'] for access: cromEnd = end, chromStart = start !!!!!!
# TODO check for BED-6 file format
bed_one = pybedtools.BedTool(args.bed1)
bed_two = pybedtools.BedTool(args.bed2)

if args.bed1 == '*.bed':
    print("snldanb")

# intersect both files
interbothBED = bed_one.intersect(bed_two, s=True, r=True)  # s for overlap on same stand, r for A and B each overlap 90%


# put results in output file
# TODO here file names
# outfile = interbothBED.saveas('intersection-output.bed', trackline='track name="intersection of both files')
# print("saved as: ", outfile.fn)


# approach: number ob overlapped found relative to file
# by raw count number or percent
#   use jacquard or make own
# for own: intersection(A,B) / union(A,B) [union = A + B - intersection]


# --------------------------------------------
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


# calculate overlap quotient with log:
def overlap_quotient_log(bed_a, bed_b, intersection_ab):
    union_ab = len_seq(bed_a) + len_seq(bed_b) - len_seq(intersection_ab)
    # print("normal overlap", len(intersection_ab) / union_ab)
    return np.math.log(len_seq(intersection_ab)) / np.math.log(union_ab)


# calculate overlap quotient lazy:
def overlap_quotient_regular(bed_a, bed_b, intersection_ab):
    union_ab = len_seq(bed_a) + len_seq(bed_b) - len_seq(intersection_ab)
    # print("normal overlap", len(intersection_ab) / union_ab)
    return len_seq(intersection_ab) / union_ab


# degree of overlap for each file with log:
def overlap_file_log(bed, intersection_ab):
    return np.math.log(len(intersection_ab)) / np.math.log(len(bed))


# degree of overlap for each file lazy:
def overlap_file_regular(bed, intersection_ab):
    return len_seq(intersection_ab) / len_seq(bed)


# print("overlap quotient for both:    ", overlap_quotient(bed_one, bed_two, interbothBED))
# print("overlap quotient for bed_one :", overlap_file(bed_one, interbothBED))
# print("overlap quotient for bed_two :", overlap_file(bed_two, interbothBED))

# example output(normal and log):
# both:  0.2555 and 0.8717
# first: 0.4232 and 0.9151
# second 0.3921 and 0.9083

# TODO remove absolute values later ?
both_files_lazy = overlap_quotient_regular(bed_one, bed_two, interbothBED)
first_file_lazy = overlap_file_regular(bed_one, interbothBED)
second_file_lazy = overlap_file_regular(bed_two, interbothBED)

# assign output values, relatives with log2
both_files_log = overlap_quotient_log(bed_one, bed_two, interbothBED)
first_file_log = overlap_file_log(bed_one, interbothBED)
second_file_log = overlap_file_log(bed_two, interbothBED)
# TODO: add values of total sequence length, percent of coverage. BUT not as part of graph figure


# --------------------------------------------

# chi_square test x: a=covered and b=not covered from 1 vs y: c=covered und d=not covered from 2
len_one = len_seq(bed_one)
len_two = len_seq(bed_two)
a = len_one * first_file_log  # = overlap_file_log(bed_one, interbothBED)
b = len_one - a  # a + b = len_one
c = len_two * second_file_log  # overlap_file_log(bed_two, interbothBED)
d = len_two - c


# X² = (n(a*d - c*b)²)/((a+c)(b+d)(a+b)(c+d)) has to be smaller then 3,841 to be accepted
# chi_result = ((len_one + len_two)*(a * d - c * b)**2) / ((a + c)*(b + d)*(a + b)*(c + d))

# validated_own = 'Values'
# validated_own += str(chi_result)  # TODO value for display only, remove later

# if chi_result >= 3.841:  # chi²_(0.95,1) --> 3.841
#    validated_own += ' are '
# else:
#    validated_own += ' may not '
# validated_own += 'statistical significant different'

# scipy version:
# https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html

validated = 'Values'
sci_out = stats.chi2_contingency([[a, b], [c, d]])
# output x², p-value and degree_of_freedom

if sci_out[1] < 0.05:  # only for chi²_(0.95)
    validated += ' are '
else:
    validated += ' may not '
validated += 'statistical significant different'
validated += '::'
validated += str(sci_out)  # TODO value for display only, remove later

# --------------------------------------------
# visualisation of output data

# dataframe for output
df = pd.DataFrame({
    "Overlap": ["Quotient", "First File", "Second File", "Quotient", "First File", "Second File"],
    "Coverage": [both_files_lazy, first_file_lazy, second_file_lazy,
                 both_files_log, first_file_log, second_file_log],
    "Size": ["Absolute", "Absolute", "Absolute",
             "Log 2", "Log 2", "Log 2"]
})

# figure
fig = px.bar(df, x="Overlap", y="Coverage", color="Size", barmode="group")

# app layout
app.layout = html.Div(children=[

    # title for the webpage
    html.H1(children="Overlap between both BED-files", style={'text-align': 'center'}),

    html.Div(id='output_container', children=[]),

    html.Br(),

    dcc.Graph(id='overlap_files', figure=fig),

    html.Br(),

    # html.H3(children=validated_own, style={'text-align': 'left'}),
    html.H3(children=validated, style={'text-align': 'left'})

])

# run app
app.run_server(debug=True)  # debug=true   means update browser by code change
# see in browser http://127.0.0.1:8050/ ONLY development server
