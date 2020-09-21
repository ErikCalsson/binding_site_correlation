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

# --------------------------------------------

# start dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# --------------------------------------------
# adding arguments to parser
parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')

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
# TODO validate output of Data: Overlap in A, B and between both

# Chi-Quadrat, Willcocs, Fisher


# TODO remove simple output test later
# Test for output validation HERE overlap in log >= 0.9
is_validate = 'Following files are statistical significant: '
not_validate = 'Following files are statistical not significant: '
if both_files_log >= 0.9:
    is_validate += 'quotient, '
else:
    not_validate += 'quotient, '
if first_file_log >= 0.9:
    is_validate += 'first file, '
else:
    not_validate += 'first file, '
if second_file_log >= 0.9:
    is_validate += 'second file, '
else:
    not_validate += 'second file, '

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
fig = px.bar(df, x="Overlap", y="Coverage", color="Size",  barmode="group")

# app layout
app.layout = html.Div(children=[

    # title for the webpage
    html.H1(children="Overlap between both BED-files", style={'text-align': 'center'}),

    html.Div(id='output_container', children=[]),

    html.Br(),

    dcc.Graph(id='overlap_files', figure=fig),

    html.Br(),

    html.H3(children=is_validate, style={'text-align': 'left'}),
    html.H3(children=not_validate, style={'text-align': 'left'})

])

# run app
app.run_server(debug=True)  # debug=true   means update browser by code change
# see in browser http://127.0.0.1:8050/ ONLY development server
