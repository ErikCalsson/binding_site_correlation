# imports extern
import argparse  # parsing arguments from terminal

# imports intern

parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')  # inactive
parser.add_argument("display_type", help=" may be dash for graphics or console for plane text", nargs='?')  # inactive
# TODO chi² arguments for p value and freedom degrees

# for accessing parsed arguments
args = parser.parse_args()
