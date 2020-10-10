# imports extern
import argparse  # parsing arguments from terminal

# imports intern


parser = argparse.ArgumentParser()
parser.add_argument("--bed1", help="first input file as name or path", required=True)
parser.add_argument("--bed2", help="second input file as name or path", required=True)
parser.add_argument("--outfile", help="output file name or path", nargs='?')
parser.add_argument("--display", help="may be True for graphics in browser else always console output", nargs='?')
parser.add_argument("--alpha", help="custom level of significance, default is 0.05", nargs='?')
parser.add_argument("--freedom", help="for specific degree of freedom chi square test")

# for accessing parsed arguments
args = parser.parse_args()
