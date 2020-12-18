# imports extern
import pybedtools
from pathlib import Path

# imports intern
import src.argument_parser as arg

# opening BED file or file from path
# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'] for access: cromEnd = end, chromStart = start !!!!!!
# TODO check for BED-6 file format
bed_one = pybedtools.BedTool(arg.args.bed1)
bed_two = pybedtools.BedTool(arg.args.bed2)
ref_fasta = arg.args.fasta
#ref_fasta = pybedtools.BedTool(arg.args.fasta)


# break if both files are the same
if bed_one == bed_two:
    print("same file can't be tested with themselves")
    exit()

# get file names
name_first = str(Path(arg.args.bed1).stem)
name_second = str(Path(arg.args.bed2).stem)
