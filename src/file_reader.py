# imports extern
import pybedtools

# imports intern
import src.argument_parser as arg

# opening BED file or file from path
# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'] for access: cromEnd = end, chromStart = start !!!!!!
# TODO check for BED-6 file format
bed_one = pybedtools.BedTool(arg.args.bed1)
bed_two = pybedtools.BedTool(arg.args.bed2)
