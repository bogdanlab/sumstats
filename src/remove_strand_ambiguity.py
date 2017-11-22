import pandas as pd
import sys, gzip

# define strand ambiguous alleles
ambiguous = set(["AT", "CG", "TA", "GC"])

# get input and output file names
input_file = sys.argv[1]
output_file = sys.argv[2]

# load gwas data
sumstats = pd.read_table(input_file, )
