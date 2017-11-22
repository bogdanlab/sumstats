import pandas as pd
import sys, gzip
import numpy as np

# get input and output file names
input_file = sys.argv[1]
output_file = sys.argv[2]

# load gwas data
sumstats = pd.read_table(input_file, delim_whitespace=True)

# filter out snps with strand ambiguous alleles
ambiguous = ["AT", "CG", "TA", "GC"]
alleles = sumstats['A1']+sumstats['A2']
sa_idx = np.where(alleles.isin(ambiguous))[0]
sumstats = sumstats.drop(sa_idx)
sumstats = sumstats.reset_index(drop=True)
print 'Droped {} SNPs with strand ambiguous alleles'.format(sa_idx.shape[0])

# if info is in data, only use info greater than 0.9
if 'INFO' in sumstats.columns:
    info_thres = 0.9
    avg_info = np.mean(sumstats['INFO'])
    print 'Average INFO is {}'.format(avg_info)
    info_idx = np.where(sumstats['INFO'] <= info_thres)[0]
    sumstats = sumstats.drop(info_idx)
    sumstats = sumstats.reset_index(drop=True)
    print 'Droped {} SNPs with INFO less than {}'.format(sa_idx.shape[0],
        info_thres)

# filter snp with too litle or too much sample size
mean_n = np.mean(sumstats['N'])
sd_n = np.std(sumstats['N'])
print 'Sample size before: {} (SD {})'.format(mean_n, sd_n)
n_idx = np.where((sumstats['N']>mean_n+2.0*sd_n) | \
                 (sumstats['N']>mean_n-2.0*sd_n))[0]
sumstats = sumstats.drop(n_idx)
sumstats = sumstats.reset_index(drop=True)
print 'Droped {} SNPs with very large or small sample size'.format(
    n_idx.shape[0])
mean_n = np.mean(sumstats['N'])
sd_n = np.std(sumstats['N'])
print 'Sample size after: {} (SD {})'.format(mean_n, sd_n)

# write output
sumstats.to_csv(output_file, sep='\t', index=False, compression='gzip')
