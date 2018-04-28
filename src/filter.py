import pandas as pd
import sys, gzip, time
import numpy as np

# get input and output file names
input_file = sys.argv[1]
output_file = sys.argv[2]

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print( "Huwenbo Shi")
print( "Command started at", cur_time)

# load gwas data
sumstats = pd.read_table(input_file, delim_whitespace=True,
    dtype={'CHR': np.unicode})
sumstats['CHR'] = sumstats['CHR'].astype(np.unicode)

##############################################################################
# drop non-autosome snps
nsnp_orig = sumstats.shape[0]
sumstats = sumstats[sumstats.CHR.apply(lambda x: x.isnumeric())]
sumstats = sumstats.reset_index(drop=True)
ndropped = -(sumstats.shape[0]-nsnp_orig)
print( 'Dropped {} non-autosommal SNPs'.format(ndropped))

#############################################################################
# drop snps without rs id
nsnp_orig = sumstats.shape[0]
sumstats = sumstats[sumstats.SNP.apply(lambda x: x[0:2]=='rs')]
sumstats = sumstats.reset_index(drop=True)
ndropped = -(sumstats.shape[0]-nsnp_orig)
print( 'Dropped {} SNPs without rs ID'.format(ndropped))

#############################################################################
# remove duplicates
nsnp_orig = sumstats.shape[0]
sumstats['CHR'] = sumstats['CHR'].astype(int)
sumstats['BP'] = sumstats['BP'].astype(int)
sumstats = sumstats.drop_duplicates('SNP', keep=False)
sumstats = sumstats.drop_duplicates('BP', keep=False)
sumstats = sumstats.sort_values(by=['CHR', 'BP'])  
sumstats = sumstats.reset_index(drop=True)
ndropped = -(sumstats.shape[0]-nsnp_orig)
print( 'Dropped {} SNPs with duplicated ID or BP'.format(ndropped))
sumstats['CHR'] = sumstats['CHR'].astype(int)
sumstats['BP'] = sumstats['BP'].astype(int)

############################################################################
# filter out snps with strand ambiguous alleles
ambiguous = ["AT", "CG", "TA", "GC"]
alleles = sumstats['A1']+sumstats['A2']
sa_idx = np.where(alleles.isin(ambiguous))[0]
sumstats = sumstats.drop(sa_idx)
sumstats = sumstats.reset_index(drop=True)
print( 'Droped {} SNPs with strand ambiguous alleles'.format(sa_idx.shape[0]))

###########################################################################
# if info is in data, only use info greater than 0.9
if 'INFO' in sumstats.columns:
    info_thres = 0.9
    avg_info = np.mean(sumstats['INFO'])
    print( 'Average INFO is {}'.format(avg_info))
    info_idx = np.where(sumstats['INFO'] <= info_thres)[0]
    sumstats = sumstats.drop(info_idx)
    sumstats = sumstats.reset_index(drop=True)
    print( 'Droped {} SNPs with INFO less than {}'.format(sa_idx.shape[0],
        info_thres))

##########################################################################
# filter snp with too litle or too much sample size
print( 'Max sample size: {}'.format(np.max(sumstats['N'])))
print( 'Min sample size: {}'.format(np.min(sumstats['N'])))
mean_n = np.mean(sumstats['N'])
sd_n = np.std(sumstats['N'])
std_thres = 5.0
print( 'Mean sample size before: {} (SD {})'.format(mean_n, sd_n))
n_idx = np.where((sumstats['N'] > mean_n + std_thres*sd_n) | \
                 (sumstats['N'] < mean_n - std_thres*sd_n))[0]
sumstats = sumstats.drop(n_idx)
sumstats = sumstats.reset_index(drop=True)
print( 'Droped {} SNPs with very large or small sample size'.format(
    n_idx.shape[0]))
mean_n = np.mean(sumstats['N'])
sd_n = np.std(sumstats['N'])
print( 'Sample size after: {} (SD {})'.format(mean_n, sd_n))

#########################################################################
# write output
sumstats.to_csv(output_file, sep='\t', index=False, compression='gzip')

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print( "Command finished at", cur_time)
