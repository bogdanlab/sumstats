import gzip, time
import numpy as np
from utils import *

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'RA_TE_2014'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = '{}.txt.gz'.format(trait)

ncase = 19234
ncontrol = 61565
ntotal = ncase + ncontrol

# create output file
out = gzip.open('./'+out_fnm, 'w')

# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\tN_CASE\tN_CONTROL\t\n')

# iterate through the file
flr = False
sumstats_f = open(sumstats_fnm, 'r')
for line in sumstats_f:

    # skip first line
    if flr == False:
        flr = True
        continue

    # split up the line into columns
    cols = line.strip().split()

    # specify indices of the fields
    chrom_idx = 1
    pos_idx = 2
    snp_id_idx = 0
    effect_allele_idx = 3
    non_effect_allele_idx = 4
    or_idx = 5
    pval_idx = 8

    # parse out the fields
    chrom = cols[chrom_idx]
    pos = cols[pos_idx]
    snp_id = cols[snp_id_idx]
    effect_allele = cols[effect_allele_idx]
    non_effect_allele = cols[non_effect_allele_idx]
    odds_ratio = cols[or_idx]
    pval = cols[pval_idx]

    # check for sanity of alleles
    if len(effect_allele) != 1 or len(non_effect_allele) != 1:
        print 'Removing SNP {} with alleles {}, {}'.format(snp_id,
            effect_allele, non_effect_allele)
        continue

    # check for sanity of pval
    if pval == 'NA':
        print 'Removing SNP {} with pval {}'.format(snp_id, pval)
        continue
    
    # get z score
    sign = '-'
    log_or = np.log(np.float(odds_ratio))
    if log_or > 0: sign = '+'
    zscore = ptoz(float(pval), sign)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N P N_CASE N_CONTROL
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntotal,
        pval,
        ncase,
        ncontrol
    )

    # write the output
    out.write(outline)

# close input file
sumstats_f.close()

# close the output
out.close()

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Command finished at", cur_time
