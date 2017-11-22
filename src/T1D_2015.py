import gzip, time
import numpy as np
from utils import *

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'T1D_2015'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

ntotal = 24341

# create output file
out = gzip.open('./'+out_fnm, 'w')

# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\tOR\n')

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
    snp_id_idx = 0
    chrom_idx = 1
    pos_idx = 2
    pval_idx = 3
    or_idx = 4
    allele_idx = 7
    
    # get other information
    snp_id = cols[snp_id_idx]
    chrom = cols[chrom_idx]
    pos = cols[pos_idx]
    alleles = cols[allele_idx].split('>')
    effect_allele = alleles[1].upper()
    non_effect_allele = alleles[0].upper()
    odds_ratio = cols[or_idx]
    pval = cols[pval_idx]

    # check for sanity of alleles
    if len(effect_allele) != 1 or len(non_effect_allele) != 1:
        print 'Removing SNP {} with alleles {}, {}'.format(snp_id,
            effect_allele, non_effect_allele)
        continue

    # check for sanity of beta
    if not isfloat(odds_ratio):
        print 'Removing SNP {} with OR and SE {}, {}'.format(snp_id,
            odds_ratio)
        continue
    
    # get z score
    if not isfloat(odds_ratio) or not isfloat(pval):
        continue
    sign = '-'
    if np.log(np.float(odds_ratio)) > 0: sign = '+'
    zscore = ptoz(np.float(pval), sign)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N P OR
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntotal,
        odds_ratio,
        pval,
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
