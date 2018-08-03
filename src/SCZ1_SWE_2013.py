import gzip, time
import numpy as np
from utils import *

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'SCZ1_SWE_2013'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt.gz'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

# CHR POS RSID EA NEA BETA SE P N

# load legend
legend = dict()
legend_fnm = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/2_Final_Allele_Fixed/other/1000G_SNP_CHR_BP.txt'
legend_f = open(legend_fnm, 'r')
for line in legend_f:
    cols = line.strip().split()
    snp = cols[0]
    chrom = cols[1]
    bp = cols[2]
    legend[snp] = (chrom, bp)
legend_f.close()
print '{} SNPs in legend'.format(len(legend))

# create output file
out = gzip.open('./'+out_fnm, 'w')

# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\tOR\tSE\tINFO\n')

# sample size information
ntot = 32143

# iterate through the file
flr = False
sumstats_f = gzip.open(sumstats_fnm, 'r')
for line in sumstats_f:

    # skip first line
    if flr == False:
        flr = True
        continue

    # split up the line into columns
    cols = line.strip().split()

    # specify indices of the fields
    snp_id_idx = 0
    eff_allele_idx = 3
    noneff_allele_idx = 4
    or_idx = 6
    se_idx = 7
    pval_idx = 8
    info_idx = 5

    # parse out the fields
    snp_id = cols[snp_id_idx]
    
    # check if snp is in legend
    if snp_id not in legend:
        print 'SNP {} not in 1000G legend'.format(snp_id)
        continue
    
    # get other information
    chrom = legend[snp_id][0]
    pos = legend[snp_id][1]
    effect_allele = cols[eff_allele_idx].upper()
    non_effect_allele = cols[noneff_allele_idx].upper()
    odd_r = cols[or_idx]
    se = cols[se_idx]
    pval = cols[pval_idx]
    info = cols[info_idx]

    # check for sanity of alleles
    if len(effect_allele) != 1 or len(non_effect_allele) != 1:
        print 'Removing SNP {} with alleles {}, {}'.format(snp_id,
            effect_allele, non_effect_allele)
        continue

    # check for sanity of beta
    if not isfloat(odd_r):
        print 'Removing SNP {} with OR and SE {}, {}'.format(snp_id,
            odd_r)
        continue
    
    # get z score
    if not isfloat(se) or not isfloat(pval):
        continue
    zscore = np.log(float(odd_r))/float(se)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N P OR SE INFO
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntot,
        pval,
        odd_r,
        se,
        info
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
