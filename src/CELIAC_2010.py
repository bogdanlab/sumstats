import gzip, time
import numpy as np
from utils import *

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'CELIAC_2010'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

ncase = 4533
ncontrol = 10750
ntotal = ncase + ncontrol

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
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\tOR\tN_CASE\tN_CONTROL\n')

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
    allele_idx = 7
    pval_idx = 3
    or_idx = 4

    # parse out the fields
    snp_id = cols[snp_id_idx]
    
    # check if snp is in legend
    if snp_id not in legend:
        print 'SNP {} not in 1000G legend'.format(snp_id)
        continue
    
    # get other information
    chrom = legend[snp_id][0]
    pos = legend[snp_id][1]
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
    # SNP CHR BP A1 A2 Z N P OR N_CASE N_CONTROL
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntotal,
        pval,
        odds_ratio,
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
