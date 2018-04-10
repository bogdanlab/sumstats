import gzip, time
import numpy as np
from utils import *

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'MDD_EUR_2012'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

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

ncase = 9240
ncontrol = 9519
ntotal = ncase + ncontrol

# create output file
out = gzip.open('./'+out_fnm, 'w')

# snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tOR\tSE\tP\tINFO\tN_CASE\tN_CONTROL\n')

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
    effect_allele_idx = 3
    non_effect_allele_idx = 4
    or_idx = 5
    se_idx = 6
    pval_idx = 7
    info_idx = 8

    # parse out the fields
    snp_id = cols[snp_id_idx]

    # check if snp is in legend
    if snp_id not in legend:
        print 'SNP {} not in 1000G legend'.format(snp_id)
        continue
    
    # get other information
    chrom = legend[snp_id][0]
    pos = legend[snp_id][1]
    effect_allele = cols[effect_allele_idx].upper()
    non_effect_allele = cols[non_effect_allele_idx].upper()
    odds_ratio = cols[or_idx]
    se = cols[se_idx]
    pval = cols[pval_idx]
    info = cols[info_idx]

    # check for sanity of alleles
    if len(effect_allele) != 1 or len(non_effect_allele) != 1:
        print 'Removing SNP {} with alleles {}, {}'.format(snp_id,
            effect_allele, non_effect_allele)
        continue

    # check for sanity of beta
    if odds_ratio == 'NA' or se == 'NA':
        print 'Removing SNP {} with beta and se {}, {}'.format(snp_id,beta,se)
        continue
    
    # get z score
    if not isfloat(odds_ratio):
        continue
    zscore = np.log(np.float(odds_ratio)) / np.float(se)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N BETA SE P FREQ
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntotal,
        odds_ratio,
        se,
        pval,
        info,
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
