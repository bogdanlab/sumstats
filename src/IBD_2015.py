import gzip, time
import numpy as np

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'IBD_2015'

root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

# specify sample size here
ncase =12924
ncontrol = 21770
ntotal = ncase + ncontrol

# create output file
out = gzip.open('./'+out_fnm, 'w')

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
    snp_id_idx = 1
    chrom_idx = 0
    pos_idx = 2
    effect_allele_idx = 3
    non_effect_allele_idx = 4
    or_idx = 8
    se_idx = 9
    pval_idx = 10
    info_idx = 7

    # parse out the fields
    snp_id = cols[snp_id_idx]
    chrom = cols[chrom_idx]
    pos = cols[pos_idx]
    effect_allele = cols[effect_allele_idx]
    non_effect_allele = cols[non_effect_allele_idx]
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
        print 'Removing SNP {} with OR and se {}, {}'.format(snp_id,
            odds_ratio, se)
        continue
    
    # get z score
    zscore = np.log(np.float(odds_ratio)) / np.float(se)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N OR SE P INFO N_CASE N_CONTROL
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'\
    .format(
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
