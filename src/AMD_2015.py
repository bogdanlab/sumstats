import gzip, time
import numpy as np

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/ALZHEIMERS_2013/ALZHEIMERS_2013.txt'
out_fnm = 'ALZHEIMERS_2013.txt.gz'

# specify sample size here
ncase = 17008
ncontrol = 37154
ntotal = ncase + ncontrol

# create output file
out = gzip.open('./'+out_fnm, 'w')

# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tBETA\tSE\tN_CASE\tN_CONTROL\t\n')

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
    chrom_idx = 0
    pos_idx = 1
    snp_id_idx = 2
    effect_allele_idx = 3
    non_effect_allele_idx = 4
    beta_idx = 5
    se_idx = 6

    # parse out the fields
    chrom = cols[chrom_idx]
    pos = cols[pos_idx]
    snp_id = cols[snp_id_idx]
    effect_allele = cols[effect_allele_idx]
    non_effect_allele = cols[non_effect_allele_idx]
    beta = cols[beta_idx]
    se = cols[se_idx]

    # check for sanity of alleles
    if len(effect_allele) != 1 or len(non_effect_allele) != 1:
        print 'Removing SNP {} with alleles {}, {}'.format(snp_id,
            effect_allele, non_effect_allele)
        continue

    # check for sanity of beta
    if beta == 'NA' or se == 'NA':
        print 'Removing SNP {} with beta and se {}, {}'.format(snp_id,beta,se)
        continue
    
    # get z score
    zscore = np.float(beta) / np.float(se)
    
    # check for sanity of z score
    if np.isnan(zscore) or np.isinf(zscore):
        print 'Removing SNP {} with Z-score {}'.format(snp_id, zscore)
        continue

    # construct the output line
    # SNP CHR BP A1 A2 Z N BETA SE N_CASE N_CONTROL
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        beta,
        se,
        ntotal,
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
