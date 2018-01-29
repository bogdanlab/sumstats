import gzip, time
import numpy as np

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

# specify path to summary stats file here
trait = 'SCZ_ASN_2017'
root_dir = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/0_Raw'
sumstats_fnm = root_dir+'/{}/{}.txt'.format(trait, trait)
out_fnm = './{}.txt.gz'.format(trait)

ncase = 7699
ncontrol = 18327
ntot = ncase + ncontrol


# load legend
legend = dict()
legend_fnm = '/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/2_Final_Allele_Fixed/other/1000G_SNP_CHR_BP.txt'
legend_f = open(legend_fnm, 'r')
for line in legend_f:
    cols = line.strip().split()
    snp = cols[0]
    chrom = cols[1]
    bp = cols[2]
    legend[chrom+':'+bp] = snp
legend_f.close()
print '{} SNPs in legend'.format(len(legend))

# create output file
out = gzip.open('./'+out_fnm, 'w')

# write the header
out.write('SNP\tCHR\tBP\tA1\tA2\tZ\tN\tBETA\tSE\tP\tNCASE\tNCONTROL\n')

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
    pos_idx = 2
    snp_id_idx = 1
    effect_allele_idx = 3
    non_effect_allele_idx = 4
    beta_idx = 6
    se_idx = 7
    pval_idx = 5

    # parse out the fields
    chrom = cols[chrom_idx]
    pos = cols[pos_idx]
    snp_key = chrom+':'+pos
    if snp_key not in legend:
        continue
    else:
        snp_id = legend[snp_key]
    effect_allele = cols[effect_allele_idx]
    non_effect_allele = cols[non_effect_allele_idx]
    beta = np.log(np.float(cols[beta_idx]))
    se = cols[se_idx]
    pval = cols[pval_idx]

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
    # SNP CHR BP A1 A2 Z N BETA SE\tP NCASE NCONTROL
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        snp_id,
        chrom,
        pos,
        effect_allele,
        non_effect_allele,
        zscore,
        ntot,
        beta,
        se,
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
