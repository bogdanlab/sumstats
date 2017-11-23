import sys, gzip, time
from utils import *

# get input and output file names
legend_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

# define equivalent alleles
equiv = dict()
equiv["AC"] = set(["TG", "AC", "TC", "AG"])
equiv["AG"] = set(["TC", "AG", "TG", "AC"])
equiv["CA"] = set(["GT", "CA", "GA", "CT"])
equiv["CT"] = set(["GA", "CT", "GT", "CA"])
equiv["TC"] = set(["AG", "TC", "AC", "TG"])
equiv["TG"] = set(["AC", "TG", "AG", "TC"])
equiv["GA"] = set(["CT", "GA", "CA", "GT"])
equiv["GT"] = set(["CA", "GT", "CT", "GA"])

# define reversed alleles
reverse = dict()
reverse["AC"] = set(["GT", "CA", "CT", "GA"])
reverse["AG"] = set(["CT", "GA", "GT", "CA"])
reverse["CA"] = set(["TG", "AC", "AG", "TC"])
reverse["CT"] = set(["AG", "TC", "TG", "AC"])
reverse["TC"] = set(["GA", "CT", "CA", "GT"])
reverse["TG"] = set(["CA", "GT", "GA", "CT"])
reverse["GA"] = set(["TC", "AG", "AC", "TG"])
reverse["GT"] = set(["AC", "TG", "TC", "AG"])

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Huwenbo Shi"
print "Command started at", cur_time

##############################################################################
# load legend
legend = dict()
legend_file = gzip.open(legend_file, 'r')
for line in legend_file:
    cols = line.strip().split()
    snp = cols[0]
    chrom = cols[1]
    bp = cols[2]
    a1 = cols[3]
    a2 = cols[4]
    legend[snp] = {'CHR': chrom, 'BP': bp, 'A1': a1, 'A2': a2}
    
    # for debugging only
    if chrom == '2':
        break
legend_file.close()
print '{} SNPs loaded from legend'.format(len(legend))

##############################################################################
# iterate through sumstats file
output_file = gzip.open(output_file, 'w')
input_file = gzip.open(input_file, 'r')
flr = False
for line in input_file:

    # figure out header
    if flr == False:
        output_file.write(line.strip()+'\n')
        header = line.strip().split()
        header_idx = dict()
        for field in header:
            header_idx[field] = header.index(field)
        flr = True
        continue

    # check the columns
    cols = line.strip().split()
    SNP = cols[0]
    CHR = cols[1]
    BP = cols[2]
    A1 = cols[3]
    A2 = cols[4]
    Z = float(cols[5])

    #########################################################################
    # make sure the snp info matched legend
    if SNP not in legend:
        print '{} not in 1000G legend'.format(SNP)
        continue
    if CHR != legend[SNP]['CHR'] or BP != legend[SNP]['BP']:
        print '{} does not match 1000G'.format(SNP)
        continue

    ########################################################################
    # match allele encoding
    A1A2 = A1+A2
    leg_A1A2 = legend[SNP]['A1']+legend[SNP]['A2']
    
    #######################################################################
    # in case they are equivalent
    if A1A2 in equiv[leg_A1A2]:
        
        # keep the allele
        cols[3] = legend[SNP]['A1']
        cols[4] = legend[SNP]['A2']
        
        # write output
        output_file.write('\t'.join(cols)+'\n')

    # in case they are reversed
    elif A1A2 in reverse[leg_A1A2]:

        # flip the allele
        cols[3] = legend[SNP]['A2']
        cols[4] = legend[SNP]['A1']

        # flip z score
        cols[5] = str(-1.0*Z)

        # flip beta if in header
        if 'BETA' in header:
            idx = header_idx['BETA']
            beta = cols[idx]
            if isfloat(beta):
                cols[idx] = str(-1.0*float(beta))

        # flip or if in header
        if 'OR' in header:
            idx = header_idx['OR']
            OR = cols[idx]
            if isfloat(OR):
                cols[idx] = str(1.0/float(OR))

        # flip freq if in header
        if 'FREQ' in header:
            idx = header_idx['FREQ']
            FREQ = cols[idx]
            if isfloat(FREQ):
                cols[idx] = str(1.0-float(FREQ))

        # write output
        output_file.write('\t'.join(cols)+'\n')

    # doesn't match
    else:
        print "{} dosn't have matched allele".format(SNP)
        continue

# close files
input_file.close()
output_file.close()

# print out time info
cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "Command finished at", cur_time
