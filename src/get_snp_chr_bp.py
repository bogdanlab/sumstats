import gzip
import pandas as pd

legend_dir = '/u/project/pasaniuc/pasaniucdata/DATA/1000GP_Phase3'

chrom_snp = []
chrom_chr = []
chrom_bp = []

for i in range(1,23):
    print i    
    legend_fnm = '{}/1000GP_Phase3_chr{}.legend.gz'.format(legend_dir, i)
    flr = False
    legend_f = gzip.open(legend_fnm, 'r')
    for line in legend_f:
        if flr == False:
            flr = True
            continue
        cols = line.strip().split()
        snp = cols[0].split(':')[0]
        if snp[0:2] != 'rs': continue

        chrom_snp.append(snp)
        chrom_chr.append(i)
        chrom_bp.append(cols[1])
    
    legend_f.close()

all_chrom = pd.DataFrame({'SNP': chrom_snp, 'CHR': chrom_chr, 'BP': chrom_bp})
all_chrom = all_chrom[['SNP', 'CHR', 'BP']]
all_chrom = all_chrom.drop_duplicates('SNP', keep=False)
all_chrom = all_chrom.drop_duplicates('BP', keep=False)
all_chrom.to_csv('1000G_SNP_CHR_BP.txt', '\t', header=False, index=None)
