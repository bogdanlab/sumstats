import gzip
import pandas as pd

legend_dir = '/u/project/pasaniuc/pasaniucdata/DATA/1000GP_Phase3'



for i in range(1,23):
    print i    
    
    chrom_snp = []
    chrom_chr = []
    chrom_bp = []
    chrom_effect_allele = []
    chrom_non_effect_allele = []

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
        effect_allele = cols[3]
        non_effect_allele = cols[2]
        if len(effect_allele) != 1 or len(non_effect_allele) != 1:
            continue

        chrom_snp.append(snp)
        chrom_chr.append(i)
        chrom_bp.append(cols[1])
        chrom_effect_allele.append(effect_allele)
        chrom_non_effect_allele.append(non_effect_allele)

    legend_f.close()

    all_chrom = pd.DataFrame({'SNP': chrom_snp,
                              'CHR': chrom_chr,
                              'BP': chrom_bp,
                              'A1': chrom_effect_allele,
                              'A2': chrom_non_effect_allele})
    all_chrom = all_chrom[['SNP', 'CHR', 'BP', 'A1', 'A2']]
    all_chrom = all_chrom.drop_duplicates('SNP', keep=False)
    all_chrom = all_chrom.drop_duplicates('BP', keep=False)
    all_chrom.to_csv('1000G_SNP_CHR_BP_A1_A2.txt', '\t', mode='a', header=False, index=None)
