import gzip

legend_dir = '/u/project/pasaniuc/pasaniucdata/DATA/1000GP_Phase3'

for i in range(1,23):
    legend_fnm = '{}/1000GP_Phase{}_chr2.legend.gz'.format(legend_dir, i)
    print legend_fnm
