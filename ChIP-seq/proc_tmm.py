import pandas as pd
import sys
import re
import subprocess as sp
from os import listdir


# read factor
def main_process(file):
    try:
        histone = re.search('(H3K\d+\w+\.).', file).group(1)
    except:
        histone = 'noknown'

    factors = {}
    with open('tmm_factor.txt', 'r') as f:
        for line in f:
            a = line.rstrip().split('\t')
            factors[a[0]] = float(a[1])

    if re.search('k4me3', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k4me3']
        else:
            factor = factors['factor.k4me3']
    elif re.search('k36', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k36']
        else:
            factor = factors['factor.k36']
    elif re.search('k4me1', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k4me1']
        else:
            factor = factors['factor.k4me1']
    elif re.search('k27ac', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k27ac']
        else:
            factor = factors['factor.k27ac']
    elif re.search('k9me3', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k9me3']
        else:
            factor = factors['factor.k9me3']
    elif re.search('k27me3', histone, re.I):
        if re.search('d1', file, re.I):
            factor = factors['factor.k27me3']
        else:
            factor = factors['factor.k27me3']
    else:
        factor = input('please input tmm factor:\n')
        try:
            factor = float(factor)
        except:
            exit('cannot tell which TMM factor to use')


    dat = pd.read_csv(file, sep='\t', header=None)
    dat.iloc[:,-1] = round(dat.iloc[:,-1] * factor,2)

    # filter the data
    dat = dat[(dat.iloc[:,-1]>1) | (dat.iloc[:,-1]<-1)]
    dat.to_csv(re.sub('.bedgraph', '.tmmNorm.bedgraph', file), sep='\t', header=False, index=False, mode='w')
    return re.sub('.bedgraph', '.tmmNorm.bedgraph', file)

def tobw(file):
    cmd = 'bedGraphToBigWig {} hg19.clean.sizes.txt {}'.format(file, re.sub('.bedgraph', '.bw', file))
    sp.call(cmd, shell=True)


for i in listdir(sys.argv[1]):
    if re.search('subtract.bedgraph', i):
        sys.stderr.write(i+'\n')
        if re.search('rep\d+_O', i):
            f = main_process(i)
            tobw(f)
        elif re.search('rep\d+_Y', i):
            out = re.sub('.bedgraph', '.tmmNorm.bedgraph', i)
            sp.call('cp {} {}'.format(i, out), shell=True)
            tobw(out)
