# SES normalization
# Normalization, bias correction, and peak calling for ChIP-seq
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342857/

# a package using this method
# Differential peak calling of ChIP-seq signals with replicates with THOR
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5175345/#B42

# -----------------------------------------------
# this is used for TL-seq. Calculate SES factor then apply it, then calculate
# TMM factor, then apply ---> generate final bedgraph file

# width = 2kb. Donot change this unless you know what you are doing
# things may need to change:
# 1. current setup is used for human hg19, including blacklist and chromosome list
# 2. counting readParam setup: if pe; if dedup; if run on windows, need to change multiprocess platform
# 3. input bams and working dir
# -----------------------------------------------

source('csaw_functions.r')
library(csaw)
library(readxl)
standard.chr <- paste0("chr", c(1:22, "X", "Y"))

bl='/Volumes/LACIE/Human_database/hg19/blacklist.bed'
bl=read.table(bl, sep='\t')
blacklist=makeGRangesFromDataFrame(bl, seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3')

# setup parallel env
#multicoreParam <- MulticoreParam(workers = 10)
# setup parameters

# this is for TL-seq
param <- readParam(minq=30,discard=blacklist, restrict=standard.chr,
                   max.frag=1000, pe="none", dedup=FALSE, BPPARAM=MulticoreParam(workers = 8))
# -----

# -----
# input bam
# -----
setwd('star_mapping/filter')
files=list.files()
#ori_bam_rev=files[grepl('onlyR1.rev.filtered.bam$', files)]
#ori_bam_fwd=files[grepl('onlyR1.fwd.filtered.bam$', files)]

files[grepl('P15_negAligned.*onlyR1.*.filtered.bam$', files)] -> control_bam
files[grepl('P15Aligned.*onlyR1.*.filtered.bam$', files)] -> p15_bam
files[grepl('P5Aligned.*onlyR1.*.filtered.bam$', files)] -> p5_bam


# ------
# count
# ------

width = 2000
counts.control = windowCounts(control_bam, width=width, bin=TRUE, param=param, filter=0)
counts.p15 = windowCounts(p15_bam, width=width, bin=TRUE, param=param, filter=0)
counts.p5 = windowCounts(p5_bam, width=width, bin=TRUE, param=param, filter=0)

# --------
# process data
# --------

region = data.frame(rowRanges(counts.control))[,1:3]
counts.control.df = data.frame(assay(counts.control))
colnames(counts.control.df) = c('fwd','rev')
counts.p15.df = data.frame(assay(counts.p15))
colnames(counts.p15.df) = c('fwd','rev')
counts.p5.df = data.frame(assay(counts.p5))
colnames(counts.p5) = c('fwd','rev')


#------------
# p5

l = list()
for(i in 1:2){
    dat = data.frame(counts.p5.df[,i], counts.control.df[,i])
    factor = calculate_factor(dat, plot=F)
    l[[i]] = factor
}

f5 = data.frame(unlist(sapply(l, "[", 1)), unlist(sapply(l, "[", 2)))
colnames(f5) = c('ses','sds')
rownames(f5) = c('fwd','rev')
write.table(f5, 'star_mapping/filter/new_combine/p5.SESfactor.txt', sep='\t')


#------------
# p15

l = list()
for(i in 1:2){
    dat = data.frame(counts.p15.df[,i], counts.control.df[,i])
    factor = calculate_factor(dat, plot=F)
    l[[i]] = factor
}

f15 = data.frame(unlist(sapply(l, "[", 1)), unlist(sapply(l, "[", 2)))
colnames(f15) = c('ses','sds')
rownames(f15) = c('fwd','rev')
write.table(f15, 'star_mapping/filter/new_combine/p15.SESfactor.txt', sep='\t')


#------------
# substract
#------------

p5.sub = list()
for(i in 1:ncol(counts.p5.df)){
    sub = counts.p5.df[,i] - (f5$ses[i]*counts.control.df[,i])
    p5.sub[[i]] = sub
}

p5.sub = data.frame(p5.sub)
colnames(p5.sub) = colnames(counts.p5.df)

p15.sub = list()
for(i in 1:ncol(counts.p15.df)){
    sub = counts.p15.df[,i] - (f15$ses[i]*counts.control.df[,i])
    p15.sub[[i]] = sub
}

p15.sub = data.frame(p15.sub)
colnames(p5.sub) = colnames(counts.p5.df)

dat = data.frame(p5=rowSums(p5.sub), p15=rowSums(p15.sub))

dat[dat<0] = 0
keep=rowSums(dat>0)>0
dat = dat[keep,]

y = DGEList(dat)
y <- calcNormFactors(y)

tmm_factor = y$samples$lib.size*y$samples$norm.factors/mean(y$samples$lib.size*y$samples$norm.factors)

data.frame(TMM_factor=tmm_factor, row.names = c('Y','O')) -> tmm

write.table(tmm, 'tmm_afterSES.txt')
