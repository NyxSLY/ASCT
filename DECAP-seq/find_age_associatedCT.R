library(csaw)
library(edgeR)
library(tidyr)
library(dplyr)
library(GenomicFeatures)


# -----------------
# genome annotation
# -----------------
gtff = 'gencode.v19.annotation.gtf'
gtf = makeTxDbFromGFF(gtff, format='gtf')
pro=promoters(gtf,upstream=1000, downstream=1000)
gene=genes(gtf)
# 3'UTR region
utr3o <- threeUTRsByTranscript(gtf)
utr3 <- makeGRangesFromDataFrame(data.frame(utr3o))
# -----------------

# -------------------
# function
# -------------------
# calculate.old = function(dat.c){
    rn = names(dat.c)
    if(!match('logCPM', rn)){
        stop("input data.frame don't have logCPM colnum")
    }

    if(ncol(dat.c)!=6){
        stop("input data.frame should have 6 colnum")
    }

    dat.c$PValue=ppois(dat.c[,5]+1, dat.c[,4]+1, F)
    dat.c$fc=(dat.c[,5]+1)/(dat.c[,4]+1)
    dat.c$logFC = log2(dat.c$fc)
    merged.c <- mergeWindows(makeGRangesFromDataFrame(dat.c),
                             tol=100,
                             max.width = 1000)
    # out.c = data.frame(merged.c$region)[,1:3]
    # write.table(out, 'mergedWindow.bed', quote=F, col.names = F, row.names = F)
    tab.best.c <- getBestTest(merged.c$id, dat.c, by.pval = F)
    tab.best.all.c = cbind(data.frame(merged.c$region), tab.best.c)
    tab.best.all.c$logFDR =-log10(tab.best.all.c$FDR)
    tab.best.all.c[tab.best.all.c$FDR==0,]$logFDR = max(tab.best.all.c[tab.best.all.c$FDR!=0,]$logFDR)+1
    tab.best.all.c=data.frame(tab.best.all.c)
    return(tab.best.all.c)
}

calculate = function(dat.c, reci){
    if(reci==T){
        backup=dat.c
        dat.c[,5] = backup[,4]
        dat.c[,4] = backup[,5]
    }

    rn = names(dat.c)
    if(!match('logCPM', rn)){
        stop("input data.frame don't have logCPM colnum")
    }

    if(ncol(dat.c)!=6){
        stop("input data.frame should have 6 colnum")
    }

    dat.c$PValue=ppois(dat.c[,5]+1, dat.c[,4]+1, F)
    dat.c$fc=(dat.c[,5]+1)/(dat.c[,4]+1)
    dat.c$logFC = log2(dat.c$fc)
    merged.c <- mergeWindows(makeGRangesFromDataFrame(dat.c),
                             tol=100,
                             max.width = 1000)
    # out.c = data.frame(merged.c$region)[,1:3]
    # write.table(out, 'mergedWindow.bed', quote=F, col.names = F, row.names = F)
    tab.best.c <- getBestTest(merged.c$id, dat.c, by.pval = F)
    tab.best.all.c = cbind(data.frame(merged.c$region), tab.best.c)
    tab.best.all.c$logFDR =-log10(tab.best.all.c$FDR)
    tab.best.all.c[tab.best.all.c$FDR==0,]$logFDR = max(tab.best.all.c[tab.best.all.c$FDR!=0,]$logFDR)+1
    tab.best.all.c=data.frame(tab.best.all.c)
    return(tab.best.all.c)
}

count <- function(ori_bam, name){
    standard.chr <- paste0("chr", c(1:22, "X", "Y"))
    # setup parallel env
    multicoreParam <- MulticoreParam(workers = 10)
    # setup parameters
    param <- readParam(minq=30,restrict=standard.chr,
                       dedup=FALSE, BPPARAM=MulticoreParam(workers = 8))

    counts = windowCounts(ori_bam, width=10, bin=TRUE, param=param, ext=NA)
    counts.pure=assay(counts)
    counts.pure1 = cbind(data.frame(rowRanges(counts))[,1:3], data.frame(counts.pure))
    names(counts.pure1)[4:6] = ori_bam
    write.csv(data.frame(counts.pure1),paste0('new_combine/', name, '.count.csv'), quote = F, row.names = F)

    keep = rowSums(counts.pure>10) > 1
    filter.counts=counts.pure[keep,]
    filter.region = rowRanges(counts)[keep,]
    dat = data.frame(filter.counts)
    s = data.frame(filter.region)[,1:3]
    row.names(dat) = paste(s$seqnames, s$start, s$end, sep='_')
    dat1=dat[,c(3,2)]
    age=factor(c('Y','O'))
    age <- relevel(age, ref="Y")
    y <- DGEList(dat1,group=age)
    y <- calcNormFactors(y)
    norm.lib.size = y$samples$lib.size * y$samples$norm.factors
    norm.factor = max(norm.lib.size)/min(norm.lib.size)
    norm.counts = dat1
    colnames(norm.counts) = age
    norm.counts[,which.max(norm.lib.size)] = norm.counts[,which.max(norm.lib.size)]/sqrt(norm.factor)
    norm.counts[,which.min(norm.lib.size)] = norm.counts[,which.min(norm.lib.size)]*sqrt(norm.factor)
    norm.counts = round(norm.counts,2)
    region=separate(data.frame(row.names(norm.counts)),1,c('chr','start','end'), sep='_')
    control.bg=cbind(region, norm.counts[,1])
    treatment.bg=cbind(region, norm.counts[,2])
    write.table(control.bg, paste0('new_combine/', name,'.Y.norm.bg'), quote=F, row.names = F, col.names = F, sep='\t')
    write.table(treatment.bg, paste0('new_combine/', name,'.O.norm.bg'), quote=F, row.names = F, col.names = F, sep='\t')
}
# -------------------


# -------------------------------------
# input data needed and set parameters
# -------------------------------------

# input bam files and count
standard.chr <- paste0("chr", c(1:22, "X", "Y"))

# input blacklist. HG19
bl='/Human_database/hg19/blacklist.bed'
bl=read.table(bl, sep='\t')
blacklist=makeGRangesFromDataFrame(bl, seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3')

# count parameter
multicoreParam <- MulticoreParam(workers = 10)
param <- readParam(minq=30,discard=blacklist, restrict=standard.chr,
                   dedup=FALSE, BPPARAM=MulticoreParam(workers = 8))

# input bams
setwd('star_mapping/filter')
files=list.files()
ori_bam_rev=files[grepl('onlyR1.rev.filtered.bam$', files)]
ori_bam_fwd=files[grepl('onlyR1.fwd.filtered.bam$', files)]
ori_bam = files[grepl('nomulti.bam$', files)]


# -------------------
# main function
# -------------------


main_function = function(ori_bam, fwd, output_dir, reci){
    setwd('star_mapping/filter')
    # output name
    if(fwd==T){
        prefix2='fwd'
    } else{
        prefix2='rev'
    }
    if(reci==T){
        prefix=paste('reciprocal', prefix2, sep='.')
    } else{
        prefix=prefix2
    }

    # use 10k bin to calculate TMM factor
    counts = windowCounts(ori_bam, width=10000, bin=TRUE, param=param, ext=NA)
    age=factor(c('C','O', 'Y'))
    x = data.frame(assay(counts))
    x.filter= dplyr::filter(x, X1>0, X2>0, X3>0)
    y = DGEList(x.filter, group=age)
    y = calcNormFactors(y)

    # do the analysis again
    x = read.csv(count_csv)
    dato = x[,4:6]
    row.names(dato) = paste(x$seqnames, x$start, x$end, sep='_')
    y.new <- DGEList(dato, norm.factors=y$samples$norm.factors, group=age, lib.size = colSums(dato))

    keep = rowSums(cpm(y.new) > 1) > 1
    cpm.filter = cpm(y.new)[keep,]

    # output folder
    setwd(output_dir)

    # ----------------------------------------
    # filter based on pre-called peak -- call peak first between negControl and P15
    # ----------------------------------------
    norm.count.c = data.frame(cpm.filter * max(colSums(dato)) /10^6)
    norm.dat.c = data.frame(norm.count.c[,1],
                            norm.count.c[,2])
    row.names(norm.dat.c) = row.names(norm.count.c)
    region.c=separate(data.frame(row.names(norm.dat.c)),1,c('seqnames','start','end'), sep='_')
    dat.c = cbind(region.c, norm.dat.c)
    colnames(dat.c)[4:5] = c('ctl', 'p15')
    dat.c$logCPM=as.array(log2(rowSums(data.frame(cpm.filter[,c(1,2)]))/2))
    pre_peak = calculate(dat.c, reci=F)

    # P5 preCall
    norm.dat.c = data.frame(norm.count.c[,1],
                            norm.count.c[,3])
    row.names(norm.dat.c) = row.names(norm.count.c)
    region.c=separate(data.frame(row.names(norm.dat.c)),1,c('seqnames','start','end'), sep='_')
    dat.c = cbind(region.c, norm.dat.c)
    colnames(dat.c)[4:5] = c('ctl', 'p5')
    dat.c$logCPM=as.array(log2(rowSums(data.frame(cpm.filter[,c(1,2)]))/2))
    p5_pre_peak = calculate(dat.c, reci=F)

    # call peak
    peak = dplyr::filter(pre_peak, FDR<0.05, logFC>0) # used for preCall method
    peak_p15 = dplyr::filter(pre_peak, FDR<0.05, fc>1.2)
    peak_region_p15 = makeGRangesFromDataFrame(peak)

    peak_p5 = dplyr::filter(p5_pre_peak, FDR<0.05, fc>1.2)
    peak_region_p5 = makeGRangesFromDataFrame(peak_p5)

    #
    s = overlapsAny(peak_region_p15, pro)
    s1 = overlapsAny(peak_region_p15, gene)
    s2 = overlapsAny(peak_region_p15, utr3)
    final_peak_p15=data.frame(peak_region_p15[!s&s1&!s2,])

    s = overlapsAny(peak_region_p5, pro)
    s1 = overlapsAny(peak_region_p5, gene)
    s2 = overlapsAny(peak_region_p5, utr3)
    final_peak_p5=data.frame(peak_region_p5[!s&s1&!s2,])

    write.csv(final_peak_p15, paste(prefix,'P15_peak.csv',sep='.'), quote=F, row.names=F)
    write.csv(final_peak_p5, paste(prefix,'P5_peak.csv', sep='.'), quote=F, row.names=F)

    # make data frame
    norm.dat = data.frame(norm.count.c[,3],
                          norm.count.c[,2])
    row.names(norm.dat) = row.names(norm.count.c)
    region=separate(data.frame(row.names(norm.dat)),1,c('seqnames','start','end'), sep='_')
    dat.new = cbind(region, norm.dat)

    # filter by peak region
    m = overlapsAny(makeGRangesFromDataFrame(dat.new), peak_region_p15) # may have problem, can not find peak_region; changed
    dat.filtered = dat.new[m,]
    dat.filtered$logCPM = log2(rowSums(cpm.filter[,c(3,2)])/2)[m]
    out = calculate(dat.filtered, reci=reci)
    filtered = dplyr::filter(out, fc>1.2, FDR<0.05)

    # filter promoter and 3'UTR
    colnames(filtered) -> name
    colnames(filtered)[1:3] <- name[7:9]
    colnames(filtered)[7:9] <- name[1:3]

    s = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), pro)
    s1 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), gene)
    s2 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), utr3)

    final=data.frame(filtered[!s&s1&!s2,])
    final_precall = filtered[!s&s1&!s2,]
    dim(final)
    write.csv(final, paste(prefix,'final.preCallPeak.csv', sep='.'), quote=F, row.names = F)
    write.table(final[,c(7,8,9,17)], paste(prefix,'final.preCallPeak.bedgraph', sep='.'), quote=F, row.names = F, col.names = F, sep = '\t')
    #-----------------


    # -------------------
    # filter -- globle mode, basically regions with greater 3 fold of medium(negControl) were kept
    # -------------------
    sh = cpm.filter[,1] %>% median()
    keep1 = rowMeans(cpm.filter[,-1])/sh > 3
    cpm.filter1 = cpm.filter[keep1,]

    norm.count = data.frame(cpm.filter1 * max(colSums(dato)) /10^6)
    norm.dat = data.frame(norm.count[,3],
                          norm.count[,2])
    row.names(norm.dat) = row.names(norm.count)
    norm.dat = round(norm.dat,2)

    region=separate(data.frame(row.names(norm.dat)),1,c('seqnames','start','end'), sep='_')
    dat = cbind(region, norm.dat)
    dat$logCPM = as.array(log2(rowSums(data.frame(cpm.filter1[,c(3,2)]))/2))
    colnames(dat)[4:5] = c('p5', 'p15')
    out1 = calculate(dat, reci=reci)
    filtered = dplyr::filter(out1, fc>1.2, FDR<0.05)
    #write.table(out1[,c(1,2,3,13,15,5)], 'all.fwd.result.bed', quote = F, row.names = F, col.names = F)

    # filter promoter and 3'UTR
    colnames(filtered) -> name
    colnames(filtered)[1:3] <- name[7:9]
    colnames(filtered)[7:9] <- name[1:3]

    s = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), pro)
    s1 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), gene)
    s2 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), utr3)

    final=data.frame(filtered[!s&s1&!s2,])
    final_global = filtered[!s&s1&!s2,]
    dim(final)
    write.csv(final, paste(prefix,'final.global3fold.csv', sep='.'), quote=F, row.names = F)
    write.table(final[,c(7,8,9,17)], paste(prefix,'final.global3fold.bedgraph', sep='.'), quote=F, row.names = F, col.names = F, sep = '\t')

    # precall and global overlap
    if(nrow(final_precall)>nrow(final_global)){
        kk = overlapsAny(makeGRangesFromDataFrame(final_precall),
                         makeGRangesFromDataFrame(final_global))
        finall = data.frame(final_precall[kk,])
    } else{
        kk = overlapsAny(makeGRangesFromDataFrame(final_global),
                         makeGRangesFromDataFrame(final_precall))
        finall = data.frame(final_global[kk,])
    }
    write.csv(finall, paste(prefix,'final.2methodOverlapped.csv', sep='.'), quote=F, row.names = F)
    write.table(finall[,c(7,8,9,17)], paste(prefix,'final.2methodOverlapped.bedgraph', sep='.'), quote=F, row.names = F, col.names = F, sep = '\t')


    # generate bedgraphd
    if(reci!=T){
        write.table(dat[,c(1,2,3,4)], paste0(prefix2, '.p5.norm.bedgraph'), quote = F, col.names = F, row.names = F, sep='\t')
        write.table(dat[,c(1,2,3,5)], paste0(prefix2, '.p15.norm.bedgraph'), quote=F, col.names=F, row.names=F, sep='\t')
        system(paste0('/Users/Luyang/Documents/Software/UCSC_tools/bedGraphToBigWig ',paste0(prefix2, '.p5.norm.bedgraph'),
                      ' /Volumes/LACIE/Human_database/hg19/hg19.clean.sizes.txt ',
                      paste0(prefix2, '.p5.norm.bw')))
        system(paste0('/Users/Luyang/Documents/Software/UCSC_tools/bedGraphToBigWig ',paste0(prefix2, '.p15.norm.bedgraph'),
                      ' /Volumes/LACIE/Human_database/hg19/hg19.clean.sizes.txt ',
                      paste0(prefix2, '.p5.norm.bw')))
    }

}

# -------------------
# parameters may need to be changed
# -------------------

fwd = F
output_dir='star_mapping/filter/new_combine/test_precall'
reci = F
#ori_bam = ori_bam_rev

if(fwd == TRUE){
    ori_bam = ori_bam_fwd
    count_csv = 'star_mapping/filter/new_combine/fwd.count.csv'
} else{
    ori_bam = ori_bam_rev
    count_csv = 'star_mapping/filter/new_combine/rev.count.csv'
}


# -------------------
# Dnot need to run again

count(ori_bam_fwd, 'fwd')
count(ori_bam_rev, 'rev')

# --------------------------------------------------------------------


main_function(ori_bam, fwd=fwd, output_dir = output_dir, reci=reci)

# ---------

# ---------------------------------
# make SES subtraction signal track
# ---------------------------------
# factor files are in star_mapping/filter/new_combine/

count_csv='star_mapping/filter/new_combine/fwd.count.csv'
x = read.csv(count_csv)

ses_p15 = data.frame(chr=x$seqnames,
                     start=x$start,
                     end=x$end,
                     ses_p15 = (x$star_P15Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam-
                         x$star_P15_negAligned.sorted.nomulti.onlyR1.fwd.filtered.bam*1.91)*1.05)

ses_p5 = data.frame(chr=x$seqnames,
                    start=x$start,
                    end=x$end,
                    ses_p5 = (x$star_P5Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam-
                        x$star_P15_negAligned.sorted.nomulti.onlyR1.fwd.filtered.bam*1.61)*0.95)

ses_p5[ses_p5$ses_p5<0,]$ses_p5 = 0
ses_p15[ses_p15$ses_p15<0,]$ses_p15 = 0

# rev
count_csv='star_mapping/filter/new_combine/rev.count.csv'
x = read.csv(count_csv)

ses_p15_rev = data.frame(chr=x$seqnames,
                     start=x$start,
                     end=x$end,
                     ses_p15 = (x$star_P15Aligned.sorted.nomulti.onlyR1.rev.filtered.bam-
                                    x$star_P15_negAligned.sorted.nomulti.onlyR1.rev.filtered.bam*1.88)*1.05)

ses_p5_rev = data.frame(chr=x$seqnames,
                    start=x$start,
                    end=x$end,
                    ses_p5 = (x$star_P5Aligned.sorted.nomulti.onlyR1.rev.filtered.bam-
                                  x$star_P15_negAligned.sorted.nomulti.onlyR1.rev.filtered.bam*1.59)*0.95)

ses_p5_rev[ses_p5_rev$ses_p5<0,]$ses_p5 = 0
ses_p15_rev[ses_p15_rev$ses_p15<0,]$ses_p15 = 0

setwd('star_mapping/filter/new_combine/SES_TMM')
write.table(ses_p5, 'p5.fwd.SES_TMM.bin10.bedgraph', sep='\t', row.names = F,
            col.names = F, quote = F)
write.table(ses_p15, 'p15.fwd.SES_TMM.bin10.bedgraph', sep='\t', row.names = F,
            col.names = F, quote = F)

write.table(ses_p5_rev, 'p5.rev.SES_TMM.bin10.bedgraph', sep='\t', row.names = F,
            col.names = F, quote = F)
write.table(ses_p15_rev, 'p15.rev.SES_TMM.bin10.bedgraph', sep='\t', row.names = F,
            col.names = F, quote = F)


# ---------------------------------
# filter promoter and 3'UTR region -- already included. Here is just for backup
# ---------------------------------

# # find ones donot overlap with the last exon
# # xx = data.frame(exonsBy(gtf, by='tx', use.names=T))
# # b = tapply(row.names(xx), xx$group_name, function(x) tail(x, n=1) )
# # last_exon = makeGRangesFromDataFrame(xx[b,])
#
# colnames(filtered) -> name
# colnames(filtered)[1:3] <- name[7:9]
# colnames(filtered)[7:9] <- name[1:3]
#
# # filter out mix
# # tab.com <- combineTests(merged$id, dat)
# # tab.com = cbind(data.frame(merged$region), tab.com)
# #
# # mix = tab.com[grepl('mix', tab.com$direction),]
# # m = overlapsAny(makeGRangesFromDataFrame(mix), makeGRangesFromDataFrame(final))
#
# s = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), pro)
# s1 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), gene)
# s2 = overlapsAny(makeGRangesFromDataFrame(filtered[,c(7:9)]), utr3)
#
# final=data.frame(filtered[!s&s1&!s2,])
#
# write.table(final[,c(7,8,9,17)], 'fwd.final.preCallPeak.bedgraph', quote=F, row.names = F, col.names = F)
# write.table(final[,c(7,8,9,17)], 'fwd.final.global3fold.bedgraph', quote=F, row.names = F, col.names = F)

# ---------------
