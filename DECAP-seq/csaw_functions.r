# Now, this script is used as a function supplemental file.

# this script is used to find different binding regions of ChIP-seq data using csaw
# essentially, it counts number of reads in definded window, filter, and run edgeR
# in experimental design, batch/donor effect was counted along with age difference
# needed input files:
# 1. bam files  -- must have
# 2. black list bed file
# 3. hg19 GTF annotation file

# multi window: combine multiple window using csaw function does not works well,
# this method loss of results. Combine it my way
# use H3 total as control to calculate norm factor using TMM method -- 0.2

library(csaw)
library(dplyr)
library(readxl)
library(edgeR)
library(GenomicFeatures)
library(jsonlite)
library(rtracklayer)
library(tidyr)


calculate_factor <- function(dat, plot, out){
    colnames(dat)[c(1,2)] = c('V1','V2')
    dat = dat[rowSums(dat==0)<2,]  # fitler both are 0
    dat = dat[order(dat[,1]),]
    dat$rank = seq(1:nrow(dat))
    dat$bin_pct = dat$rank/nrow(dat)
    dat$dat1_tag_pct = cumsum(dat[,1])/sum(dat[,1])
    dat$dat2_tag_pct = cumsum(dat[,2])/sum(dat[,2])
    dat$max = abs(dat$dat1_tag_pct-dat$dat2_tag_pct)
    maxpq = dat[which.max(abs(dat$dat1_tag_pct-dat$dat2_tag_pct)),]
    k = maxpq$rank
    factor = sum(dat[,1][1:k])/sum(dat[,2][1:k])
    total_scale = sum(dat[,1])/sum(dat[,2])
    if(plot==TRUE){
		pdf(out)
        plot(dat$bin_pct, dat$dat1_tag_pct, type='line',xlim=c(0,1), col='red', lwd=1)
        points(dat$bin_pct, dat$dat2_tag_pct, type='line', col='black', lwd=1)
        points(dat$bin_pct, dat$max, type='line', lty=2, lwd=0.5)
        abline(v=maxpq$bin_pct,lty=2, lwd=0.5)
		dev.off()
    }
    return(list('ses'=factor, 'sds'=total_scale))
    #
}

workflow <- function(bam.files, width, binned, param=param, pro=pro, gene=gene){
    # the whole workflow from count to DE detect and saving results files
    # the function is ideal for try with different width.
    # binned object should be counted first
    data <- windowCounts(bam.files, width=width, param=param)
    save_count(data, paste0('/Volumes/Data1/ChIP_seq_11_2018/csaw/nonPromoterH3K36me3diff_bin',
                            width,'_count.json'))

    res <- main(data, binned)
    save(res, paste0('/Volumes/Data1/ChIP_seq_11_2018/csaw/nonPromoterH3K36me3diff_bin',
                     width,'_compositionNorm.json'))
    final <- merge_filter(res, pro, gene)
    l=list('count_data'=data, 'res'=res, 'final'=final)
    return(l)
}

myCombine <- function(datlist){
    names(datlist) <- paste0('l_', seq(length(datlist)))
    length_list <- lapply(datlist, nrow)
    ind <- names(sort(unlist(length_list), decreasing = T))
    comb <- datlist[[ind[1]]]
    for(i in 2:length(ind)) {
        sp <- !overlapsAny(makeGRangesFromDataFrame(comb), makeGRangesFromDataFrame(datlist[[ind[i]]]))
        comb <- rbind(comb, datlist[[ind[i]]])
    }
    return(comb)
}

getsig <- function(lis, pvalue=0.01){
    # input is a list, the output from main or main_single function
    # this function will only keep windows with pvalue < 0.01
    print(pvalue)
    lis$tab -> dat
    sig.dat <- dat[dat$PValue<pvalue,]
    s = row.names(sig.dat)
    newdat <- data.frame(separate(data.frame(s), s, c('chr','start','end'), sep='_'))
    return(list('region'=newdat,'tab'=sig.dat))
}

main=function(data, binned){
    normfacs=binned$norm.factors
    # filter by global enrichment
    filter.stat <- filterWindows(data, binned, type="global")
    min.fc <- 3
    keep <- filter.stat$filter > log2(min.fc)
    filtered.data <- data[keep,]
    age=factor(c('Y','O','Y','O'))
    age <- relevel(age, ref="Y")
    donor=factor(c('d1','d1','d3','d3'))
    dat=assay(filtered.data)
    colnames(dat)=c('Y1','O1','Y3','O3')
    temp=data.frame(rowRanges(filtered.data))
    rownames(dat)=paste(temp$seqnames, temp$start, temp$end, sep = '_')
    y <- DGEList(dat, norm.factors=normfacs, group=age, lib.size = filtered.data$totals)
    design <- model.matrix(~donor+age)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit)
    return(list('res'=res, 'filtered.data'=filtered.data, 'normfacs'=normfacs))
}

save_count <- function(data, outfile){
    # save the csaw object directly. This is most convenient way
    save(data, file = outfile)
}

save_count1 <- function(data, outfile){
    dat = assay(data)
    region = data.frame(rowRanges(data))[,1:3]
    dat = cbind(region,dat)
    colnames(dat) = c('chr','start','end',data$bam.files)
    write.table(dat, outfile, quote=F, row.names = F, col.names = T)
}


merge_filter=function(res, pro, gene){
    merged <- mergeWindows(rowRanges(res$filtered.data), tol=width(rowRanges(res$filtered.data))[1],
                           max.width = 5000)
    tabcom <- combineTests(merged$id, res$res$table)
    tabbest <- getBestTest(merged$id, res$res$table)
    tabcom=cbind(data.frame(merged$region), tabcom)
    tabbest=cbind(data.frame(merged$region), tabbest)
    is.sig <- tabcom$FDR <= 0.05
    is.sig1 = tabbest$FDR <= 0.05
    summary(is.sig)
    summary(is.sig1)
    tabcom.sig=tabcom[is.sig,]
    s = overlapsAny(makeGRangesFromDataFrame(tabcom.sig), pro)
    s1 = overlapsAny(makeGRangesFromDataFrame(tabcom.sig), gene)
    final=data.frame(tabcom.sig[!s&s1,])
    tabcom.sig$promoter='No'
    tabcom.sig$promoter[s]='Yes'
    tabcom.sig$genic='No'
    tabcom.sig$genic[s1]='Yes'
    return(list('final'=final, 'allsig'=tabcom.sig))
}

save_edgeR=function(res, outfile){
    # save to json and can be read by fromJSON function
    range=data.frame(rowRanges(res$filtered.data))
    tab=data.frame(res$res$table)
    new=cbind(range, tab)
    normfacs=res$normfacs
    newlist=list('region'=range, 'tab'=tab, 'normfacs'=normfacs)
    write_json(newlist, outfile)
    #jsontest = fromJSON('newlist.json')
}

ma_plot <- function(count_data, ref_num, case_num, normfacs){
    par(mfrow=c(1,1))
    #normfacs=count_data$norm.factors
    bin.ab <- scaledAverage(count_data)
    adjc <- calculateCPM(count_data, use.norm.factors=FALSE)
    print(count_data$norm.factors)
    smoothScatter(bin.ab, adjc[,ref_num]-adjc[,case_num], ylim=c(-6, 6),
                  xlab="Average abundance", ylab="Log-ratio")
    abline(h=log2(normfacs[ref_num]/normfacs[case_num]), col="black")
    par(mfrow=c(1,1))
}


ma_plot_edgeR <- function(count_data, ref_num, case_num, normfacs){
    par(mfrow=c(1,1))
    #normfacs=count_data$norm.factors
    bin.ab <- aveLogCPM(count_data)
    adjc <- cpm(count_data, normalized.lib.sizes=FALSE, log=TRUE)
    print(count_data$norm.factors)
    smoothScatter(bin.ab, adjc[,ref_num]-adjc[,case_num],
                  xlab="Average abundance", ylab="Log-ratio")
    abline(h=log2(normfacs[ref_num]/normfacs[case_num]), col="black")
    par(mfrow=c(1,1))
}






main_single=function(data, binned){
    normfacs=binned$norm.factors
    # filter by global enrichment
    filter.stat <- filterWindows(data, binned, type="global")
    min.fc <- 3
    keep <- filter.stat$filter > log2(min.fc)
    filtered.data <- data[keep,]
    age=factor(c('Y','O'))
    age <- relevel(age, ref="Y")
    # donor=factor(c('d1','d1','d3','d3'))
    dat=assay(filtered.data)
    colnames(dat)=c('Y','O')
    temp=data.frame(rowRanges(filtered.data))
    rownames(dat)=paste(temp$seqnames, temp$start, temp$end, sep = '_')
    y <- DGEList(dat, norm.factors=normfacs, group=age, lib.size = filtered.data$totals)
    design <- model.matrix(~0+age)
    res <- exactTest(y, dispersion=0.05)
    #fit <- glmFit(y, design, dispersion=0.05)
    #res <- glmLRT(fit)
    return(list('res'=res, 'filtered.data'=filtered.data, 'normfacs'=normfacs))
}


# input blacklist. HG19
.backgroud_files <- function(){
    bl='blacklist.bed'
    bl=read.table(bl, sep='\t')
    blacklist=makeGRangesFromDataFrame(bl, seqnames.field = 'V1',
                                       start.field = 'V2',
                                       end.field = 'V3')
    # read sample label
    label=read_excel('sequencing\ libraries\ info.xlsx')
    k4me3=label[grepl('H3K4me3', label$sample),]

    # input bam files and count
    standard.chr <- paste0("chr", c(1:22, "X", "Y"))
    bam.files=list.files(dir)[grepl('.bam.gz$', list.files(dir))]

    # find out h3k4me3
    k4me3.bam.files=bam.files[grepl('H3K4me3', label$sample)]
    h3.bam.files=bam.files[grepl('H3 total', label$sample)]


    # read h3k36me3 data
    dir='/Volumes/Data1/ChIP_seq_11_2018/csaw'
    setwd(dir)
    files=list.files()
    h36count=files[grepl('H3K36me3.*count.*json', files)]

}
