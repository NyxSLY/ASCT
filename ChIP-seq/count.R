library(csaw)
library(readxl)
# source("/Users/Luyang/scripts/R/csaw_home_template.r")


#standard.chr <- paste0("chr", c(1:22, "X", "Y"))
#mouse.chr <- c(1:19, "X", "Y")
# input blacklist. HG19
#bl <- "/Volumes/LACIE/Human_database/hg19/blacklist.bed"
#bl <- read.table(bl, sep = "\t")
#blacklist <- makeGRangesFromDataFrame(bl,
#    seqnames.field = "V1",
#    start.field = "V2",
#    end.field = "V3"
#)

# setup parallel env
#multicoreParam <- MulticoreParam(workers = 10)
# setup parameters

# this is for PE ChIP-seq, not for RNA-seq or TL-seq

#param <- readParam(
#    minq = 30, discard = blacklist, restrict = standard.chr,
#    max.frag = 1000, pe = "both", dedup = TRUE, BPPARAM = MulticoreParam(workers = 8)
#)


.config = function(genome='hg19', dedup=TRUE, PE=TRUE){
    if(genome == 'hg19'){
        blacklist = '/Users/Luyang/scripts/R/csaw/blacklist/hg19_blacklist.bed'
        blacklist = readbed(blacklist)
        standard.chr <- paste0("chr", c(1:22, "X", "Y"))
    }
    if(PE==TRUE){
        param <- readParam(minq=30,discard=blacklist, restrict=standard.chr,
                   max.frag=1000, pe="both", dedup=dedup, BPPARAM=MulticoreParam(workers = 8))
    } else{
        #param <- readParam(minq=30,discard=blacklist, restrict=standard.chr,
        #           max.frag=1000, pe="None", dedup=dedup, BPPARAM=MulticoreParam(workers = 8))??
    }
    return(param)
}


count_hist <- function(bam.files, prefix, dir, width = 2000) {
    param = .config()
    if (width >= 1000) {
        width_s <- paste0(as.character(width / 1000), "k")
    } else {
        width_s <- width
    }

    outfile <- file.path(dir, paste0(prefix, "_count.bin", width_s, ".rds"))
    if (!file.exists(outfile)) {
        # print(123)
        counts <- windowCounts(bam.files, width = width, bin = TRUE, param = param, filter = 0)
        saveRDS(counts, file = outfile)
    }
    else {
        print(paste0("file already exists: ", outfile))
    }
    return(outfile)
}

get_SES <- function(counts, counts.ctl, outdir, out, sample_names, plot=TRUE) {
    region <- data.frame(rowRanges(counts.ctl))[, 1:3]
    counts.ctl.df <- data.frame(assay(counts.ctl))
    counts.df <- data.frame(assay(counts))
    colnames(counts.df) <- sample_names

    l <- list()
    for (i in 1:ncol(counts.df)) {
        dat <- data.frame(counts.df[, i], counts.ctl.df[, i])
        factor <- calculate_factor(dat,
            plot = TRUE,
            out = paste0(
                outdir, "/",
                out, ".SES2kfactor", colnames(counts.df)[i], ".pdf"
            )
        )
        l[[i]] <- factor
    }

    f <- data.frame(unlist(sapply(l, "[", 1)), unlist(sapply(l, "[", 2)))
    colnames(f) <- c("ses", "sds")
    rownames(f) <- sample_names
    write.table(f, paste0(outdir, "/", out, ".SES2kfactor.txt"), sep = "\t")
    return(paste0(outdir, "/", out, ".SES2kfactor.txt"))
}

get_tmm <- function(counts, counts.ctl, factor_file, outdir, out, sample_names) {
    region <- data.frame(rowRanges(counts.ctl))[, 1:3]
    counts.ctl.df <- data.frame(assay(counts.ctl))
    counts <- data.frame(assay(counts))
    colnames(counts) <- sample_names
    factor <- read.table(factor_file)
    new_sub <- list()
    for (i in 1:ncol(counts)) {
        sub <- counts[, i] - (factor$ses[i] * counts.ctl.df[, i])
        ratio <- (counts[, i] + 1) / (factor$ses[i] * counts.ctl.df[, i] + 1)
        new_sub[[i]] <- sub
    }
    new_sub <- data.frame(new_sub)
    colnames(new_sub) <- colnames(counts)
    new_sub[new_sub < 0] <- 0
    keep <- rowSums(new_sub > 0) > 0
    new_sub <- new_sub[keep, ]
    y <- DGEList(new_sub)
    y <- calcNormFactors(y)
    norm.lib.size <- y$samples$lib.size * y$samples$norm.factors
    factor <- norm.lib.size[1] / norm.lib.size[2]
    o <- data.frame(list(factor = factor))
    write.table(o, paste0(outdir, "/", out, ".TMM2kfactor.txt"),
        row.names = F, sep = "\t"
    )

    write.table(data.frame(y$samples), paste0(outdir, "/", out, ".TMM_sampleInfo.txt"),
        row.names = F, sep = "\t"
    )

    return(paste0(outdir, "/", out, ".TMM_sampleInfo.txt"))
}

calculate <- function(dat.c, reci) {
    # calculate possion pvalue
    if (reci == T) {
        backup <- dat.c
        dat.c[, 5] <- backup[, 4]
        dat.c[, 4] <- backup[, 5]
    }

    rn <- names(dat.c)
    if (!match("logCPM", rn)) {
        stop("input data.frame don't have logCPM colnum")
    }

    if (ncol(dat.c) != 6) {
        stop("input data.frame should have 6 colnum")
    }

    dat.c$PValue <- ppois(dat.c[, 5] + 1, dat.c[, 4] + 1, F)
    dat.c$fc <- (dat.c[, 5] + 1) / (dat.c[, 4] + 1)
    dat.c$logFC <- log2(dat.c$fc)
    merged.c <- mergeWindows(makeGRangesFromDataFrame(dat.c),
        tol = width(makeGRangesFromDataFrame(dat.c))[1],
        max.width = max(width(makeGRangesFromDataFrame(dat.c))[1] * 10, 5000)
    )
    # out.c = data.frame(merged.c$region)[,1:3]
    # write.table(out, 'mergedWindow.bed', quote=F, col.names = F, row.names = F)

    # combineTests merge
    # tabcom <- combineTests(merged.c$id, dat.c)
    # tabcom.all.c = cbind(data.frame(merged.c$region), tabcom)
    # tabcom.all.c$logFDR =-log10(tabcom.all.c$FDR)
    # tryCatch({tabcom.all.c[tabcom.all.c$FDR==0,]$logFDR = max(tabcom.all.c[tabcom.all.c$FDR!=0,]$logFDR)+1},
    #         warning = function(w) {},
    #         error = function(e){},
    #         finally = {}
    # )
    # tabcom.all.c=data.frame(tabcom.all.c)
    # return(tabcom.all.c)

    # getBest merge
    tab.best.c <- getBestTest(merged.c$id, dat.c, by.pval = T)
    tab.best.all.c <- cbind(data.frame(merged.c$region), tab.best.c)
    tab.best.all.c$logFDR <- -log10(tab.best.all.c$FDR)
    tryCatch({
        tab.best.all.c[tab.best.all.c$FDR == 0, ]$logFDR <- max(tab.best.all.c[tab.best.all.c$FDR != 0, ]$logFDR) + 1
    },
    warning = function(w) {},
    error = function(e) {},
    finally = {}
    )
    tab.best.all.c <- data.frame(tab.best.all.c)
    return(tab.best.all.c)

    # no merge
    # dat.c$FDR <- p.adjust(dat.c$PValue, "fdr")
    # return(dat.c)
}

calculate_merge <- function(dat.c, reci) {
    # merge the window first, then calculate possion pvalue based on merged value
    if (reci == T) {
        backup <- dat.c
        dat.c[, 5] <- backup[, 4]
        dat.c[, 4] <- backup[, 5]
    }

    rn <- names(dat.c)
    if (!match("logCPM", rn)) {
        stop("input data.frame don't have logCPM colnum")
    }

    if (ncol(dat.c) != 6) {
        stop("input data.frame should have 6 colnum")
    }

    merged.c <- mergeWindows(makeGRangesFromDataFrame(dat.c),
        tol = width(makeGRangesFromDataFrame(dat.c))[1],
        max.width = max(width(makeGRangesFromDataFrame(dat.c))[1] * 10, 5000)
    )

    input <- csaw:::.check_test_inputs(merged.c$id, dat.c, weight = NULL)
    input$tab$group <- input$ids
    merged.dat <- data.frame(merged.c$region)
    merged.dat$new1 <- tapply(input$tab[, 4], input$tab$group, sum)
    merged.dat$new2 <- tapply(input$tab[, 5], input$tab$group, sum)
    colnames(merged.dat)[c(6, 7)] <- colnames(input$tab)[c(4, 5)]

    merged.dat$PValue <- ppois(merged.dat[, 7] + 1, merged.dat[, 6] + 1, F)
    merged.dat$fc <- (merged.dat[, 7] + 1) / (merged.dat[, 6] + 1)
    merged.dat$logFC <- log2(merged.dat$fc)

    merged.dat$FDR <- p.adjust(merged.dat$PValue, "fdr")
    return(merged.dat)
}

calculate_nomerge <- function(dat.c, reci) {
    # calculate possion pvalue
    if (reci == T) {
        backup <- dat.c
        dat.c[, 5] <- backup[, 4]
        dat.c[, 4] <- backup[, 5]
    }

    rn <- names(dat.c)
    if (!match("logCPM", rn)) {
        stop("input data.frame don't have logCPM colnum")
    }

    if (ncol(dat.c) != 6) {
        stop("input data.frame should have 6 colnum")
    }

    dat.c$PValue <- ppois(dat.c[, 5] + 1, dat.c[, 4] + 1, F)
    dat.c$fc <- (dat.c[, 5] + 1) / (dat.c[, 4] + 1)
    dat.c$logFC <- log2(dat.c$fc)
    merged.c <- mergeWindows(makeGRangesFromDataFrame(dat.c),
        tol = width(makeGRangesFromDataFrame(dat.c))[1],
        max.width = max(width(makeGRangesFromDataFrame(dat.c))[1] * 10, 5000)
    )

    # no merge
    dat.c$FDR <- p.adjust(dat.c$PValue, "fdr")
    return(dat.c)
}

myCombine <- function(datlist){
    #datlist = sapply(datlist, '[', 1)
    names(datlist) <- paste0('l_', seq(length(datlist)))
    length_list <- lapply(datlist, nrow)
    ind <- names(sort(unlist(length_list), decreasing = T))
    comb <- datlist[[ind[1]]]
    for(i in 2:length(ind)) {
        sp <- !overlapsAny(makeGRangesFromDataFrame(datlist[[ind[i]]]), makeGRangesFromDataFrame(comb))
        comb <- rbind(comb, datlist[[ind[i]]][sp,])
    }
    return(comb)
}

find_diffpeak <- function(file1, file2, ses_factor, tmm_sample_info, sample_name,
                          peak_region, mode, reci=F){
    print(mode)
    load(file1)
    case = counts
    load(file2)
    control = counts
    case.df = data.frame(assay(case))
    region = data.frame(rowRanges(case))[,1:3]
    colnames(case.df) = sample_name
    control.df = data.frame(assay(control))
    new = list()
    new_sub = list()
    ses_factor = read.table(ses_factor)
    for(i in 1:ncol(case.df)){
        sub = case.df[,i] - (ses_factor$ses[i]*control.df[,i])
        ratio = (case.df[,i]+1) / (ses_factor$ses[i]*control.df[,i]+1)
        new[[i]] = ratio
        new_sub[[i]] = sub
    }

    # TMM factor
    new = data.frame(new)
    new_sub = data.frame(new_sub)
    colnames(new) = colnames(case.df)
    colnames(new_sub) = colnames(case.df)
    #keep = rowSums(new_sub>0)>0
    new_sub[new_sub<0] = 0
    #keep=rowSums(new_sub>0)>0
    # use subtraction for following analysis
    #new_sub = new_sub[keep,]
    group = sample_names
    sample_info = read.csv(tmm_sample_info, sep='\t')
    y = DGEList(new_sub, lib.size=sample_info$lib.size, norm.factors=sample_info$norm.factors,
                group=group)  # need to change if want to make function

    if(mode=='median'){
        keep1 = rowSums(cpm(y)>0)>length(sample_name)/2
        yy = y[keep1,]
        keep2 = rowSums(cpm(y)>median(cpm(yy)))>length(sample_name)/2
        y = y[keep2,]
    } else if(mode=='upper'){
        keep1 = rowSums(cpm(y)>0)>length(sample_name)/2
        yy = y[keep1,]
        keep2 = rowSums(cpm(y)>quantile(cpm(yy))[4])>length(sample_name)/2
        y = y[keep2,]
    } else if(mode=='fpkm'){
        keep2 = rowSums(cpm(y)/(getWidths(case)[1]/1000)>1)>length(sample_name)/2
        y = y[keep2,]
    } else if(mode=='cpm'){
        keep2 = rowSums(cpm(y)>1)>length(sample_name)/2
        y = y[keep2,]
    } else{
        stop('filter mode is not median!')
    }
    region.c = region[keep2,]
    # only used for 2 samples. Will use poisson expression
    cpm.filter = cpm(y)

    # sh = cpm.filter[,2] %>% median()
    # keep1 = rowMeans(cpm.filter)/sh > 3
    # cpm.filter1 = cpm.filter[keep1,]
    #
    # norm.count = data.frame(cpm.filter1 * max(sample_info$lib.size)/10^6)
    # norm.dat = round(norm.count,2)
    # region.c2 = region.c[keep1,]
    #
    # dat = cbind(region.c2, norm.dat)
    # dat$logCPM = as.array(log2(rowSums(cpm.filter1))/2)
    # colnames(dat)[4:5] = sample_name
    # out1 = calculate(dat, reci=reci)
    # filtered = dplyr::filter(out1, fc>1.2, FDR<0.05)

    norm.count.c = data.frame(cpm.filter * max(sample_info$lib.size) /10^6)  # as long as it's a constant, use which libsize doesn't matter
    norm.dat.c = norm.count.c
    dat.c = cbind(region.c, norm.dat.c)
    colnames(dat.c)[4:5] = sample_names
    if(reci==T){
        dat.c = dat.c[,c(1,2,3,5,4)]
    }
    m = overlapsAny(makeGRangesFromDataFrame(dat.c), peak_region)
    dat.c = dat.c[m,]
    dat.c$logCPM=as.array(log2(rowSums(cpm.filter[,c(1,2)])/2))[m]
    pre_peak = calculate(dat.c, reci=F)
    diffpeak = dplyr::filter(pre_peak, FDR<0.05, fc>1) # may change this
    return(diffpeak)
}

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
