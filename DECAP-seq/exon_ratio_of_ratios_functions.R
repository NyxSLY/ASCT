# use tapply
library(ggplot2)
library(tidyr)

getratio = function(dat, n, ifSecBase){
    if(length(dat)<3){
        return(NA)  # must have at least 3 exons
    }
    if(length(dat)<n){
        return(NA)
    }
    
    if(n <0){
        n = -n
        if(ifSecBase == F){
            ratio = dat[length(dat)-n+1]/dat[1]
        } else{
            ratio = dat[length(dat)-n+1]/dat[2]
        }
    } else{
        # ratio = (dat[n]+1)/(dat[1]+1)
        if(ifSecBase == F){
            ratio = dat[n]/dat[1]
        } else{
            ratio = dat[n]/dat[2]
        }
    }
    return(ratio)
}

# get the value of total expression and expression of the exon, to feed 
# chi's squre test
getvalue = function(dat, n){
    if(length(dat)<3){
        return(NA)  # must have at least 3 exons
    }
    if(length(dat)<n){
        return(NA)
    }
    
    if(n <0){
        n = -n
        # add 1 for presudo count
        l = list(exon=dat[length(dat)-n+1]+1,other=sum(dat)-(dat[length(dat)-n+1])+1)
    } else{
        l = list(exon=dat[n]+1,other=sum(dat)-(exon=dat[n])+1)
    }
    return(l)
}

getratio_other_exons=function(dat, anno, ifSecBase){
    # should be intermediate exons. Remove the last exon
    name_list = anno$normalized_run
    s = tapply(dat$Length, dat$tx_name, length)
    keep = s>=3
    #dat = dat[keep,]
    if(ifSecBase == F){
        length_sum=tapply(dat$Length, dat$tx_name, function(x) sum(x[2:(length(x)-1)]))
        }
    else{
        length_sum=tapply(dat$Length, dat$tx_name, function(x) sum(x[3:(length(x)-1)]))
        }
    new = dat[,c(1:6)]
    n = nrow(anno[anno$group=='Y',])
    m = nrow(anno)
    #n = length(name_list)/2
    
    if(length(name_list)>2){
        new$y = rowMeans(data.frame(dat[,c(7:(7+n-1))]))
        new$o = rowMeans(data.frame(dat[,c((7+n):(7+n+m-1))]))
    } else{
        new$y = dat[,7]
        new$o = dat[,8]
    }
    new$group_name = dat$group_name
    new$TPM = dat$TPM
    if(ifSecBase == F){
        y = tapply(new$y, new$tx_name, function(x) sum(x[2:(length(x)-1)]))
        o = tapply(new$o, new$tx_name, function(x) sum(x[2:(length(x)-1)]))
    }
    else{
        y = tapply(new$y, new$tx_name, function(x) sum(x[3:(length(x)-1)]))
        o = tapply(new$o, new$tx_name, function(x) sum(x[3:(length(x)-1)]))
    }
    
    length_sum = length_sum[keep]
    y=y[keep]
    o=o[keep]
    if(ifSecBase == F){
        y_ratio = y/length_sum/(tapply(new$y, new$tx_name, function(x) x[1])[keep]/tapply(new$Length, new$tx_name, function(x) x[1])[keep])
        o_ratio = o/length_sum/(tapply(new$o, new$tx_name, function(x) x[1])[keep]/tapply(new$Length, new$tx_name, function(x) x[1])[keep])
    }
    else{
        y_ratio = y/length_sum/(tapply(new$y, new$tx_name, function(x) x[2])[keep]/tapply(new$Length, new$tx_name, function(x) x[2])[keep])
        o_ratio = o/length_sum/(tapply(new$o, new$tx_name, function(x) x[2])[keep]/tapply(new$Length, new$tx_name, function(x) x[2])[keep])
    }
    
    k1 = is.na(y_ratio)
    k2 = is.na(o_ratio)
    k3 = is.infinite(y_ratio)
    k4 = is.infinite(o_ratio)
    k5 = y_ratio==0
    k6 = o_ratio==0
    y_ratio = y_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    o_ratio = o_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    # filter gene with < 3 exons
    return(data.frame(y=y_ratio, o=o_ratio))
}  # validated

boxp = function(dat, n, ifSecBase){
    exon.p5 = tapply(dat$P5.nc, dat$Name, getratio, n, ifSecBase)
    exon.p15 = tapply(dat$P15.nc, dat$Name, getratio, n, ifSecBase)
    k1=is.na(exon.p5)
    k2=is.infinite(exon.p5)
    k3 = exon.p5==0
    kk1=is.na(exon.p15)
    kk2=is.infinite(exon.p15)
    kk3 = exon.p15==0
    exon.p5 = exon.p5[!k1&!k2&!kk1&!kk2&!kk3&!k3]
    exon.p15 = exon.p15[!k1&!k2&!kk1&!kk2&!kk3&!k3]
    return(list(p5=exon.p5, p15=exon.p15))
}

ratio_boxp_cal = function(dat,n){
    exon.p5 = tapply(dat$P5.nc, dat$Name, getratio, n)
    exon.p15 = tapply(dat$P15.nc, dat$Name, getratio, n)
    k1=is.na(exon.p5)
    k2=is.infinite(exon.p5)
    k3 = exon.p5==0
    kk1=is.na(exon.p15)
    kk2=is.infinite(exon.p15)
    kk3 = exon.p15==0
    exon.p5 = exon.p5[!k1&!k2&!kk1&!kk2&!kk3&!k3]
    exon.p15 = exon.p15[!k1&!k2&!kk1&!kk2&!kk3&!k3]
    return(log2(exon.p15)-log2(exon.p5))
}

boxplot_filter = function(x){
    l = boxplot.stats(x$value)$stats[c(1, 5)]
    if(l[1]<0){
        sign1 = -1
    }
    else{
        sign1 = 1
    }
    if(l[2]<0){
        sign2 = -1
    }
    else{
        sign2 = 1
    }
        
    if(nrow(x[x$value>l[2],])>0){
        x[x$value>l[2],]$value=abs(l[2])*1.2*sign2
    }
    if(nrow(x[x$value<l[1],])>0){
        x[x$value<l[1],]$value=abs(l[1])*1.2*sign1
    }
    return(x)
}

# ggplot2 boxplot
ggplot_boxplot1 = function(dat, name){
    dat = data.frame(log2(dat))
    dat.m = gather(dat)
    dat.m = boxplot_filter(dat.m)
    dat.m$key = relevel(as.factor(dat.m$key), ref='y')
    png(paste0(name, '_boxplot.png'), width=400, height = 700, units = 'px')
    s = ggplot(dat.m, aes(x=key, y=value, fill=key)) +
        geom_boxplot(outlier.shape=NA) + 
        scale_fill_manual(values=c("#F6402A", "#4F4DB3")) +
        theme_minimal()
    print(s)
    dev.off()
}

ggplot_boxplot = function(dat, name){
    dat = data.frame(log2(dat))
    dat.m = gather(dat)
    dat.m = boxplot_filter(dat.m)
    dat.m$key = relevel(as.factor(dat.m$key), ref='y')
    pdf(paste0(name, '_boxplot.pdf'), width=4, height = 7)
    s = ggplot(dat.m, aes(x=key, y=value, fill=key)) +
        geom_boxplot(outlier.shape=NA, lwd=1) + 
        scale_fill_manual(values=c("#F6402A", "#4F4DB3")) +
        theme_minimal()+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=18, color="black", family="ArialMT"),
            axis.ticks.y = element_line(color="black", size = 1),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            axis.ticks.length.y = unit(0.25, "cm"),
            #legend.title = element_blank(),
            legend.position = "none"
        ) +
        scale_y_continuous(breaks=seq(-2,5,1))
    print(s)
    dev.off()
}

ggplot_violin = function(dat, name, save){
    dat = data.frame(log2(dat))
    dat.m = gather(dat)
    dat.m = boxplot_filter(dat.m)
    dat.m$key = relevel(as.factor(dat.m$key), ref='y')
    s = ggplot(dat.m, aes(x=key, y=value, fill=key)) +
        geom_violin(outlier.shape=NA, lwd=1) + 
        scale_fill_manual(values=c("#F6402A", "#4F4DB3")) +
        theme_minimal()+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=18),
            axis.ticks.y = element_line(color="black", size = 1),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            #legend.title = element_blank(),
            legend.position = "none"
        ) +
        scale_y_continuous(breaks=seq(-2,5,1))
    if(save==F){
        print(s)
    }
    else{
        png(paste0(name, '_boxplot.png'), width=400, height = 700, units = 'px')
        print(s)
        dev.off()
    }

}

ggplot_hist = function(dat, name, save, lim=1.5, breaks=50){
    # log2 ratio of ratios  - histplot
    dat = data.frame(log2(dat))
    dat = data.frame(logratio = dat$o - dat$y)
    a = boxplot.stats(dat$logratio)$stats[c(1, 5)]
    wid = (a[2]-a[1])/breaks
    p = ggplot(dat, aes(x=logratio)) +
        geom_histogram(binwidth=wid, fill="dark blue", colour="white", alpha=0.9) +
        xlim(-lim,lim) +
        geom_vline(xintercept=0, linetype="dashed", color = "white", lwd=0.8)  + 
        theme_minimal()+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color="black", size = 1),
            axis.text.x = element_text(size=18, face='bold', color='black'),
            axis.text.y = element_text(size=18, face='bold', color='black'),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            #legend.title = element_blank(),
            legend.position = "none"
        )
    
    if(save==F){
        print(p)
    }
    else{
        pdf(paste0(name, '_histogram.pdf'), width=10, height = 5.5)
        print(p)
        dev.off()
    }
}



# for quantile TPM
group_ratio_plot = function(dat, name, name_list, ifSecBase){
    dat.fwd = dplyr::filter(dat, Strand=='+')
    dat.rev = dplyr::filter(dat, Strand=='-')
    dat.fwd = arrange(dat.fwd, tx_name, Start)
    dat.rev = arrange(dat.rev, tx_name, desc(Start))
    dat = rbind(dat.fwd, dat.rev)
    
    r2 = cal(name_list, dat, 2, ifSecBase)
    r3 = cal(name_list, dat, 3, ifSecBase)
    r4 = cal(name_list, dat, 4, ifSecBase)
    rl = cal(name_list, dat, -1, ifSecBase)
    other = getratio_other_exons(dat, name_list)
    
    if(ifSecBase==F){
        plot_dat = rbind(gather(data.frame(second=log2(r2$o)-log2(r2$y))),
                         gather(data.frame(third=log2(r3$o)-log2(r3$y))),
                         gather(data.frame(fourth=log2(r4$o)-log2(r4$y))),
                         gather(data.frame(last=log2(rl$o)-log2(rl$y))),
                         gather(data.frame(others=log2(other$o)-log2(other$y))))
        
        plot_dat$key = factor(plot_dat$key, levels=c('second', 'third', 'fourth', 'last', 'others'))
    }else {
        plot_dat = rbind(gather(data.frame(second=log2(r2$o)-log2(r2$y))),
                         gather(data.frame(third=log2(r3$o)-log2(r3$y))),
                         gather(data.frame(fourth=log2(r4$o)-log2(r4$y))),
                         gather(data.frame(fifth=log2(r5$o)-log2(r5$y))),
                         gather(data.frame(last=log2(rl$o)-log2(rl$y))))
        
        plot_dat$key = factor(plot_dat$key, levels=c('second', 'third', 'fourth', 'fifth', 'last'))
    }

    plot_dat = boxplot_filter(plot_dat)
    png(paste0(name, '_boxplot.png'), width=10, height = 5.5, units = 'in', res=1200)
    p = ggplot(plot_dat, aes(x=key, y=value, fill=value)) +
        geom_boxplot(outlier.shape=NA, fill="#FDF6E4", size=1) + 
        theme_minimal() + theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=18, face="bold", color='black'),
            axis.ticks.y = element_line(color="black", size = 1),
            axis.ticks.length=unit(.25, "cm"),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            #legend.title = element_blank(),
            legend.position = "none"
        ) +
        xlab('')
    
    print(p)
    dev.off()
}

cal = function(anno, dat, exon_num, ifSecBase){
    s = tapply(dat$Length, dat$tx_name, length)
    keep = s>=3
    dat = dat[keep,]

    name_list = anno$normalized_run
    
    new = lapply(name_list, function(x) dat[,x])
    new = data.frame(new)
    colnames(new) = name_list
    new_y = new[, anno[anno$group=='Y',]$normalized_run]
    new_o = new[, anno[anno$group=='O',]$normalized_run]
    
    y = rowMeans(data.frame(new_y))
    o = rowMeans(data.frame(new_o))

    new_dat = cbind(dat[,1:6], y=y, o=o)
    y_ratio = tapply(new_dat$y, new_dat$tx_name, getratio, exon_num, ifSecBase)
    o_ratio = tapply(new_dat$o, new_dat$tx_name, getratio, exon_num, ifSecBase)
    
    k1 = is.na(y_ratio)
    k2 = is.na(o_ratio)
    k3 = is.infinite(y_ratio)
    k4 = is.infinite(o_ratio)
    k5 = y_ratio==0
    k6 = o_ratio==0
    y_ratio = y_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    o_ratio = o_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    return(data.frame(y=y_ratio, o=o_ratio))
}

cal_return_ref = function(name_list, dat, exon_num, ifSecBase){
    s = tapply(dat$Length, dat$tx_name, length)
    keep = s>=3
    dat = dat[keep,]
    name_list = anno$normalized_run
    new = lapply(name_list, function(x) dat[,x])
    new = data.frame(new)
    colnames(new) = name_list
    n = nrow(anno[anno$group=='Y',])
    if(ncol(new)>2){
        y = rowMeans(data.frame(new[,1:n]))
        o = rowMeans(data.frame(new[,(n+1):ncol(new)]))
    } else{
        y = new[,1]
        o = new[,2]
    }
    
    new_dat = cbind(dat[,1:6], y=y, o=o)
    y_ratio = tapply(new_dat$y, new_dat$tx_name, getratio, exon_num, ifSecBase)
    o_ratio = tapply(new_dat$o, new_dat$tx_name, getratio, exon_num, ifSecBase)
    
    k1 = is.na(y_ratio)
    k2 = is.na(o_ratio)
    k3 = is.infinite(y_ratio)
    k4 = is.infinite(o_ratio)
    k5 = y_ratio==0
    k6 = o_ratio==0
    y_ratio = y_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    o_ratio = o_ratio[!k1&!k2&!k3&!k4&!k5&!k6]
    out = data.frame(y=y_ratio, o=o_ratio)
    filtered.new_dat = dplyr::filter(new_dat, new_dat$tx_name %in% rownames(out))
    out$ref_exon_y = tapply(filtered.new_dat$y, filtered.new_dat$tx_name, function(x) x[1])
    out$ref_exon_o = tapply(filtered.new_dat$o, filtered.new_dat$tx_name, function(x) x[1])
    return(out)
}

cal_pvalue = function(anno, dat, exon_num){
    s = tapply(dat$Length, dat$tx_name, length)
    keep = s>=3
    dat = dat[keep,]
    name_list = anno$normalized_run
    name_list = gsub('.nc', '',name_list)
    
    new = lapply(name_list, function(x) dat[,x])
    new = data.frame(new)
    colnames(new) = name_list
    n = nrow(anno[anno$group=='Y',])
    if(ncol(new)>2){
        y = rowMeans(new[,1:n])
        o = rowMeans(new[,(n+1):ncol(new)])
    } else{
        y = new[,1]
        o = new[,2]
    }
    
    new_dat = cbind(dat[,1:6], y=y, o=o)
    y_l = tapply(new_dat$y, new_dat$tx_name, getvalue, exon_num)
    o_l = tapply(new_dat$o, new_dat$tx_name, getvalue, exon_num)
    
    exons = unlist(sapply(y_l,'[',1))
    other = unlist(sapply(y_l,'[',2))
    y = data.frame(exon=exons, other=other)

    exons = unlist(sapply(o_l,'[',1))
    other = unlist(sapply(o_l,'[',2))
    o = data.frame(exon=exons, other=other)
    
    k1 <- rowSums(is.na(y))==0
    k2 <- rowSums(is.na(o))==0

    y=y[k1&k2,]
    o=o[k1&k2,]
    
    if(nrow(y)==0 || nrow(o)==0){
        #print(13)
        return(NA)
    }
    
    l=list()
    for(i in 1:nrow(y)){
        s = rbind(y[i,],o[i,])
        l[[gsub('.exon','',row.names(y[i,]))]] = chisq.test(s)$p.value
    }
    
    out = data.frame(unlist(l))
    colnames(out)[1] = 'chisq_pvalue'
    return(out)
}

pieP = function(fc, name){
    png(paste0(name, '_pie.png'), width=5, height = 5, units = 'in', res=1200)
    print(pie_chart(fc))
    dev.off()
    
}

pie_chart = function(fc){
    library(scales)
    fc = log2(fc$o)-log2(fc$y)
    x = data.frame('FC GT 1'=sum(fc>1),
                   'FC LT -1'=sum(fc< -1),
                   'middle'=length(fc)-sum(fc>1)-sum(fc< -1)
    )
    x = gather(x)
    s = sum(x$value)
    ggplot(x, aes(x="", y=value, fill=key)) + geom_bar(width = 1, stat = "identity") + 
        coord_polar("y", start=0) + 
        blank_theme + theme(axis.text.x=element_blank(), legend.position="none")
}


normalize_factor = function(files,id_tab, ignoreTxVersion){
    library(tximport)
    txi = tximport(files, type='salmon', tx2gene=id_tab,ignoreTxVersion=ignoreTxVersion)
    library(edgeR)
    y=DGEList(counts=txi$counts)
    y=calcNormFactors(y)
    factor = (y$samples$lib.size * y$samples$norm.factors)/mean(y$samples$lib.size)
    names(factor) = row.names(y$samples)
    return(factor)
}

get_cpm = function(files,id_tab, ignoreTxVersion, db){
    library(tximport)
    txi = tximport(files, type='salmon', tx2gene=id_tab,ignoreTxVersion=ignoreTxVersion)
    library(edgeR)
    y=DGEList(counts=txi$counts)
    y=calcNormFactors(y)
    factor = (y$samples$lib.size * y$samples$norm.factors)/mean(y$samples$lib.size)
    s=cpm(y)
    s = data.frame(s)
    s$id=rownames(s)
    ids = getGeneAnno(geneID=s$id, annoDb = db)
    colnames(ids)[1]='id'
    if(grepl("\\.",s$id)){
        s = separate(s,'id', c('id', 'na'))
        s = s[,-ncol(s)]
    }
    s = merge(s, ids, by='id')
    return(s)
}

get_abundance = function(files,id_tab, ignoreTxVersion, db){
    library(tximport)
    txi = tximport(files, type='salmon', tx2gene=id_tab,ignoreTxVersion=ignoreTxVersion)
    s = data.frame(txi$abundance)
    s$id=rownames(s)
    ids = getGeneAnno(geneID=s$id, annoDb = db)
    colnames(ids)[1]='id'
    if(grepl("\\.",s$id)){
        s = separate(s,'id', c('id', 'na'))
        s = s[,-ncol(s)]
    }
    s = merge(s, ids, by='id')
    return(s)
}

calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
    inputVector <- sort(inputVector)
    inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
    slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
    xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
    y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
    
    if(drawPlot){  #if TRUE, draw the plot
        plot(1:length(inputVector), inputVector,...)
        b <- y_cutoff-(slope* xPt)
        abline(v= xPt,h= y_cutoff,lty=2,col=8)
        points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
        abline(coef=c(b,slope),col=2)
        title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
        axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    }
    return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector), xPt=xPt, slope=slope))
}

numPts_below_line <- function(myVector,slope,x){
    yPt <- myVector[x]
    b <- yPt-(slope*x)
    xPts <- 1:length(myVector)
    return(sum(myVector<=(xPts*slope+b)))
}

blank_theme <- theme_minimal()+
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
    )


my_theme <- theme_minimal()+
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, face='bold', color='black'),
        axis.ticks.y = element_line(color="black", size = 1),
        axis.ticks.length=unit(.25, "cm"),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        #legend.title = element_blank(),
        legend.position = "none"
    )


theme=theme_classic() + theme(legend.position="none", 
                              element_line(colour = 'black', size = 1),
                              text = element_text(size=30, colour='black'),
                              axis.ticks=element_blank())
