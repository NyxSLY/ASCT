## -------------------------------------------------------------------------------------------------------------
library(GenomicFeatures)
library(tidyr)
library(dplyr)
library(data.table)
source('exon_ratio_of_ratios_functions.R')
source('annotation.R')
args = commandArgs(trailingOnly=TRUE)

#------------------
# load gtf
#------------------
#gtff_mm = '/Volumes/LACIE/mouse_database/Mus_musculus.GRCm38.80.gtf'
#gtff = '/Volumes/LACIE/Human_database/hg19/gencode.v19.annotation.gtf'
#dir='/Volumes/Luyang/work/more_pData/forTL/E-MTAB-4652_02'
#fil='/Volumes/Luyang/work/more_pData/forTL/E-MTAB-4652_02/featureCountExon.txt'
#annof='/Volumes/Luyang/work/more_pData/forTL/E-MTAB-4652_02/anno.txt'
#species='human'

dir=args[1]
fil=args[2]
annof=args[3]
species=args[4]
outdir=args[5]

if(is.na(args[5])){
    out_dir = file.path(dir, 'analysis')
} else{
    out_dir = file.path(dir, outdir)
}

dir.create(out_dir)


if(species=='human'){
    library(EnsDb.Hsapiens.v75)
    hg19=EnsDb.Hsapiens.v75
    tx = transcriptsBy(hg19)
    db = 'EnsDb.Hsapiens.v75'
} else if(species=='mouse'){
    library(EnsDb.Mmusculus.v79)
    mm10 = EnsDb.Mmusculus.v79
    tx = transcriptsBy(mm10)
    db = 'EnsDb.Mmusculus.v79'
} else if(species=='zebrafish'){
    library(org.Dr.eg.db)
    zebrafish = org.Dr.eg.db
    tx = transcriptsBy(zebrafish)
    db = 'org.Dr.eg.db'
} else{
    print('need to check species')
}

tx = data.frame(tx)
tx$seqnames = paste0('chr', tx$seqnames)
id_tab_1 = tx[,c('tx_id', 'group_name')]
id_tab =tx[,c('tx_id', 'group_name', 'width')]

colnames(id_tab_1)[1] = 'tx_name'
colnames(id_tab)[1] = 'tx_name'

anno=read.table(annof, header=T, stringsAsFactors = F, sep='\t')
anno = anno[anno$group %in% c('Y','O'),]
anno$group = factor(anno$group, levels=c('Y','O'))
anno=anno[order(anno$group, anno$run),]

setwd(dir)
dirs = list.dirs(recursive=F)
dirs = dirs[grepl('_salmon', dirs)]
files = sapply(anno$run,function(x) file.path(dir,basename(dirs[grepl(x, dirs)]), 'quant.sf'))
fact = normalize_factor(files, id_tab_1, ignoreTxVersion=T)

write.csv(fact, file.path(out_dir,'tmm_factor.csv'), quote=F)

# # find the main tx

#------------------
# read in salmon result
#------------------

# used for nature SETD2 KD
#dirs = dirs[c(1,2)]

#salmons = lapply(files, read.table, header = T)
salmons = lapply(files, fread, header = T)
x.salmons = data.frame(salmons[[1]][,c(1:3)])

if(grepl("\\.",x.salmons$Name[1])){
    x.salmons = separate(x.salmons,'Name', c('tx_name', 'na'))
    x.salmons = x.salmons[,-2]
} else{
    colnames(x.salmons)[1] = 'tx_name'
}

colnames(x.salmons)[1] = 'tx_name'
for(i in 1:length(salmons)){
    x.salmons = cbind(x.salmons,salmons[[i]]$TPM)
    colnames(x.salmons)[ncol(x.salmons)] = paste0(names(salmons)[i],'_TPM')
}

# calculate average TPM
x.salmons$TPM = rowSums(x.salmons[,c(4:ncol(x.salmons))])/(ncol(x.salmons)-3)
#head(x.salmons)
xx = merge(x.salmons, id_tab, by='tx_name')

#head(xx)
xx = arrange(xx, group_name, desc(TPM))
#head(xx)
xx.main = xx[!duplicated(xx$group_name),]

# filter low
xx.main.fil = dplyr::filter(xx.main, TPM>1)
#dim(xx.main.fil)
fil.list = xx.main.fil[,c('tx_name', 'group_name', 'TPM')]
#head(fil.list)
#dim(fil.list)
#quantile(xx.main$TPM)

# -------------------------------------------------
# -------------------------------------------------
# get major tx coord and anno
# tx > 3kb
major_tx = merge(fil.list, tx, by='tx_name')
major_tx = dplyr::filter(major_tx, width>3000)
major_tx_bed = major_tx[,c('seqnames', 'start', 'end', 'tx_name', 'TPM', 'strand')]
write.table(major_tx_bed, file.path(out_dir, 'major_transcript.bed'), sep='\t', quote = F, col.names = F,
            row.names = F)

# get all tx coord - sorted by TPM
# tx > 3kb
# TPM > 1, about 50%
all_tx = merge(xx, tx, by='tx_name')
all_tx = dplyr::filter(all_tx, width.x>3000, TPM>1)
all_tx_bed = all_tx[,c('seqnames', 'start', 'end', 'tx_name', 'TPM', 'strand')]
write.table(all_tx_bed, file.path(out_dir, 'all_transcript.bed'), sep='\t', quote = F, col.names = F,
            row.names = F)

# get longest tx of each gene coord - sorted by TPM
# tx > 3kb
# TPM > 1, about 50%
all_tx = merge(xx, tx, by='tx_name')
all_tx = arrange(all_tx, group_name.x, desc(Length))
long_tx = all_tx[!duplicated(all_tx$group_name.x),]
long_tx = dplyr::filter(long_tx, width.x>3000, TPM>1)
long_tx_bed = long_tx[,c('seqnames', 'start', 'end', 'tx_name', 'TPM', 'strand')]
write.table(long_tx_bed, file.path(out_dir, 'long_transcript.bed'), sep='\t', quote = F, col.names = F,
            row.names = F)


#------------------
# read in feature count exon result
#------------------

read_count = function(fil){
  if(grepl(':', fil)){
    nfil=unlist(strsplit(fil,':'))
  } else{
    nfil = fil
  }

  fil0 = nfil[1]
  con<-file(fil0)
  open(con)
  rn = read.table(con,skip=1,nrow=1, stringsAsFactors = F)
  close(con)

  x = read.table(fil0, skip=2, stringsAsFactors = F)
  colnames(x) = as.character(t(rn))
  colnames(x)[1] = 'tx_name'

  if(length(nfil) > 1){
    for(i in 2:length(nfil)){
      fil1 = nfil[i]
      con<-file(fil1)
      open(con)
      rn1 = read.table(con,skip=1,nrow=1, stringsAsFactors = F)
      close(con)
      rn1v = as.character(rn1)
      rn1v = rn1v[7:length(rn1v)]

      # remove original data first
      x = x[, !(names(x) %in% rn1v)]

      # add bu data
      x1 = read.table(fil1, skip=2, stringsAsFactors = F)
      cn = colnames(x)
      cn = c(cn, rn1v)
      x = cbind(x, x1[,7:ncol(x1)])
      colnames(x) = cn
    }
  }
  return(x)
}

x = read_count(fil)

if(grepl("\\.",x$tx_name[1])){
    x = separate(x,'tx_name', c('tx_name', 'na'))
    x = x[,-2]
}

# change column order to anno file
new = x[,c(1:6)]
for(i in 1:nrow(anno)){
    new = cbind(new,x[,grepl(anno$run[i],colnames(x))])
    colnames(new)[i+6] = anno$run[i]
}
x = new
#head(x)
# ----------------
# apply TMM normalize factor.
# ----------------

sa = data.frame(x[,7:ncol(x)])
x.sa = data.frame(x[,1:6])

for(i in 1:ncol(sa)){
    print(i)
    x.sa = cbind(x.sa,data.frame(sa[,i] / fact[i]))
}
colnames(x.sa)[7:ncol(x.sa)] = colnames(x)[7:ncol(x)]
x = x.sa
#head(x)

name_list = c()
for(i in 7:ncol(x)){
    x = cbind(x,x[,i]/x$Length)
    colnames(x)[ncol(x)]=paste0(colnames(x)[i],'.nc')
    name_list = c(name_list, paste0(colnames(x)[i],'.nc'))
}

#head(x)
#name_list

dat = merge(x, fil.list, by='tx_name')
#head(dat)

dat.fwd = dplyr::filter(dat, Strand=='+')
dat.rev = dplyr::filter(dat, Strand=='-')
dat.fwd = arrange(dat.fwd, group_name, Start)
dat.rev = arrange(dat.rev, group_name, desc(Start))
dat = rbind(dat.fwd, dat.rev)
dat$tx_name = as.character(dat$tx_name)
dat$Chr = as.character(dat$Chr)
setwd(out_dir)
write.csv(dat, 'filtered data used for calculating exon ratios.csv')
#------------------
# functions
#------------------



# get sample names
#name_list
#deparse(substitute(name_list[1]))



## -------------------------------------------------------------------------------------------------------------
setwd(out_dir)
anno$normalized_run = paste0(anno$run, '.nc')

r2 = cal(anno, dat, 2, F)
r3 = cal(anno, dat, 3, F)
r4 = cal(anno, dat, 4, F)
r5 = cal(anno, dat, 5, F)
rl = cal(anno, dat, -1, F)
rl2 = cal(anno, dat, -2, F)
rl3 = cal(anno, dat, -3, F)
other = getratio_other_exons(dat, anno, F)

r2x = cal(anno, dat, 2, T)
r3x = cal(anno, dat, 3, T)
r4x = cal(anno, dat, 4, T)
r5x = cal(anno, dat, 5, T)
rlx = cal(anno, dat, -1, T)
otherx = getratio_other_exons(dat, anno, T)

ggplot_boxplot(r2, 'second_exon')
ggplot_boxplot(r3, 'third_exon')
ggplot_boxplot(r4, '4th_exon')
ggplot_boxplot(rl, 'last_exon')
ggplot_boxplot(other, 'all_other_exon')

#------------------------------------
# vs the second exon
#------------------------------------


ggplot_boxplot(r3x, 'secondVs.third_exon')
ggplot_boxplot(r4x, 'secondVs.4th_exon')
ggplot_boxplot(rlx, 'secondVs.last_exon')
ggplot_boxplot(otherx, 'secondVs.all_other_exon')

# pvalue
r2.p=wilcox.test(r2[,1], r2[,2], paired=T)
r3.p=wilcox.test(r3[,1], r3[,2], paired=T)
r4.p=wilcox.test(r4[,1], r4[,2], paired=T)
rl.p=wilcox.test(rl[,1], rl[,2], paired=T)

r3x.p=wilcox.test(r3x[,1], r3x[,2], paired=T)
r4x.p=wilcox.test(r4x[,1], r4x[,2], paired=T)
rlx.p=wilcox.test(rlx[,1], rlx[,2], paired=T)

dat.pvalue = data.frame(pvalue=c(r2.p$p.value, r3.p$p.value, r4.p$p.value, rl.p$p.value,
                                 r3x.p$p.value, r4x.p$p.value, rlx.p$p.value))
dat.pvalue$name = c('2ndVsFirst', '3rdVsFirst', '4thVsFirst', 'LastVsFirst',
                    '3rdVsSecond', '4thVsSecond', 'LastVsSecond')

write.csv(dat.pvalue, 'pvalues.csv')



plot_dat = rbind(gather(data.frame(second=log2(r2$o)-log2(r2$y))),
                 gather(data.frame(third=log2(r3$o)-log2(r3$y))),
                 gather(data.frame(fourth=log2(r4$o)-log2(r4$y))),
                 gather(data.frame(last=log2(rl$o)-log2(rl$y))))
plot_dat$type='first'

plot_dat1 = rbind(
                  gather(data.frame(thirdx=log2(r3x$o)-log2(r3x$y))),
                  gather(data.frame(fourthx=log2(r4x$o)-log2(r4x$y))),
                  gather(data.frame(lastx=log2(rlx$o)-log2(rlx$y))))
plot_dat1$type='second'

plot_dat = rbind(plot_dat, plot_dat1)


# pvalue
r2.p1=wilcox.test(log2(r2$o)-log2(r2$y))
r3.p1=wilcox.test(log2(r3$o)-log2(r3$y))
r4.p1=wilcox.test(log2(r4$o)-log2(r4$y))
rl.p1=wilcox.test(log2(rl$o)-log2(rl$y))

r3x.p1=wilcox.test(log2(r3x$o)-log2(r3x$y))
r4x.p1=wilcox.test(log2(r4x$o)-log2(r4x$y))
rlx.p1=wilcox.test(log2(rlx$o)-log2(rlx$y))

dat.pvalue = data.frame(pvalue=c(r2.p1$p.value, r3.p1$p.value, r4.p1$p.value, rl.p1$p.value,
                                 r3x.p1$p.value, r4x.p1$p.value, rlx.p1$p.value))
dat.pvalue$name = c('2ndVsFirst', '3rdVsFirst', '4thVsFirst', 'LastVsFirst',
                    '3rdVsSecond', '4thVsSecond', 'LastVsSecond')
dat.pvalue$average = c(mean(log2(r2$o)-log2(r2$y)),
                       mean(log2(r3$o)-log2(r3$y)),
                       mean(log2(r4$o)-log2(r4$y)),
                       mean(log2(rl$o)-log2(rl$y)),
                       mean(log2(r3x$o)-log2(r3x$y)),
                       mean(log2(r4x$o)-log2(r4x$y)),
                       mean(log2(rlx$o)-log2(rlx$y)))

write.csv(dat.pvalue, 'pvalues of the main boxplot.csv')

# output data for main boxplot

max.len = max(length(r2$o), length(r2$y), length(r3$o), length(r3$y),length(r4$o), length(r4$y),length(rl$o), length(rl$y))


outdatahere = data.frame(ratio_exon2_y = c(r2$y, rep(NA, max.len - length(r2$y))),
ratio_exon2_o = c(r2$o, rep(NA, max.len - length(r2$o))),
ratio_exon3_y=c(r3$y, rep(NA, max.len - length(r3$y))),
ratio_exon3_o=c(r3$o, rep(NA, max.len - length(r3$o))),
ratio_exon4_y=c(r4$y, rep(NA, max.len - length(r4$y))),
ratio_exon4_o=c(r4$o, rep(NA, max.len - length(r4$o))),
ratio_exonL_y=c(rl$y, rep(NA, max.len - length(rl$y))),
ratio_exonL_o=c(rl$o, rep(NA, max.len - length(rl$o))))

write.csv(outdatahere, 'data used in the main boxplot - vs 1stExon.csv')

max.len = max(length(r2x$o), length(r2x$y), length(r3x$o), length(r3x$y),length(r4x$o), length(r4x$y),length(rlx$o), length(rlx$y))

outdatahere2 = data.frame(ratio_exon2_y = c(r2x$y, rep(NA, max.len - length(r2x$y))),
  ratio_exon2_o = c(r2x$o, rep(NA, max.len - length(r2x$o))),
ratio_exon3_y=c(r3x$y, rep(NA, max.len - length(r3x$y))),
ratio_exon3_o=c(r3x$o, rep(NA, max.len - length(r3x$o))),
ratio_exon4_y=c(r4x$y, rep(NA, max.len - length(r4x$y))),
ratio_exon4_o=c(r4x$o, rep(NA, max.len - length(r4x$o))),
ratio_exonL_y=c(rlx$y, rep(NA, max.len - length(rlx$y))),
ratio_exonL_o=c(rlx$o, rep(NA, max.len - length(rlx$o))))

write.csv(outdatahere2 , 'data used in the main boxplot - vs 2ndExon.csv')



plot_dat$key = factor(plot_dat$key, levels=c('second', 'third', 'fourth', 'last','thirdx','fourthx','lastx'))
plot_dat = boxplot_filter(plot_dat)
name='ratio'
#fill="#F8AC02"
p = ggplot(plot_dat, aes(x=key, y=value, fill=type)) +
    geom_boxplot(outlier.shape=NA, size=1) +
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
    xlab('') +
    # scale_y_continuous(breaks=seq(-1.5,1.5,0.5)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    scale_fill_manual(values=c(alpha("#F8AC02",0.6), "#FCF9CE"))

pdf(paste0(name, '_boxplot.pdf'), width=10, height = 5.5)
print(p)
dev.off()

#------------------------------------
# quartile TPM
#------------------------------------
qtl=quantile(dat$TPM)
q1 = dplyr::filter(fil.list, TPM<qtl[2])
q2 = dplyr::filter(fil.list, TPM<qtl[3], TPM>=qtl[2])
q3 = dplyr::filter(fil.list, TPM<qtl[4], TPM>=qtl[3])
q4 = dplyr::filter(fil.list, TPM>=qtl[4])

dat.q1 = dplyr::filter(dat, tx_name %in% q1$tx_name)
dat.q2 = dplyr::filter(dat, tx_name %in% q2$tx_name)
dat.q3 = dplyr::filter(dat, tx_name %in% q3$tx_name)
dat.q4 = dplyr::filter(dat, tx_name %in% q4$tx_name)

make_boxplot_4quantiles = function(dat){
    r2 = cal(anno, dat, 2, F)
    r3 = cal(anno, dat, 3, F)
    r4 = cal(anno, dat, 4, F)
    r5 = cal(anno, dat, 5, F)
    rl = cal(anno, dat, -1, F)
    rl2 = cal(anno, dat, -2, F)
    rl3 = cal(anno, dat, -3, F)

    plot_dat = rbind(gather(data.frame(second=log2(r2$o)-log2(r2$y))),
                     gather(data.frame(third=log2(r3$o)-log2(r3$y))),
                     gather(data.frame(fourth=log2(r4$o)-log2(r4$y))),
                     gather(data.frame(last=log2(rl$o)-log2(rl$y))))
    plot_dat$type='first'

    plot_dat$key = factor(plot_dat$key, levels=c('second', 'third', 'fourth', 'last','secondx','thirdx','fourthx','lastx'))
    plot_dat = boxplot_filter(plot_dat)
    maxx = max(plot_dat$value)*1.25
    minn = min(plot_dat$value)*1.25
    name='ratio'
    #fill="#F8AC02"
    p = ggplot(plot_dat, aes(x=key, y=value, fill=type)) +
        geom_boxplot(outlier.shape=NA, size=1) +
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
        xlab('') +
        # scale_y_continuous(breaks=seq(-1.5,1.5,0.5)) +
        geom_hline(yintercept=0, linetype="dashed", color = "black") +
        scale_fill_manual(values=c(alpha("#F8AC02",0.6), "#FCF9CE"))

    # pvalue
    r2.p1=wilcox.test(log2(r2$o)-log2(r2$y))
    r3.p1=wilcox.test(log2(r3$o)-log2(r3$y))
    r4.p1=wilcox.test(log2(r4$o)-log2(r4$y))
    rl.p1=wilcox.test(log2(rl$o)-log2(rl$y))
    dat.pvalue = data.frame(pvalue=c(r2.p1$p.value, r3.p1$p.value, r4.p1$p.value, rl.p1$p.value))
    dat.pvalue$name = c('2ndVsFirst', '3rdVsFirst', '4thVsFirst', 'LastVsFirst')
    dat.pvalue$average = c(mean(log2(r2$o)-log2(r2$y)),
                       mean(log2(r3$o)-log2(r3$y)),
                       mean(log2(r4$o)-log2(r4$y)),
                       mean(log2(rl$o)-log2(rl$y)))
    return(list(p=p, max=maxx, min=minn, pvalue=dat.pvalue))
}

dat.q1.p = make_boxplot_4quantiles(dat.q1)
dat.q2.p = make_boxplot_4quantiles(dat.q2)
dat.q3.p = make_boxplot_4quantiles(dat.q3)
dat.q4.p = make_boxplot_4quantiles(dat.q4)

maxx = max(dat.q1.p[['max']],dat.q2.p[['max']],dat.q3.p[['max']],dat.q4.p[['max']])
minn = min(dat.q1.p[['min']],dat.q2.p[['min']],dat.q3.p[['min']],dat.q4.p[['min']])

# write p-value
write.csv(dat.q1.p[['pvalue']], 'pvalues of Q1 boxplot.csv')
write.csv(dat.q2.p[['pvalue']], 'pvalues of Q2 boxplot.csv')
write.csv(dat.q3.p[['pvalue']], 'pvalues of Q3 boxplot.csv')
write.csv(dat.q4.p[['pvalue']], 'pvalues of Q4 boxplot.csv')


pdf('ratio boxplot of 4quantiles.pdf', width=10, height=11)
cowplot::plot_grid(
    dat.q1.p[['p']],
    dat.q2.p[['p']],
    dat.q3.p[['p']],
    dat.q4.p[['p']],
    nrow=2
)
dev.off()
#------------------------------------
# pie chart
#------------------------------------

library(cowplot)
pdf('piechart.pdf')
cowplot::plot_grid(
    pie_chart(r2),
    pie_chart(r3),
    pie_chart(r4),
    pie_chart(rl),
    pie_chart(other),
    nrow=1,
    labels = c('exon2','exon3','exon4','exonl','all_other')
)
dev.off()

fclog = function(fc, name){
    fc = data.frame(fc=log2(fc$o)-log2(fc$y))
    x = data.frame('More'=sum(fc>1),'Less'=sum(fc< -1))
    x = gather(x)
    x$name = name
    return(x)
}
#------------------------------------
# barplot
#------------------------------------
dat.plot = rbind(fclog(r2, 'r2'),
                 fclog(r3, 'r3'),
                 fclog(r4, 'r4'),
                 fclog(rl, 'rl'))
pdf('barplot.pdf', height=5, width=7)
ggplot(dat.plot, aes(x=name, y=value, fill=key)) + geom_bar(color='black', stat = "identity",
                                                            position='dodge',width=0.8) +
    theme_classic() +
    theme(axis.text.y=element_text(size=20, colour = 'black'),
          axis.title.y = element_text(size=20, colour = 'black'),
          axis.line = element_blank(),
          axis.ticks.length.y = unit(0.25, "cm"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position="right") +
    labs(y='number of genes', fill='Group')+
    scale_fill_manual(values = c('#18BDC2','#F3766E'))
dev.off()


#------------------------------------
# quantile TPM - 10 groups
#------------------------------------

qtl=quantile(dat$TPM, seq(0, 1,0.1))
l = list()
for(i in 1:(length(qtl)-1)){
    a = dplyr::filter(fil.list, TPM<qtl[i+1], TPM>=qtl[i])
    l[[i]] = dplyr::filter(dat, tx_name %in% a$tx_name)
}
new =lapply(l, getratio_other_exons, anno=anno, ifSecBase=F)
for(i in 1:length(new)){
    new[[i]]$type=i
    new[[i]]$logratio = log2(new[[i]]$o/new[[i]]$y)
}
new = do.call(rbind, new)
colnames(new)[4] = 'value'
new = boxplot_filter(new)
p = ggplot(new, aes(x=as.factor(type), y=value, group=)) +
    geom_boxplot(outlier.shape=NA, size=1) +
    stat_summary(fun.y=mean, colour="darkred", geom="point",
                 shape=18, size=3,show.legend = FALSE)+
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
    #xlab('') +
    #scale_y_continuous(breaks=seq(-0.2,0.6,0.2)) +
    #ylim(-2,2) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    scale_fill_manual(values=c(alpha("#F8AC02",0.6), "#FCF9CE"))


pdf('10group.pdf', width=10, height = 5.5)
print(p)
dev.off()

# getratio_other_exons(other, 'secondVs.all_other_exon', F)

# generate CPM table

cpm.tab = get_cpm(files, id_tab, T, db)

# get abundance
abs.tab = get_abundance(files, id_tab, T, db)

write.csv(abs.tab, 'abundance_tab.csv', quote = F)


## -------------------------------------------------------------------------------------------------------------
l = list()
n = 2  # don't change this
nn = n  # store original n value

while(1){
    res = cal(anno, dat, n, F)
    m = nrow(res)
    if(m>0){
        x = log2(res$o) - log2(res$y)
        print(n)
        names(x) = rownames(res)
        l[[as.character(n)]] = x
        n = n + 1
    } else{
        # get ref exon data
        out = cal_return_ref(anno, dat, nn, F)
        ref_exon_data = out[,c(3,4)]
        ref_exon_data$tx_name = rownames(ref_exon_data)
        break
    }
}

tx_name = dat$tx_name%>%unique
base = data.frame(tx_name)
row.names(base) = tx_name
base = merge(base, ref_exon_data, by='tx_name', all=T)

for(i in 1:length(l)){
    print(i)
    input = data.frame(l[[i]])
    colnames(input) = paste0('exon',i+1)
    input$tx_name = row.names(input)
    base = merge(base, input, by='tx_name', all.x=T)
}

# make a table
temp = xx.main.fil
colnames(temp)[1] = 'tx_name'
tab = merge(temp, base, by='tx_name')

gene_anno = getGeneAnno(tab$group_name, annoDb = db)
colnames(gene_anno)[1] = 'group_name'
if(grepl("\\.", tab$group_name)[1]){
    tab = separate(tab,'group_name', c('group_name', 'na'))
    tab = tab[,-(which(colnames(tab)=='na'))]
}

write.csv(tab, 'exon ratio of ratios.csv',quote = F)

#----------------------------------------
# try with second exon

l = list()
n = 2
nn = n  # store original n value
l_o = list()
l_y = list()

while(1){
    res = cal(anno, dat, n, T) # use the second exon as base
    m = nrow(res)
    if(m>0){
        x = log2(res$o) - log2(res$y)
        y = log2(res$y)
        o = log2(res$o)

        print(n)
        names(x) = rownames(res)
        names(y) = rownames(res)
        names(o) = rownames(res)
        l[[as.character(n)]] = x
        l_o[[as.character(n)]] = o
        l_y[[as.character(n)]] = y
        n = n + 1
    } else{
        # get ref exon data
        out = cal_return_ref(anno, dat, nn, T)
        ref_exon_data = out[,c(3,4)]
        ref_exon_data$tx_name = rownames(ref_exon_data)
        break
    }
}

tx_name = dat$tx_name%>%unique
base = data.frame(tx_name)
row.names(base) = tx_name
base = merge(base, ref_exon_data, by='tx_name', all=T)
base_y = base
base_o = base

for(i in 1:length(l)){
    print(i)
    input = data.frame(l[[i]])
    input_y = data.frame(l_y[[i]])
    input_o = data.frame(l_o[[i]])

    colnames(input) = paste0('exon',i+1)
    colnames(input_y) = paste0('exon',i+1)
    colnames(input_o) = paste0('exon',i+1)

    input$tx_name = row.names(input)
    input_y$tx_name = row.names(input_y)
    input_o$tx_name = row.names(input_o)

    base = merge(base, input, by='tx_name', all.x=T)
    base_y = merge(base_y, input_y, by='tx_name', all.x=T)
    base_o = merge(base_o, input_o, by='tx_name', all.x=T)
}

# make a table
temp = xx.main.fil
colnames(temp)[1] = 'tx_name'
tab = merge(temp, base, by='tx_name')
tab_y = merge(temp, base_y, by='tx_name')
tab_o = merge(temp, base_o, by='tx_name')


anno = getGeneAnno(tab$group_name, annoDb = db)
colnames(anno)[1] = 'group_name'
if(grepl("\\.", tab$group_name)[1]){
    tab = separate(tab,'group_name', c('group_name', 'na'))
    tab = tab[,-(which(colnames(tab)=='na'))]
    tab_y = separate(tab_y,'group_name', c('group_name', 'na'))
    tab_y = tab_y[,-(which(colnames(tab_y)=='na'))]
    tab_o = separate(tab_o,'group_name', c('group_name', 'na'))
    tab_o = tab_o[,-(which(colnames(tab_o)=='na'))]
}


tab = merge(gene_anno, tab, by='group_name')
tab_y = merge(gene_anno, tab_y, by='group_name')
tab_o = merge(gene_anno, tab_o, by='group_name')
write.csv(tab, 'exon ratio of ratios_VsSecondExon.csv',quote = F)
write.csv(tab_y, 'exon ratio_VsSecondExon_Y.csv',quote = F)
write.csv(tab_o, 'exon ratio_VsSecondExon_O.csv',quote = F)


## -------------------------------------------------------------------------------------------------------------
setwd(out_dir)
#x = read.csv('exon ratio of ratios.csv')
x = read.csv('exon ratio of ratios_VsSecondExon.csv', stringsAsFactors = F) # get vs second exon data
# for test
#x = read.csv('/Volumes/LACIE/HSC/readcount/exon_ratios/exon ratio of ratios_VsSecondExon.csv')
#x = read.csv('/Volumes/LACIE/MSC_RNA_seq/ucMSC/salmon/quants/exon_ratios/exon ratio of ratios.csv')
#xx = read.csv('exon chisq test.csv')

# TODO: try to find ones that O > Y in the ref exon
first_exon = names(x)[grepl('exon\\d+',names(x))][1]
ind = which(names(x)==first_exon)
h = x[,c(2:6)]
h$TPM = x$TPM
h$ref_exon_y = x$ref_exon_y
h$ref_exon_o = x$ref_exon_o

# -----------
# get maximum ROR of each gene and define a threshold
# -----------
dat1 = x[,c((ind+1):ncol(x))]
#dat1 = x[,c((ind+1):ncol(x))]  ## use ind+1, because I am using second exon data, need to start from exon3
lastExonNa = function(test){
    test[!is.na(test)][length(test[!is.na(test)])] = NA
    test
}
# remove the ratio of last exon, because usually last exon is larger, and CT is
# unusual to happen at last exon

dat1 = apply(dat1, 1, lastExonNa)
dat1 = data.frame(t(dat1))

h$max_ratio = apply(dat1, 1, max, na.rm=T)
h = h[h$max_ratio!=-Inf,]

plot(sort(h$max_ratio), type='l')  # I can make a heatmap here - not here, in all others
abline(h=0)
sort(h$max_ratio[h$max_ratio>0]) %>% plot
out = calculate_cutoff(h$max_ratio[h$max_ratio>=0])
out.neg = calculate_cutoff(abs(h$max_ratio[h$max_ratio<=0]))  # all 0 because I add exon2 vs exon2, it's 0, so in a all negtive gene, max is 0
all = sort(h$max_ratio)
pos = sort(h$max_ratio[h$max_ratio>=0])
neg = sort(abs(h$max_ratio[h$max_ratio<=0]))

# -----------
# plot
# -----------
pdf('rank ror and define threshold.pdf', h=9.6, w=8)
par(mfrow=c(1,1))
plot(all,pch=16,cex=0.5,col='#27AAE1', ylim=c(-6,6))
#points(all,type='l')
#abline(h=0, lty=2)
#plot(all, type='l')
x1 = which(all==pos[out$xPt])
y1 = all[x1]
points(x1,y1,pch=16,cex=0.9,col=2)
x2 = which(all==-neg[out.neg$xPt])[1]
y2 = all[x2][1]
points(x2,y2,pch=16,cex=0.9,col=2)
abline(v= x1,h= y1,lty=2,col=8, lw=1)
abline(v= x2,h= y2,lty=2,col=8, lw=1)
b <- out$absolute-(out$slope* x1)
abline(coef=c(b,out$slope),col=2)
b <- -(out.neg$absolute)-(out.neg$slope * x2)
abline(coef=c(b,out.neg$slope),col=2)
abline(h=0, lty=2, lw=1.5)
par(mfrow=c(1,1))
dev.off()

# -----------
# start to process ct -
# -----------

ct = h[h$max_ratio>out$absolute,]

which_exon = function(dat, x){
    out=list()
    for(i in 1:nrow(dat)){
        r = dat[i,2]
        id = dat[i,]$tx_name
        inx = which(x[x$tx_name==id,]==r)[1]
        out[[i]] = names(x[x$tx_name==id,][inx])
    }
    return(out)
}

exon_name = which_exon(ct[,c('tx_name', "max_ratio")],x)
ct$which_exon = unlist(exon_name)
exon_name_all = which_exon(h[,c('tx_name', "max_ratio")],x)
h$which_exon = unlist(exon_name_all)
get_coord = function(ct, dat){
    res = list()
    for(i in 1:nrow(ct)){
        new = dat[dat$tx_name==ct[i,]$tx_name,]
        exon_num = as.numeric(gsub('exon','',ct[i,]$which_exon))
        coord = new[exon_num,][c(2,3,4)]
        coord = paste(as.character(coord[1]), coord[2], coord[3], sep='_')
        res[[i]] = coord
    }
    return(res)
}

coords = get_coord(ct, dat)
ct$region = unlist(coords)
coords_all = get_coord(h, dat)
h$region = unlist(coords_all)
# first filter. Only count following exons' ROR all >= 0
l = list()
follow = list()  # following exon, how many has log2 RoR < 0
for(i in 1:nrow(ct)){
    exon_num = as.numeric(gsub('exon','',ct[i,]$which_exon))
    next_exon = paste0('exon',exon_num+1)
    ind = which(names(x)==next_exon)
    left = x[x$tx_name==ct[i,]$tx_name,][ind:ncol(x)]
    left = left[!is.na(left)]
    l[[i]] = x[x$tx_name==ct[i,]$tx_name,][,next_exon]
    follow[[i]] = sum(left<0)
}

ct$next_exon = sapply(l, '[', 1)
ct$follow_exons = sapply(follow, '[', 1)
colnames(ct)[2] = 'gene_name'
ct.f = dplyr::filter(ct, follow_exons==0)
ct.f = ct.f[!ct.f$region=='NA_NA_NA',]

# filter ref exon ratio between Y and O
ref_exon_ratio = log2(h$ref_exon_o/h$ref_exon_y)
ref_exon_ratio = ref_exon_ratio[!is.na(ref_exon_ratio)]
ref_exon_ratio = ref_exon_ratio[!is.infinite(ref_exon_ratio)]
plot(sort(ref_exon_ratio))
abline(h=0)
ref.out = calculate_cutoff(ref_exon_ratio[ref_exon_ratio>0])
ref.out.neg = calculate_cutoff(abs(ref_exon_ratio[ref_exon_ratio<0]))
interval = c(-ref.out.neg$absolute,ref.out$absolute)
# end

ct.ff = dplyr::filter(ct.f, log2(ref_exon_o/ref_exon_y) > interval[1])
#log2(ref_exon_o/ref_exon_y) < interval[2])

write.csv(ct.f, 'ct_gene_listVsSecondExon.csv',quote=F)
write.csv(ct.ff, 'ct_gene_listVsSecondExon.ff.csv',quote=F)

tx_bed = tx[,c('tx_name', 'seqnames', 'start', 'end', 'strand')]

ct.bed = merge(tx_bed, ct.f, by='tx_name')
ct.bed = ct.bed[,c(2,3,4,1,12,5)]
write.table(ct.bed,'ct_gene_listVsSecondExon.bed', quote = F, sep='\t', row.names = F, col.names = F)

ct.ff.bed = merge(tx_bed, ct.ff, by='tx_name')
ct.ff.bed = ct.ff.bed[,c(2,3,4,1,12,5)]
write.table(ct.ff.bed,'ct_gene_listVsSecondExon.ff.bed', quote = F, sep='\t', row.names = F, col.names = F)


## -------------------------------------------------------------------------------------------------------------
ind = which(names(ct.ff[i,])=='which_exon')
get_seperate_ratio = function(x, dat){
    exon = x['which_exon']
    tx = x['tx_name']
    n = as.integer(gsub('exon','', exon))
    ss = cal(anno, dat, n, T)
    out = ss[grepl(tx, rownames(ss)),]
    return(out)
}

y_ratio = read.csv('exon ratio_VsSecondExon_Y.csv', stringsAsFactors = F)
o_ratio = read.csv('exon ratio_VsSecondExon_O.csv', stringsAsFactors = F)

ct.ff = read.csv('ct_gene_listVsSecondExon.ff.csv', stringsAsFactors = F)

get_seperate_ratio1 = function(x, dat_y, dat_o){
    exon = x['which_exon']
    tx = x['tx_name']
    y = dat_y[dat_y$tx_name==tx,]
    o = dat_o[dat_o$tx_name==tx,]
    out = data.frame('y'=y[,exon],'o'=o[,exon], 'tx'=tx)
    rownames(out) = tx
    return(out)
}

out = apply(ct.ff, 1, get_seperate_ratio1, y_ratio, o_ratio)
out = do.call(rbind.data.frame, out)

out.ct = apply(ct, 1, get_seperate_ratio1, y_ratio, o_ratio)
out.ct = do.call(rbind.data.frame, out.ct)

#out = apply(ct.ff, 1, get_seperate_ratio, dat)
#out = do.call(rbind.data.frame, out)
h = h[!h$region=='NA_NA_NA',]
out_all = apply(h, 1, get_seperate_ratio1, y_ratio, o_ratio)
out_all = do.call(rbind.data.frame, out_all)
#hist((out$o/out$y), xlim=c(-10,10))

#mat = log2(out)
mat = out
colnames(mat) = c('Y', 'O')
mat$chr='chr1'
mat$start='1'
mat$end='2'
mat$id=rownames(mat)
mat$sc=0
mat$strand='+'
mat=mat[,c('chr', 'start', 'end', 'id', 'sc', 'strand', 'Y', 'O')]
mat.ff = dplyr::filter(mat, id %in% ct.ff.bed$tx_name)
print(head(mat.ff))
#hist(mat.ff$o-mat.ff$y)

make_dp_mat <- function(mat, out, ratio){
    write.table(mat, 'intermediate_exon_vs_first_exon.res',sep='\t', quote=F,
                row.names = F, col.names = F)

    if(ratio==F){
        headline=paste0('@{"upstream":[0,0],"downstream":[0,0],"body":[100,100],"bin size":[100,100],"ref point":[null,null],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":true,"nan after end":false,"proc number":8,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0,0],"unscaled 3 prime":[0,0],"group_labels":["genes"],"group_boundaries":[0,',
                        as.character(nrow(mat)),
                        '],"sample_labels":["Y","O"],"sample_boundaries":[0,1,2]}')
    }
    else{
        headline=paste0('@{"upstream":[0],"downstream":[0],"body":[100],"bin size":[100],"ref point":[null],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":true,"nan after end":false,"proc number":8,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0],"unscaled 3 prime":[0],"group_labels":["genes"],"group_boundaries":[0,',
                        as.character(nrow(mat)),
                        '],"sample_labels":["ratio"],"sample_boundaries":[0,1]}')
    }


    write.table(headline, 'headline.txt', quote = F, row.names = F, col.names = F)
    system(paste0('cat headline.txt intermediate_exon_vs_first_exon.res > ',out))
    system(paste0('gzip -f ', out))
    system('rm intermediate_exon_vs_first_exon.res headline.txt')
}

make_dp_mat(mat.ff, 'selected.ff.mat',F)
make_dp_mat(mat, 'intermediate_exon_vs_first_exon.mat',F)

left = out_all[!out_all$tx %in% out.ct$tx,]

pdf(paste0(name, '_smoothScatter.pdf'), width=7, height = 7)
smoothScatter(data.frame(left), xlim=c(-10,10), ylim=c(-10,10),
              colramp = colorRampPalette(c("white", "grey4")),
              xlab='ratio in Y', ylab='ratio in O')
points(data.frame(out.ct), pch=20, col='blue', cex=1)
lines(x=c(-10,10),y=c(-10,10), lty=2)
dev.off()

system('/Volumes/luyang-1/work/Figures/Fig1/heatmap.sh selected.ff.mat.gz')
