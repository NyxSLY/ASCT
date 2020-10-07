library(csaw)
library(edgeR)
library(GenomicFeatures)
library(rtracklayer)
source('count.R')
source('commonAPI.R')


#setwd('C:\\Users\\Luyang\\Desktop\\new\\nn\\THOR_method')
#setwd('/Volumes/Data1/ChIP_seq_11_2018/csaw/THOR_method/counts')

find_diffpeak <- function(file1, file2, ses_factor, tmm_sample_info, sample_name,
                          peak_region, reci=F,debug=F){
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
    new_sub[new_sub<0] = 0
    group = sample_names
    sample_info = read.csv(tmm_sample_info, sep='\t')
    y = DGEList(new_sub, lib.size=sample_info$lib.size, norm.factors=sample_info$norm.factors,
                group=group)  # need to change if want to make function

    keep1 = rowSums(cpm(y)>0)>length(sample_name)/2
    y = y[keep1,]

    region.c = region[keep1,]
    # only used for 2 samples. Will use poisson expression
    cpm.filter = cpm(y)

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
    # if k36, it is wide pick, merge is not appropiate
    if(grepl('36', sample_name[1])){
        pre_peak = calculate_nomerge(dat.c, reci=F)
    } else{
        pre_peak = calculate_merge(dat.c, reci=F)
    }

    #--debug--
    if(debug==T){
        write.csv(pre_peak, paste0(basename(file1),'_allTestResult.csv'))
    }

    diffpeak = dplyr::filter(pre_peak, FDR<0.05, fc>1) # may change this
    return(diffpeak)
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

toGR <- function(data){
    dat = makeGRangesFromDataFrame(data, seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3',
                                   keep.extra.columns = T)
    return(dat)
}

# hg19 annotation
if(!exists('gene') | !exists('pro')){
    gtff = '/Volumes/LACIE/Human_database/hg19/gencode.v19.annotation.gtf'
    gtf = makeTxDbFromGFF(gtff, format='gtf')
    pro=promoters(gtf,upstream=2000, downstream=2000)
    gene=genes(gtf)
}


# ----------input zone--------
# this part need to change based on input
label=read_excel('/sequencing\ libraries\ info.xlsx')
dir='remove_duplication'
setwd(dir)
bam.files=list.files(dir)[grepl('.bam.gz$', list.files(dir))]
# find out specific marker files
h3.bam.files=bam.files[grepl('-3 H3 total', label$sample)]
k4me3.bam.files=bam.files[grepl('-3 H3K4me3', label$sample)]
k36me3.bam.files=bam.files[grepl('-3 H3K36me3', label$sample)]
k27ac.bam.files=bam.files[grepl('-3 H3K27ac', label$sample)]


output_dir = '/counts'
genernal_output_dir = '/Volumes/Data1/ChIP_seq_11_2018/csaw/THOR_method/'
prefix = 'H3K4me3'
control_prefix = 'H3_total'
bams = k27ac.bam.files
sample_names = c('Y3_H3K4me3','O3_H3K4me3')
input_bams = h3.bam.files
p5_peak = readbed('Y.H3K27ac.MUSICpeak.bed')
p15_peak = readbed('O.H3K27ac.MUSICpeak.bed')

# ----count------
# count
for(i in c(200,400,2000)){
    count_hist(bams, prefix, output_dir, width=i)
}

for(i in c(200,400,2000)){
    count_hist(input_bams, control_prefix, output_dir, width=i)
}

#f = file.path(output_dir,list.files(output_dir))
f = list.files(output_dir)
case = file.path(output_dir,f[grepl(paste0(prefix,'.*RData'), f)])
control = file.path(output_dir,f[grepl(paste0(control_prefix,'.*RData'), f)])

# get SES and TMM factor using 2k bin
case2k = case[grepl('2k',case)]
control2k = control[grepl('2k',control)]
load(case2k)
counts.case2k = counts
load(control2k)
counts.control2k = counts

ses_factor = get_SES(counts.case2k, counts.control2k, genernal_output_dir, prefix, sample_names)

tmm_sampleInfo = get_tmm(counts.case2k, counts.control2k, ses_factor,
                         genernal_output_dir, prefix, sample_names)
f
#
l=list()
for(i in 1:length(case)){
    print(case[i])
    l[[i]] = find_diffpeak(case[i], control[i+1], ses_factor, tmm_sampleInfo, sample_names,peak_region=p15_peak, reci=F)
}

final = myCombine(l)
# 2kb is to long for TPB binding sites. remove that
#final = final[final$width!=2000,]
dim(final)
write.csv(final, file.path(genernal_output_dir, paste0(sample_names[2],'_diffpeak.csv')), quote=F, row.names=F)

# reciprocal here means find Young diff peak
l=list()
for(i in 1:length(case)){
    print(case[i])
    l[[i]] = find_diffpeak(case[i], control[i+1], ses_factor, tmm_sampleInfo, sample_names,peak_region=p5_peak, reci=T)
}

final = myCombine(l)
# 2kb is to long for TPB binding sites. remove that
#final = final[final$width!=2000,]
dim(final)
write.csv(final, file.path(genernal_output_dir, paste0(sample_names[1],'_diffpeak.csv')), quote=F, row.names=F)
