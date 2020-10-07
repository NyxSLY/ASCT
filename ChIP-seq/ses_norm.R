# SES normalization
# Normalization, bias correction, and peak calling for ChIP-seq
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342857/

# a package using this method
# Differential peak calling of ChIP-seq signals with replicates with THOR
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5175345/#B42


library(readxl)

# -----
# input bam
# -----


label=read_excel('sequencing\ libraries\ info.xlsx')

dir='remove_duplication'
setwd(dir)
bam.files=list.files(dir)[grepl('.bam.gz$', list.files(dir))]
h3.bam.files=bam.files[grepl('H3 total', label$sample)]
marker='H3K27ac'

# find out marker file
bam.files=bam.files[grepl(marker, label$sample)]
names = c(paste0('rep1_Y_',marker),paste0('rep1_O_',marker),paste0('rep3_Y_',marker), paste0('rep3_O_',marker))

# read SES factor
factor_dir='csaw/SES_norm'
f = file.path(factor_dir, paste0(tolower(marker),'.factor.txt'))
factor = read.table(f)



l=list()
for(i in 1:length(bam.files)){
    cmd <- paste0('bamCompare --outFileFormat bedgraph --scaleFactors 1:',factor$ses[i],' -p max --operation subtract --pseudocount 1 -bs 10 --smoothLength 30 --minMappingQuality 30 --blackListFileName blacklist.bed --extendReads --ignoreDuplicates -b1 ',
                  bam.files[i],' -b2 ',h3.bam.files[i],' -o signal_tracks/SES/',
                  names[i],'.ses.subtract.bedgraph')
    l[[i]] = cmd
}


file=paste0('run_',marker,'.sh')
sink(file)

for(i in 1:4){
    cat(l[[i]])
    cat("\n")
}


sink()
