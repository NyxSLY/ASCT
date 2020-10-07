library(GenomicRanges)
library(GenomicFeatures)


# read bed file and convert to genomic ranges
readbed <- function(file) {
    data <- read.table(file, sep = "\t", stringsAsFactors=F)
    if(!grepl('chr', data$V1[1])){
        data$V1 = paste0('chr', data$V1)
    }
    dat <- makeGRangesFromDataFrame(data,
        seqnames.field = "V1",
        start.field = "V2",
        end.field = "V3",
        keep.extra.columns = T
    )
    return(dat)
}

# load hg19 annotation
loadhg19 <- function() {
    gtff <- "/Volumes/LACIE/Human_database/hg19/gencode.v19.annotation.gtf"
    gtf <- makeTxDbFromGFF(gtff, format = "gtf")
    pro <- promoters(gtf, upstream = 2000, downstream = 2000)
    gene <- genes(gtf)
    utr3 <- threeUTRsByTranscript(gtf)
    return(list(gtf=gtf, pro=pro, gene=gene, utr3=utr3))
}


toGR <- function(data){
    dat = makeGRangesFromDataFrame(data, seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3',
                                   keep.extra.columns = T)
    return(dat)
}
