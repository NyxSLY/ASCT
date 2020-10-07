library(readxl)


label=read_excel('/Volumes/Data1/ChIP_seq_11_2018/sequencing\ libraries\ info.xlsx')

dir='/Volumes/Data1/ChIP_seq_11_2018/bam/remove_duplication'
setwd(dir)
bam.files=list.files(dir)[grepl('.bam.gz$', list.files(dir))]

# find out specific marker files
h3.bam.files=bam.files[grepl('H3 total', label$sample)]
k4me1.bam.files=bam.files[grepl('H3K4me1', label$sample)]
k27ac.bam.files=bam.files[grepl('H3K27ac', label$sample)]
k4me3.bam.files=bam.files[grepl('H3K4me3', label$sample)]
k36me3.bam.files=bam.files[grepl('H3K36me3', label$sample)]



generate_code = function(k4me1.bam.files, h3.bam.files, his, dir){
    setwd('/Volumes/Data1/ChIP_seq_11_2018/MUSIC')
    punctate_list = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H2az', 'PolII')
    broad_list = c('H3K9me3', 'H3K36me3', 'H3K27me3', 'H3K79me2', 'H4K20me1')
    
    if(his %in% punctate_list){
        model='get_optimal_punctate_ERs'
    } else if(his %in% broad_list){
        model='get_optimal_broad_ERs'
    } else{
        stop('Something is wrong. Maybe this is a new histome marker, do not know which model should be used')
    }
    names(k4me1.bam.files) = c('rep1_Y','rep1_O','rep3_Y', 'rep3_O')
    names(h3.bam.files) = c('rep1_Y','rep1_O','rep3_Y', 'rep3_O')
    
    for(i in 1:length(k4me1.bam.files)){
        case = file.path(dir, k4me1.bam.files[i])
        ctl = file.path(dir, h3.bam.files[i])
        prefix = paste(names(k4me1.bam.files)[i], his,sep='.')
        cmd1 = paste0('mkdir ',prefix,';cd ',prefix)
        cmd2 = paste0('mkdir chip;mkdir control')
        cmd3 = paste0('samtools view ',case,' | MUSIC -preprocess SAM stdin chip/')
        cmd4 = paste0('samtools view ',ctl,' | MUSIC -preprocess SAM stdin control/')
        cmd5 = 'mkdir chip/sorted;mkdir chip/dedup;mkdir control/sorted;mkdir control/dedup
MUSIC -sort_reads chip chip/sorted
MUSIC -sort_reads control control/sorted
MUSIC -remove_duplicates chip/sorted 2 chip/dedup
MUSIC -remove_duplicates control/sorted 2 control/dedup
cd chip/dedup;rm KI*;rm GL*
head -25 chr_ids.txt > chr_ids.txt1;mv chr_ids.txt1 chr_ids.txt
cd ../../control/dedup;rm KI*;rm GL*
cd ../..'
        cmd6 = paste0('run_MUSIC.csh -',model,' ./chip/dedup ./control/dedup /Volumes/LACIE/Human_database/hg19/multi_mappability_100
cd ..')
        fileConn<-file(paste0(prefix,'.sh'))
        writeLines(c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6), fileConn)
        close(fileConn)
    }
}

generate_code(k4me3.bam.files, h3.bam.files, 'H3K4me3')
generate_code(k4me1.bam.files, h3.bam.files, 'H3K4me1')
generate_code(k27ac.bam.files, h3.bam.files, 'H3K27ac')
generate_code(k36me3.bam.files, h3.bam.files, 'H3K36me3')


#--------------------------------

# from old chip-seq data
dir1='/Volumes/luyang/Histone_Modification_ChIP_seq/hg19/clean_bam'
bam.files1 = list.files(dir1)[grepl('.bam$', list.files(dir1))]
# filter rep2 data
bam.files1 = bam.files1[!grepl('2_', bam.files1)]
k27me3.bam.files = bam.files1[grepl('H3K27me3', bam.files1)]
k9me3.bam.files = bam.files1[grepl('H3K9me3', bam.files1)][c(1,3)]  # Y3 is missing
h3.old.bam.files = bam.files1[grepl('H3_total', bam.files1)]
input.old.bam.files =  bam.files1[grepl('input', bam.files1)]


generate_code_k27me3 = function(k4me1.bam.files, h3.bam.files, his, dir){
    setwd('/Volumes/Data1/ChIP_seq_11_2018/MUSIC')
    punctate_list = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H2az', 'PolII')
    broad_list = c('H3K9me3', 'H3K36me3', 'H3K27me3', 'H3K79me2', 'H4K20me1')
    
    if(his %in% punctate_list){
        model='get_optimal_punctate_ERs'
    } else if(his %in% broad_list){
        model='get_optimal_broad_ERs'
    } else{
        stop('Something is wrong. Maybe this is a new histome marker, do not know which model should be used')
    }
    names(k4me1.bam.files) = c('rep1_O','rep3_O','rep1_Y', 'rep3_Y')
    names(h3.bam.files) = c('rep1_O','rep3_O','rep1_Y', 'rep3_Y')
    
    for(i in 1:length(k4me1.bam.files)){
        case = file.path(dir, k4me1.bam.files[i])
        ctl = file.path(dir, h3.bam.files[i])
        prefix = paste(names(k4me1.bam.files)[i], his,sep='.')
        cmd1 = paste0('mkdir ',prefix,';cd ',prefix)
        cmd2 = paste0('mkdir chip;mkdir control')
        cmd3 = paste0('samtools view ',case,' | MUSIC -preprocess SAM stdin chip/')
        cmd4 = paste0('samtools view ',ctl,' | MUSIC -preprocess SAM stdin control/')
        cmd5 = 'mkdir chip/sorted;mkdir chip/dedup;mkdir control/sorted;mkdir control/dedup
MUSIC -sort_reads chip chip/sorted
MUSIC -sort_reads control control/sorted
MUSIC -remove_duplicates chip/sorted 2 chip/dedup
MUSIC -remove_duplicates control/sorted 2 control/dedup
cd chip/dedup;rm KI*;rm GL*
head -25 chr_ids.txt > chr_ids.txt1;mv chr_ids.txt1 chr_ids.txt
cd ../../control/dedup;rm KI*;rm GL*
cd ../..'
        cmd6 = paste0('run_MUSIC.csh -',model,' ./chip/dedup ./control/dedup /Volumes/LACIE/Human_database/hg19/multi_mappability_100
cd ..')
        fileConn<-file(paste0(prefix,'.sh'))
        writeLines(c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6), fileConn)
        close(fileConn)
    }
}

generate_code_k27me3(k27me3.bam.files, input.old.bam.files, 'H3K27me3', dir1)

generate_code_k9me3 = function(k4me1.bam.files, h3.bam.files, his, dir){
    setwd('/Volumes/Data1/ChIP_seq_11_2018/MUSIC')
    punctate_list = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H2az', 'PolII')
    broad_list = c('H3K9me3', 'H3K36me3', 'H3K27me3', 'H3K79me2', 'H4K20me1')
    
    if(his %in% punctate_list){
        model='get_optimal_punctate_ERs'
    } else if(his %in% broad_list){
        model='get_optimal_broad_ERs'
    } else{
        stop('Something is wrong. Maybe this is a new histome marker, do not know which model should be used')
    }
    names(k4me1.bam.files) = c('rep1_O','rep1_Y')
    h3.bam.files = h3.bam.files[c(1,3)]
    names(h3.bam.files) = c('rep1_O','rep1_Y')
    
    for(i in 1:length(k4me1.bam.files)){
        case = file.path(dir, k4me1.bam.files[i])
        ctl = file.path(dir, h3.bam.files[i])
        prefix = paste(names(k4me1.bam.files)[i], his,sep='.')
        cmd1 = paste0('mkdir ',prefix,';cd ',prefix)
        cmd2 = paste0('mkdir chip;mkdir control')
        cmd3 = paste0('samtools view ',case,' | MUSIC -preprocess SAM stdin chip/')
        cmd4 = paste0('samtools view ',ctl,' | MUSIC -preprocess SAM stdin control/')
        cmd5 = 'mkdir chip/sorted;mkdir chip/dedup;mkdir control/sorted;mkdir control/dedup
MUSIC -sort_reads chip chip/sorted
MUSIC -sort_reads control control/sorted
MUSIC -remove_duplicates chip/sorted 2 chip/dedup
MUSIC -remove_duplicates control/sorted 2 control/dedup
cd chip/dedup;rm KI*;rm GL*
head -25 chr_ids.txt > chr_ids.txt1;mv chr_ids.txt1 chr_ids.txt
cd ../../control/dedup;rm KI*;rm GL*
cd ../..'
        cmd6 = paste0('run_MUSIC.csh -',model,' ./chip/dedup ./control/dedup /Volumes/LACIE/Human_database/hg19/multi_mappability_100
cd ..')
        fileConn<-file(paste0(prefix,'.sh'))
        writeLines(c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6), fileConn)
        close(fileConn)
    }
}
generate_code_k9me3(k9me3.bam.files, input.old.bam.files, 'H3K9me3', dir1)


