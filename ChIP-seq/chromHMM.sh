# 6 markers model from bam files

java -jar ChromHMM.jar BinarizeBam hg19.chrom.sizes.txt remove_duplication sample_table.txt binarizedbed

java -jar ChromHMM.jar LearnModel binarizedbed learnModel_6marker_Reads_donor3 10 hg19
