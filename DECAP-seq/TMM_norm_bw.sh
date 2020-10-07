# scaleFactro should be 1/TMM factor

p='star_mapping/filter/'
file='star_P15Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 0.83 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P15Aligned.sorted.nomulti.onlyR1.rev.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 0.83 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P5Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 0.90 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P5Aligned.sorted.nomulti.onlyR1.rev.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 0.90 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P15_negAligned.sorted.nomulti.onlyR1.rev.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 1.46 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P15_negAligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
bamCoverage -b $p$file --outFileName ${file/.bam/.TMM.bw} -of bigwig --scaleFactor 1.46 --binSize 10 --minMappingQuality 30 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed
