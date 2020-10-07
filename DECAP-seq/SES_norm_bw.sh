# SES normliazation by deeptools
# generate bedgraph file first. May need to use TMM factor later

# p15
p='star_mapping/filter/'

file='star_P15Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
ctrl='star_P15_negAligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
bamCompare --scaleFactors 1:1.91 --skipNonCoveredRegions -b1 $p$file -b2 $p$ctrl --outFileName ${file/.bam/.SES.bedgraph} -of bedgraph --operation subtract --binSize 10 --minMappingQuality 30 --skipZeroOverZero -p 10 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed
file='star_P15Aligned.sorted.nomulti.onlyR1.rev.filtered.bam'
ctrl='star_P15_negAligned.sorted.nomulti.onlyR1.rev.filtered.bam'
bamCompare --scaleFactors 1:1.88 --skipNonCoveredRegions -b1 $p$file -b2 $p$ctrl --outFileName ${file/.bam/.SES.bedgraph} -of bedgraph --operation subtract --binSize 10 --minMappingQuality 30 --skipZeroOverZero -p 10 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

# p5
file='star_P5Aligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
ctrl='star_P15_negAligned.sorted.nomulti.onlyR1.fwd.filtered.bam'
bamCompare --scaleFactors 1:1.61 --skipNonCoveredRegions -b1 $p$file -b2 $p$ctrl --outFileName ${file/.bam/.SES.bedgraph} -of bedgraph --operation subtract --binSize 10 --minMappingQuality 30 --skipZeroOverZero -p 10 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed

file='star_P5Aligned.sorted.nomulti.onlyR1.rev.filtered.bam'
ctrl='star_P15_negAligned.sorted.nomulti.onlyR1.rev.filtered.bam'
bamCompare --scaleFactors 1:1.59 --skipNonCoveredRegions -b1 $p$file -b2 $p$ctrl --outFileName ${file/.bam/.SES.bedgraph} -of bedgraph --operation subtract --binSize 10 --minMappingQuality 30 --skipZeroOverZero -p 10 --blackListFileName /Volumes/LACIE/Human_database/hg19/blacklist.bed
