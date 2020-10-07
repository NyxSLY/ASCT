for i in *.redup.bam.gz;
do
	bedtools intersect -v -abam $i -b /Volumes/LACIE/Human_database/hg19/blacklist.bed > ${i/.bam/.blacklistTrimmed.bam}
done
