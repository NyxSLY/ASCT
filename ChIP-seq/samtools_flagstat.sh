for bam in *bam.gz;do
samtools flagstat $bam > ${bam/.bam.gz/.flagstat.txt}
done