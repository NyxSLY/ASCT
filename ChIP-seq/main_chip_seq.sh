# trim and map

for i in *R1_001.fastq.gz;do

f2=${i/R1/R2}
echo $i
trim_galore --paired $i $f2
done

bowtie2 -p 6 -t -x /Volumes/LACIE/Human_database/hg19/bowtie2_index/hg19 -1 {} -2 {} -S {}.sam

for i in *R1_001_val_1.fq.gz;do
  f2=${i/R1_001_val_1/R2_001_val_2}
  echo $i
  bowtie2 -p 6 -t -x /Volumes/LACIE/Human_database/hg19/bowtie2_index/hg19 -1 $i -2 $f2 | samtools view -S -b -q 30 - | samtools sort - -@ 6 -o ${i/_R1_001_val_1.fq/.sorted.bam}
done

# get the number of mapped reads
sh samtools_flagstat.sh

# remove duplicates
python redup_then2bed.py .

# trim reads mapped to blacklist
sh trim_blacklist.sh

# find diff bound region
Rscript find_diff.R

# generate normalized signaltracks
bamCompare --outFileFormat bedgraph --scaleFactors <1:normfactor> -p max --operation subtract -bs 10 --minMappingQuality 30 --blackListFileName blacklist.bed --extendReads --ignoreDuplicates -b1 <chip.bam> -b2 <H3_total.bam> -o <ses.subtract.bedgraph>
python proc_tmm.py .

# chromHMM
sh chromHMM.sh
