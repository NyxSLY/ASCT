# trim adapter
trim_galore --paired --small_rna ../P15_TL_seq.read1 ../P15_TL_seq.read2

# mapping
STAR --runThreadN 12 --genomeDir /home/ec2-user/STAR_Genome/hg19 --readFilesIn /home/ec2-user/TL-seq/P15_TL_seq.read1_val_1.fq /home/ec2-user/TL-seq/P15_TL_seq.read2_val_2.fq --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix /home/ec2-user/TL-seq/star_mapping/star_P15 --outSAMtype BAM Unsorted

# sort
samtools view -@ 10 -q 30 -b paired_trim.sam | samtools sort - paired_trim.nomulti

# kepp 1bp
for file in *nomulti.bam
do
  echo 'start to separate strand'
  only_keep_R1_separate_strand.sh $file

  samtools view -m 30 -f 2 -b ${file/.bam/.onlyR1.fwd.bam} > ${file/.bam/.onlyR1.fwd.filtered.bam}
  samtools view -m 30 -f 2 -b ${file/.bam/.onlyR1.rev.bam} > ${file/.bam/.onlyR1.rev.filtered.bam}

  # convert to bedgraph
  bedtools genomecov -ibam ${file/.bam/.onlyR1.fwd.filtered.bam} -bg -split > ${file/.bam/.onlyR1.fwd.filtered.bedgraph}
  bedtools genomecov -ibam ${file/.bam/.onlyR1.rev.filtered.bam} -bg -split > ${file/.bam/.onlyR1.rev.filtered.bedgraph}

  # convert to bw
  bedGraphToBigWig ${file/.bam/.onlyR1.fwd.filtered.bedgraph} /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt ${file/.bam/.onlyR1.fwd.filtered.bw}
  bedGraphToBigWig ${file/.bam/.onlyR1.rev.filtered.bedgraph} /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt ${file/.bam/.onlyR1.rev.filtered.bw}

  # keep the first bp
  # convert to bed file
  echo 'start to convert to bed file'
  bedtools bamtobed -i ${file/.bam/.onlyR1.fwd.filtered.bam} > ${file/.bam/.onlyR1.fwd.filtered.bed}
  bedtools bamtobed -i ${file/.bam/.onlyR1.rev.filtered.bam} > ${file/.bam/.onlyR1.rev.filtered.bed}

  # keep the first bp in the fwd bam and the last bp in the rev bam
  echo 'start to keep the first bp'
  awk '$3=$2' ${file/.bam/.onlyR1.fwd.filtered.bed} | sed 's/ /\t/g' > ${file/.bam/.onlyR1.fwd.filtered.1bp.bed}
  awk '$2=$3' ${file/.bam/.onlyR1.rev.filtered.bed} | sed 's/ /\t/g' > ${file/.bam/.onlyR1.rev.filtered.1bp.bed}

  # convert to bedgraph
  bedtools genomecov -i ${file/.bam/.onlyR1.fwd.filtered.1bp.bed} -g /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt -bg > ${file/.bam/.onlyR1.fwd.filtered.1bp.bedgraph}
  bedtools genomecov -i ${file/.bam/.onlyR1.rev.filtered.1bp.bed} -g /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt -bg > ${file/.bam/.onlyR1.rev.filtered.1bp.bedgraph}

  # convert to bw
  bedGraphToBigWig ${file/.bam/.onlyR1.fwd.filtered.1bp.bedgraph} /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt ${file/.bam/.onlyR1.fwd.filtered.1bp.bw}
  bedGraphToBigWig ${file/.bam/.onlyR1.rev.filtered.1bp.bedgraph} /Volumes/LACIE/Human_database/hg19/hg19.chrom.sizes.txt ${file/.bam/.onlyR1.rev.filtered.1bp.bw}
  echo 'All finished!'
  date
done

# SES and TMM normalization - calculate factors
Rscript ses_normalization.R

# normalization by SES and TMM factors - facctors are hardcoded here
sh SES_norm_bw.sh
sh TMM_norm_bw.sh

# identify age associated CT website
Rscript find_age_associatedCT.R
