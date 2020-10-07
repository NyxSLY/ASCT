# single end strand specific read seperate

file=$1

samtools view -b -F 16 $file > ${file/.bam/.forward.bam}
samtools index ${file/.bam/.forward.bam}

samtools view -b -f 16 $file > ${file/.bam/.reverse.bam}
samtools index ${file/.bam/.reverse.bam}
