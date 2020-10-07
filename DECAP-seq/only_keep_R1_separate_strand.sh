## this script only keep the R1 reads
## the name of fwd and res may need to be switched depend on the protocol of strand specific library constraction
set -ue

# Get the bam file from the command line
DATA=$1

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#

samtools view -b -f 64 -F 32 $DATA > ${DATA/.bam/.onlyR1.rev.bam}
samtools index ${DATA/.bam/.onlyR1.rev.bam}



# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#


samtools view -b -f 96 $DATA > ${DATA/.bam/.onlyR1.fwd.bam}
samtools index ${DATA/.bam/.onlyR1.rev.bam}

#
# Combine alignments that originate on the reverse strand.
#

# delete temp files
