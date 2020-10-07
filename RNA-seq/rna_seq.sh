# workflow

# mapping
# indexed genome is under hisat2 default path
sh mapping.sh <if_paired_end_data> <species>

# transcript-level read count:
# index is made from hg19 or mm10
sh salmon.sh <if_paired_end_data> <species>

# exon-level read count:
# gtf file need to be download genecode or ensembl website
sh featureCount.sh <if_paired_end_data> <species>

# downstream analysis:
command:
Rscript downstream.R <output_dir> <featureCountExon> <annofile> <species>

#annofile example
run	cellline	group
SRR3659172	KNS42	O
SRR3659174	KNS42	O
SRR3659176	KNS42	O
SRR3659178	SF188	Y
SRR3659180	SF188	Y
SRR3659182	SF188	Y
