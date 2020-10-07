# RNA-seq analysis workflow

## Mapping
```shell
sh mapping.sh <if_paired_end_data> <species>
```

## Transcript-level read count:
index is made from hg19 or mm10

```shell
sh salmon.sh <if_paired_end_data> <species>
```
## Exon-level read count:
gtf file need to be download genecode or ensembl website

```shell
sh featureCount.sh <if_paired_end_data> <species>
```
## Downstream analysis:
```shell
Rscript downstream.R <output_dir> <featureCountExon> <annofile> <species>
```

### anno file example
```
run	cellline	group
SRR3659172	KNS42	O
SRR3659174	KNS42	O
SRR3659176	KNS42	O
SRR3659178	SF188	Y
SRR3659180	SF188	Y
SRR3659182	SF188	Y
```
