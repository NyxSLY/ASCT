paired=$1
species=$2

if [ $species = 'human' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using human hg19 GTF'
	echo '---------------------------------------------'
	echo
	gtf_file='Human_database/hg19/gencode.v19.annotation.sorted.gtf'
elif [ $species = 'mouse' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using mouse mm10 GTF'
	echo '---------------------------------------------'
	echo
	gtf_file='mouse_database/Mus_musculus.GRCm38.80.gtf'
elif [ $species = 'mouseERCC' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using mouse mm10 GTF'
	echo '---------------------------------------------'
	echo
	gtf_file='mouse_database/Mus_musculus.GRCm38.80.ercc.gtf'
elif [ $species = 'Danio_rerio' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using zebrafish GRCz11 GTF'
	echo '---------------------------------------------'
	echo
	gtf_file='GENOMES/zebrafish/Danio_rerio.GRCz11.98.gtf'
elif [ $species = 'Caenorhabditis_elegans' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using Caenorhabditis_elegans ce10'
	echo '---------------------------------------------'
	echo
	gtf_file='GENOMES/worm_genome/Caenorhabditis_elegans.WBcel235.85.gtf'
else
	echo 'Usage:'
	echo 'featureCount.sh <paired> <species>'
	echo 'species must be human or mouse or zebrafish'
	exit 1
fi


if [ $paired = 'F' ];then
	base='featureCounts -t exon -g transcript_id -f -O -T 8 -a '${gtf_file}' -o featureCountExon.txt '
	for i in *.bam;do
		base=${base}${i}" "
	done
	echo $base
	$base
elif [ $paired = 'T' ];then
	base='featureCounts -t exon -g transcript_id -f -O -T 8 -p -C -a '${gtf_file}' -o featureCountExon.txt '
	for i in *.bam;do
		base=${base}${i}" "
	done
	echo $base
	$base
else
	echo 'paired must be "T" or "F"'
	exit 1
fi
