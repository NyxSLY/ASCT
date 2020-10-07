paired=$1
species=$2

if [ $species = 'human' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using human hg19 index'
	echo '---------------------------------------------'
	echo
	index='/Volumes/LACIE/Human_database/hg19/salmon/hg19_index'
elif [ $species = 'mouse' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using mouse mm10 index'
	echo '---------------------------------------------'
	echo
	index='/Volumes/LACIE/mouse_database/salmon/mm10_index'
else
	echo 'Usage:'
	echo 'featureCount.sh <paired> <species>'
	echo 'species must be human or mouse'
	exit 1
fi

if [ $paired = 'F' ];then
	for i in *.fq.gz;do
		salmon quant -i ${index} -l A -r $i -o ${i/_trimmed.fq.gz/_salmon} -p 8 --validateMappings
	done
elif [ $paired = 'T' ];then
	for i in *_1_val_1.fq*;do
		echo $i
		if [[ $i =~ 'fq.gz'$ ]];then
			x=${i/_1_val_1.fq/_2_val_2.fq}
			salmon quant -i ${index} -l A -1 $i -2 $x -o ${i/_1_val_1.fq.gz/_salmon} -p 8 --validateMappings
		fi

		if [[ $i =~ 'fq'$ ]];then
			echo $i
			x=${i/_1_val_1.fq/_2_val_2.fq}
			salmon quant -i ${index} -l A -1 $i -2 $x -o ${i/_1_val_1.fq/_salmon} -p 8 --validateMappings
		fi
	done
else
	echo 'paired must be "T" or "F"'
	exit 1
fi	