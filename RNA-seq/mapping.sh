paired=$1
species=$2

if [ $species = 'human' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using human hg19 index'
	echo '---------------------------------------------'
	echo
	index='hg19'
elif [ $species = 'mouse' ];then
	echo
	echo '---------------------------------------------'
	echo 'Using mouse mm10 index'
	echo '---------------------------------------------'
	echo
	index='mm10'
else
	echo 'Usage:'
	echo 'mapping.sh <paired> <species>'
	echo 'species must be human or mouse'
	exit 1
fi





if [ $paired = 'F' ];then
	echo '---------------------------------------------'
	echo 'This is only used for single-end fastq files'
	echo '---------------------------------------------'
	echo ' '
	for i in *.fastq*;do
		if [[ $i =~ 'fastq'$ ]] || [[ $i =~ 'fastq.gz'$ ]];then
			echo $i
			if [ -f ${i/.fastq/_trimmed.fq} ] || [ -f ${i/.fastq/_trimmed.fq.gz} ];then
				echo 'trimmed file '${i/.fastq/_trimmed.fq}' already exist! skpping...'
			else
				trim_galore --gzip $i
				rm $i
			fi
		fi
	done

	echo '------------------'
	echo 'trimming finished'
	echo '------------------'

	for i in *_trimmed.fq*;do
		if [[ $i =~ '_trimmed.fq.gz'$ ]];then
			echo $i
			hisat2 -p 8 -x ${index} -U $i --summary-file ${i/_trimmed.fq.gz/mappingSummary.txt} -S ${i/_trimmed.fq.gz/.sam}
			if [ -f ${i/_trimmed.fq.gz/}.sam ]; then
				samtools sort -@ 2 ${i/_trimmed.fq.gz/.sam} -o ${i/_trimmed.fq.gz/.bam}
				samtools index ${i/_trimmed.fq.gz/.bam}
				bai_file=${i/_trimmed.fq.gz/.bam}.bai
				if [ -f $bai_file ]; then rm ${i/_trimmed.fq.gz/.sam}; fi
			fi
		fi
	done

	echo '------------------'
	echo 'mapping finished'
	echo '------------------'
elif [ $paired = 'T' ];then
	echo '---------------------------------------------'
	echo 'This is only used for paired-end fastq files'
	echo '---------------------------------------------'
	echo ' '
	for i in *.fastq*;do
		if [[ $i =~ 'fastq'$ ]] || [[ $i =~ 'fastq.gz'$ ]];then
			if [[ $i =~ '_1' ]];then
				echo $i
				x=${i/_1/_2}
				if [ -f ${i/_1.fastq/_1_val_1.fq} ] || [ -f ${i/_1.fastq/_1_val_1.fq.gz} ];then
					echo 'trimmed file '${i/_1.fastq/_1_val_1.fq.gz}' already exist! skpping...'
				else
					trim_galore --gzip --paired $i $x
					rm $i $x
				fi
			fi
		fi
	done

	echo '------------------'
	echo 'trimming finished'
	echo '------------------'

	for i in *_1_val_1.fq*;do
		if [[ $i =~ 'fq.gz'$ ]];then
			echo $i
			x=${i/_1_val_1.fq/_2_val_2.fq}
			hisat2 -p 8 -x ${index} -1 $i -2 $x --summary-file ${i/_1_val_1.fq.gz/mappingSummary.txt} -S ${i/_1_val_1.fq.gz/.sam}
			if [ -f ${i/_1_val_1.fq.gz/}.sam ]; then
				samtools sort -@ 2 ${i/_1_val_1.fq.gz/.sam} -o ${i/_1_val_1.fq.gz/.bam}
				samtools index ${i/_1_val_1.fq.gz/.bam}
				bai_file=${i/_1_val_1.fq.gz/.bam}.bai
				if [ -f $bai_file ]; then rm ${i/_1_val_1.fq.gz/.sam}; fi
			fi
		fi
	done

	echo '------------------'
	echo 'mapping finished'
	echo '------------------'
else
	echo 'paired must be "T" or "F"'
	exit 1
fi
