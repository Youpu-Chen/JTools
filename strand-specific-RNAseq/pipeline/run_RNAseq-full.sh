#!/bin/bash
# Define the input dir
raw_data="00.raw_data"                # path to raw sequencing datasets
seq=$1                                # input seq
name=${seq%.*}                        # prefix
anno_data="00.anno"                   # path to gtf file
clean_data="01.clean_data"
ref_data="02.refgenome_bowtieindex"   # path to reference genome and index of bowtie 1.0
align_data="02.alignments"
test_data="0x.test"
count_data="03.featureCounts"

# Triming
echo "--------------------------Now checking the existence of raw data directory--------------------------"
if [ -d $raw_data ]
then
	echo "Checked!"
	mkdir "$clean_data"
	echo "--------------------------Now create the directory of clean data directory--------------------------"
	echo "1) trim_galore and fastp is applied to trim the adapter as well as the low-quality sequences."
	cmd1="The command is:\ntrim_galore -j 8 -o $clean_data -q 15 --phred33 ${raw_data}/$1 2>${clean_data}/trim-galore.log"
	echo "1.1) trim-galore mode"
	echo -e "$cmd1"
	trim_galore -j 8 -o $clean_data -q 30 --phred33 ${raw_data}/$1 2>${clean_data}/trim_galore.log

	cmd2="The command is:\nfastp -w 64 -i ${raw_data}/$1 -o ${clean_data}/${name}.fastp.fastq -q 15 --n_base_limit 3 -j ${clean_data}/${name}.json -h ${clean_data}/${name}.html >${clean_data}/fastp.log 2>&1"
	echo "1.2) fastp mode"
	fastp -w 64 -i ${raw_data}/$1 -o ${clean_data}/${name}.fastp.fastq -q 15 --n_base_limit 3 -j ${clean_data}/${name}.json -h ${clean_data}/${name}.html >${clean_data}/fastp.log 2>&1
fi
touch trim-galore.ok
touch fastp.ok
echo "--------------------------Now triming is done!--------------------------"


# Alignment
echo "--------------------------Now conducting the alignment using bowtie 1.0--------------------------"
if [ -e "fastp.ok" ]
then
	echo "--------------------------Now checking the index file of bowtie 1.0--------------------------"
	echo "Checked!"
	echo "--------------------------Now create the alignment directory--------------------------"
	mkdir "$align_data"
	echo "2) bowtie 1.0 for the alignment."
	echo "2.1) RNA-Seq mode"
	cmd1="The command is:\nbowtie --sam -p 64 -v 3 -k 1 --best -q ${ref_data}/genome ${clean_data}/${name}.fastp.fastq ${align_data}/${name}.unsorted.sam 2>${align_data}/bowtie.log"
	echo -e "$cmd1"
	bowtie --sam -p 64 -v 3 -k 1 --best -q ${ref_data}/genome ${clean_data}/${name}.fastp.fastq ${align_data}/${name}.unsorted.sam 2>${align_data}/bowtie.log
fi
touch bowtie.ok
echo "--------------------------Now alignmen is done!--------------------------"


# Sort and index
echo "--------------------------Now using samtools to sort the output SAM file and bedtools to retrieve the read depth--------------------------"
if [ -e "bowtie.ok" ]
then
	echo "--------------------------Now sorting the SAM file--------------------------"
	echo "3) samtools 1.9 for SAM sorting"
	cmd1="The command is:\nsamtools sort -@ 64 -m 4G -O sam -o ${align_data}/${name}.sorted.sam ${align_data}/${name}.unsorted.sam"
	echo -e "${cmd1}"
	samtools sort -@ 64 -m 4G -O sam -o ${align_data}/${name}.sorted.sam ${align_data}/${name}.unsorted.sam

	echo "4) bedtoosl for retrieveing read depth"
	cmd2="The command is hidden. Details in run_RNAseq-full.sh"
	echo "${cmd2}"
	# 
	samtools faidx ${ref_data}/genome.fa
	awk '{print $1"\t"$2}' ${ref_data}/genome.fa.fai > ${ref_data}/genome.txt
	bedtools makewindows -g ${ref_data}/genome.txt -w 1000 > ${ref_data}/windows.bed
	
	# 
	samtools view -f 16 -b ${align_data}/${name}.sorted.sam > ${align_data}/${name}.sorted.rev.bam
	samtools view -F 16 -b ${align_data}/${name}.sorted.sam > ${align_data}/${name}.sorted.fwd.bam
	bedtools coverage -a ${ref_data}/windows.bed -b ${align_data}/${name}.sorted.rev.bam > ${align_data}/rev.depth.txt
	bedtools coverage -a ${ref_data}/windows.bed -b ${align_data}/${name}.sorted.fwd.bam > ${align_data}/fwd.depth.txt
fi
touch sort.ok
touch read_depth.ok
echo "--------------------------Now sorting and retrieving read depth is done!--------------------------"


# Quantification
echo "--------------------------Now conducting expression quantification using featureCounts--------------------------"
if [ -e "read_depth.ok" ]
then
mkdir $count_data
cmd="
featureCounts \
-T 16 \
-t exon \
-g gene_id \
-a ${anno_data}/sacCer2.gtf \
-o ${count_data}/rna.mat \
${align_data}/${name}.sorted.sam 2>${count_data}/featureCounts.log
"
echo -e $cmd

featureCounts \
-T 16 \
-t exon \
-g gene_id \
-a ${anno_data}/sacCer2.gtf \
-o ${count_data}/rna.mat \
${align_data}/${name}.sorted.sam 2>&1 >${count_data}/quan_rna.og.log
fi
touch count.ok
echo "--------------------------Now count is done!--------------------------"