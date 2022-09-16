# README

`Evolution`



# OrthoFinder2

## `ConvertOrthogroup2geneID.py`

this script parses the `Orthogroup.tsv` to a gene list file, a single line contains a gene id.

```shell
python3 ConvertOrthogroup2geneID.py Orthogroup.tsv rewrite.Orthogroups.tsv
```



## `ParseN0_Ortho.py`

`ParseN0_Ortho.py` could be used to parse `N0.tsv`, which is a splited Orthogroup text file (another version of `Orthogroup.tsv`), records the orthogroup info, into a gene list.

Each line of gene list contains a orthogroup ID and a gene ID.

```shell
python3 ParseN0_Ortho.py N0.tsv ortholog_list.tsv
```



## `count_organize.R` + `Statistics.py`

`count_organize.R` is a R script, which reads gene list generated from `ParseN0_Ortho.py` and count matrix from featureCounts, then output `Hierarchical_Orthogroups.Count.csv`,

```shell
Orthogroup,A,B,C
N0.HOG0000000,NA,0,0
N0.HOG0000001,NA,0,NA
N0.HOG0000002,NA,0,NA
```

and `Statistics.py` could be used to do the summary statistics,

```shell
python3 Statistics.py Hierarchical_Orthogroups.Count.csv Group_count.txt
```



## `SeekOrthogroup.py`

This script take a set of blast results to infer the gene belonging relationships.

Let me introduce a little bit about the background: let's say that there is a flow-development gene which is inclued not in the OrthoFinder2, and we want to find the gene's possible orthology,

First we ran blast, to retrieve the best-hit of that gene in each species, and using `SeekOrthogroup.py` to descripted script.

```shell
cat input_seq.config | while read id   # input_seq.config is the input seq configuration file
do
	tra1=${id#*/}
	input=${tra1%.*}
	# echo $input
	mkdir blast_results/${input}_blast_results 
	cat db.config | while read db  # db.config is the prefix of blast database
	do
		echo "blastp -num_threads 16 -query $id -db blast_database/$db -out blast_results/${input}_blast_results/${input}_${db}.ofmt6.evalue5.txt -outfmt 6 -evalue 1e-5" >> ${input}.blast.commands
	done
done

cat A*.commands >> all.commands
ParaFly -c all.commands -CPU 8 > 220916.blastp.err.log 2>&1 &


cat dir.config | while read id
do
	python3 SeekOrthogroup.py $id ../../Orthogroups.tsv >> Arab_Ortholog.txt
done
```

`Arab_Ortholog.txt` is the final output.



# PhyloMCL

## `geneid_length.py`

this script is used to prepare the `<gene_id>\t<gene_length> ` file for PhyloMCL analysis.



# General

## `FileterOrtholog.py`

this script is used to filter the <font color='yellow'>false-positive</font> orthogroup, 

Assuming no chromosomal arrangement, it means that the ortholog from A and B, it will not appear on the chromosome C.

So, using a `chr_pair.txt` file as input, and parse it as filtering conditions, it'll be able to filter the incorrectly inferred orthologs.

```shell
python3 FileterOrtholog.py chr_pair.txt Ortholog_fasta_directory
```

> Note: it hasn't been adopted into general usage.

Here is an example of `chr_pair.txt`,

```shell
A1	B3	C3
A2	B17	C4
A3	B1	C6
A4	B5	NA
A5	B4	NA
```

In the case above, let's assume the species is hexaploid, so one haplotype is going to have 3 subgenome,

When 2 subgenomes have no significant synteny relation ships with the other one, which means the evolution info is lost in one subgenome, it is defined as NA.



## `GrepOrthologFromGFF.py`

this script takes `rewrite.Orthogroups.tsv` as filtering conditions, to grep corresponding gene id, and re-write it into `BED`, which then could be used as input of [seqkit](https://bioinf.shenwei.me/seqkit/).

```shell
python3 GrepOrthologFromGFF.py rewrite.Orthogroups.tsv input.gff3 output.bed
```

