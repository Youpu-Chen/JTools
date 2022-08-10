# tools for microRNome analysis

This repo is used to deposit some scripts I used in sRNA-related analysis

Chinese Version -> [Click here](./「我的microNome组学分析流程」第1版)

## (1) fastq2fasta.py

This script is used to convert fastq file (single-ended Illumina reads) to fasta format, but with additional demands.

The sequence id of fasta file will be reformatted like below:

`seq<serial number>_<frequency>_<length>`





---

Run the code below you should install `pandas` at first.

## (2) Get the abbreviation and sequence of specified species from miRBase

> Note: the "sequence" is hairpin RNA and mature microRNA 

the corresponding script is `miRBasehandle.py`

In this script, you should specify the abbreviation of specified species in the above script, then run the code below.

```python
latin_name = LatinNamematch('organisms.txt.gz')
Organismhairpin()
Organismmature()
```



## (3) Organize the bowtie default output and count the frequency of each microRNA

In order to calculate the frequency of each microRNA, we need to prepare a `.pos` file first,

- 1st column: hpRNA ID
- 2nd column: start position of mature microRNA
- 3rd column: end position of mature microRNA
- 4th column: mature microRNA ID

Using `blasthandle.py` to prepare this file:

> Note: You should use mature microRNA as input to blast against the hpRNA (using specified species)

Example code,

```python
Blastfilter('osa_mature_hairpin_blast.txt') 
mirnaMatch('Filtered_osa_mature_hairpin_blast.txt')
Getpos('osa', 'Matched_Filtered_osa_mature_hairpin_blast.txt')
```

> Note: `osa_mature_hairpin_blast.txt` is the result of blast analysis.
>
> Using function `Getpos` to generate the `.pos` file

After you get the `.pos` file, you can now run the `organizebowtie.py` to get the summary statistics of you microRNA alignment analysis.

```python
premature2start, premature2end, premature2mature = PostionGet('osa_hairpin_mature.pos')

bowtieoutput_list = BowtieOutputNameGet('bowtie_output')
print(bowtieoutput_list)

# tmp_df = OneBowtieOutputOrganize('bowtie_output/32RT-1.2207091646.sam')
# print(tmp_df)
# tmp_df.to_csv('./tmp.32RT-1.csv')

tmp_df = MultipleBowtieOuputOrganize(bowtieoutput_list)
tmp_df.to_csv('all_samples.csv')
print(tmp_df)
```

> Note: 
>
> - `OneBowtieOutputOrganize` function is used to generate one sample output of bowtie alignment.
>
> In this analysis, you should move all the bowtie results in one dir, which named as `bowtie_output`
>
>  Then using `BowtieOutputNameGet` to get all the filename.
>
> Finally, using `MultipleBowtieOuputOrganize` to generate the summary statistics table.



## (4) 「Pandas Version」Organize the bowtie default output and count the frequency of each microRNA

Using the `pandas_organizebowtie.py`, We can now output the same kind of table as R.

First, using the `.pos` file as well.

Second, using `BowtieOutputNameGet` to get the bowtie output filename.

Third, using `MultipleBowtieOuputOrganize` to generate raw count table, which is set to `pandas_all_samples.csv` by default.



## (5) Convert the raw count table to RPM normalized table

> Note: the bwt output files used in this process are not the same as before.
>
> In this part, we should use the bwt output which generated from the alignment against your specified genome.

The script relates to this part is `RPM.py`.

First, using `BowtieOutputNameGet` to generate the bwt output filename list.

Second, using `DFRPMconvert` to get the RPM-normalized table (`RPM.normalized.csv` by default).



## (6) miRNA frequency visualization

1) Using the script `sample_freqplot.R` you could complete the basic visualization of the frequency of different length of miRNA in one sample.

All you need to is to change the specified name in the script.

## (7) DESeq2 analysis

The corresponding script is `DESeq.R`

## (8) Heatmap visualization

the script `pheatmap.R` and `ComplexHeatmap.R` are on the way ~



# To be continued

