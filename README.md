# tools for sRNA analysis

This repo is used to deposit some scripts I used in sRNA-related analysis



## (1) v1.fastq2fasta.py

This script is used to convert fastq file (single-ended Illumina reads) to fasta format, but with additional demands.

The sequence id of fasta file will be reformatted like below:

`seq<serial number>_<frequency>_<length>`





---

Run the code below you should install `pandas` at first.

## (2) Get the Latin name and sequence of specified species from miRBase

`miRBasehandle.py`

In this script, you should specify the abbreviation of specified species in the above script, then run the code below.

```python
latin_name = LatinNamematch('organisms.txt.gz')
Organismhairpin()
Organismmature()
```



## (3) Organize the default output of bowtie and count the frequency of each microRNA

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
> In this analysis, you should move all the bowtie results in one dir, which named after ``
>
>  Then using `BowtieOutputNameGet` to get all the filename.
>
> Finally, using `MultipleBowtieOuputOrganize` to generate the summary statistics table.



