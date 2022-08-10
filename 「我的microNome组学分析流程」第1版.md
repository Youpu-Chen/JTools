> Note：本篇文档用于总结来到吴老师实验室之后开展的一系列microRNA组学分析流程

# Upstream Analysis

## （1）miRNA-Seq

miRNA-Seq得到的原始测序数据，为SE（single-ended）测序数据，reads较短，一般在50 bp左右。

```shell
seqkit stats input.fa
# file       format  type   num_seqs      sum_len  min_len  avg_len  max_len
# input.fa   FASTA   DNA   8,072,734  191,145,079       18     23.7       45
```



## （2）fq转fa

1）fastq格式的数据如下，

```shell
@V350036046L1C001R0010000005/1_U1
TTTGGATTGAAGGGAGCTCTGTGGAAGGGGCATGCAGAGGAG
+
dedddefceeedbd[fdcfeeccceda[eeefdeeeedMeef
```

2）需要转换成fa格式，

fasta序列的seq ID需根据原始测序数据进行整理，格式为：`seq<order>_<frequency>_<length>`

使用脚本：`fastq2fasta.py`

使用方法：

```shell
python fastq2fasta.py --input input.fastq --output outputprefix.fa --gzip
```

> Note：该工作开展于目录`01.raw_data`

## （3）Build Index

需要进行后续的bowtie比对以及数据过滤，需要提前下载对应类型的序列 & 构建索引。

- 水稻非编码RNA序列
- 水稻参考基因组序列（reference genome）
- 水稻hairpin RNA序列（miRBase）

### 1、构建水稻非编码RNA序列索引

数据来源：`Ensembl Plants`

```shell
# download
wget -c http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/oryza_sativa/ncrna/Oryza_sativa.IRGSP-1.0.ncrna.fa.gz
gunzip *
```

提取对应类型的序列以及相关操作如下，

```shell
less -S Oryza_sativa.IRGSP-1.0.ncrna.fa | grep ">" | awk '{print $5}' | uniq
# gene_biotype:snRNA
# gene_biotype:RNase_MRP_RNA
# gene_biotype:rRNA
# gene_biotype:tRNA
# gene_biotype:pre_miRNA
# gene_biotype:SRP_RNA
# gene_biotype:snoRNA
# gene_biotype:antisense_RNA
# gene_biotype:sense_intronic

# 构建configure文件
vim id.config
# gene_biotype:snRNA
# gene_biotype:rRNA
# gene_biotype:tRNA
# gene_biotype:snoRNA
# gene_biotype:SRP_RNA

# 提取上述5种类型的序列
cat id.config | while read id
do
	type=${id#*:}
	# echo $type
	less -S Oryza_sativa.IRGSP-1.0.ncrna.fa | grep ">" | awk '{if($5~/'"$id"'/) print $0}' > Osa_${type}.seqid
done

# 合并上述5种类型的序列
cat Osa_{snRNA,rRNA,tRNA,snoRNA,SRP_RNA}.seqid > ncRNA.seqid
awk '{print $1}' ncRNA.seqid | cut -d ">" -f 2 > clean.ncRNA.seqid
seqkit grep -f clean.ncRNA.seqid Oryza_sativa.IRGSP-1.0.ncrna.fa > 4_Oryza_sativa.IRGSP-1.0.ncrna.fa

# 序列ID重编码
python reformatID.py

# 文件重命名
cp reformat_4_Oryza_sativa.IRGSP-1.0.ncrna.fa Osa_ncRNA.fa
```

构建bowtie索引，

```shell
DT=`date +"%y%m%d%H%M"`
bowtie-build Osa_ncRNA.fa Osa_ncRNA 1>bowtie.buildindex.${DT}.log
```



### 2、构建水稻参考基因组索引

数据来源：`Ensembl Plants`

```shell
# download
wget -c http://ftp.ensemblgenomes.org/pub/plants/release-53/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
gunzip -c Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz > osa.fa
```

构建bowtie索引，

```shell
DT=`date +"%y%m%d%H%M"`
bowtie-build osa.fa osa_genome 1>bowtie.buildindex.${DT}.log
```



### 3、构建水稻hairpin RNA

> Note：hpRNA，即microRNA前体序列。

#### 「稍微拓展一下」如何获取miRBase数据库中，给定物种的缩写（abbreviation）以及microRNA序列？

1）查看miRBase数据库的`organisms.txt.gz`文件，比如水稻在该文件中的metadata格式如下，

```shell
zless -S organisms.txt.gz | grep "Oryza sativa"
# osa     OSA     Oryza sativa    Viridiplantae;Magnoliophyta;monocotyledons;    4530
```



2）获取水稻的hpRNA和mature miRNA序列

使用脚本：`miRBasehandle.py`

> 脚本使用说明：
>
> - 目前阶段该脚本还未添加`sys.argv`，需使用vim/vscode/notepad++等软件自行修改`species_name`参数
> - <font color='yellow'>该脚本自动将U转为T</font>

运行方式为，

```shell
python miRBasehandle.py
```

#### 「切入正题」构建hairpin RNA索引

```shell
DT=`date +"%y%m%d%H%M"`
bowtie-build osa.hairpin.fa hairpin 1>bowtie.buildindex.${DT}.log
```



### 4、构建pos

将原始数据比对到水稻hairpin序列之后，需要根据`pos`文件对其进行过滤。

**目的**：？

pos文件包含信息如下，

- 1st：hpRNA ID
- 2nd：mature microRNA的起始位点
- 3rd：mature microRNA的终止位点
- 4th：mature microRNA ID

> Note：需先使用mature microRNA序列作为query序列，比对到hairpin RNA序列上（数据来自miRBase，可使用`miRBasehandle.py`生成）

#### 1）blast

```shell
# 工作目录：00.miRBase.osa_microRNA
# 构建hairpin RNA blastdatabase
makeblastdb -in osa.hairpin.fa -dbtype nucl -out osa_hairpin -parse_seqids 1>hairpin.makedb.log

# blast short-reads alignment
blastn -query osa.mature.fa -db osa_hairpin -out osa_mature_hairpin_blast.txt -outfmt 6 -task blastn-short -max_hsps 1
```

#### 2）使用`blasthandle.py`生成pos文件

> Note：`blasthandle.py`还未加入`sys.argv`参数，用于参数输入，需自行修改文件，
>
> ```python
> Blastfilter('osa_mature_hairpin_blast.txt') 
> mirnaMatch('Filtered_osa_mature_hairpin_blast.txt')
> Getpos('osa', 'Matched_Filtered_osa_mature_hairpin_blast.txt')
> 
> # osa_mature_hairpin_blast.txt，为上述blast分析结果
> # Blastfilter()，按identity、microRNA长度以及错配碱基数（不可大于0）进行过滤
> # mirnaMatch()，对microRNA和hpRNA的ID进行匹配
> # Getpos()，生成.pos文件
> ```

结果文件：`osa_hairpin_mature.pos`



## （4）Sequence Alignment（序列回帖）

使用软件：`bowtie version 1.0.0`

### 1、Sequence Alignment against non-coding RNA（非编码RNA的比对）

**目的**：过滤比对到non-coding RNA的序列（e.g. 过滤掉串联重复序列所产生的miRNA）

使用bowtie参数如下，

```shell
-p 4     # 线程数为4
-v 0     # 允许错配数为0
-f       # 指示输入序列为fa格式
```

示例代码，

```shell
# 工作目录：02.bowtie
# 构建configure文件
ls ../01.raw_data/*.fa | tr '\t' '\n' > fa.config

# ParaFly
cat fa.config | while read id
do
	DT=`date +"%y%m%d%H%M"`
	rawsample=$(basename $id)
	cleansample=${rawsample%.*}
	# echo $cleansample
	echo "bowtie -p 4 -v 0 -f ../00.index/Osa_ncRNA $id ./${cleansample}.${DT}.sam 2>${cleansample}.${DT}.log" >> Osa_ncrna.bowtie.commandlines
done

# Run ParaFly
ParaFly -c Osa_ncrna.bowtie.commandlines -CPU 3 &
```

### 2、Filtering aligned fasta reads

1）获取需排除的序列ID

```shell
# 工作目录：02.bowtie
# 构建configure文件
ls *.sam > sam_id.config

# 生成ID
cat sam_id.config | while read id
do
	sample=${id%%.*}
	# echo $sample
	awk '{print $1}' $id > ${sample}.IDexcluded.txt
done
```

2）过滤比对到水稻non-coding RNA上的序列

```shell
# 工作目录：03.onefilter_data
# 构建configure文件
ls ../01.raw_data/*.fa > fa.config
ls ../02.bowtie/*.txt > excludedID_text.config
cut -d "/" -f 3 fa.config | cut -d "." -f 1 > sample_name.config
paste sample_name.config fa.config excludedID_text.config > config

# 过滤比对到水稻non-coding RNA上的序列
cat config | while read id 
do
	arr=($id)
	sample=${arr[0]}
	seq=${arr[1]}
	text=${arr[2]}
	# echo $sample $seq $text
	seqkit grep -v -f $text $seq > ${sample}_onefilter.fa
done

# 结果统计
seqkit stats *.fa > onefilter.table.txt
seqkit stats -b ../01.raw_data/*.fa > rawdata.table.txt
```



### 3、Sequence Alignment against Osa reference genome

**目的**：过滤掉测序数据中的污染（e.g. 病毒、细菌等）

```shell
# 工作目录：04.bowtie
# 构建configure文件
ls ../03.onefilter_data/*.fa > fa.config
# ParaFly
cat fa.config | while read id
do
	DT=`date +"%y%m%d%H%M"`
	seq=$(basename $id)
	name=${seq%_*}
	echo "bowtie -p 4 -v 0 -f ../00.osa_genome/osa_genome $id ./${name}.${DT}.sam 2>${name}.${DT}.log" >> osa_genome.bowtie.commandlines
done
# Run ParaFly
ParaFly -c osa_genome.bowtie.commandlines -CPU 3 &
```



### 4、Filtering unaligned reads

1）获取比对到水稻参考基因组上的seq ID

```shell
# 工作目录：04.bowtie
# 构建configure文件
ls *.sam > sam_id.config
# 生成ID
cat sam_id.config | while read id
do
	sample=${id%%.*}
	# echo $sample
	awk '{print $1}' $id > ${sample}.IDincluded.txt
done
```

2）过滤未比对到水稻参考基因组上的序列

```shell
# 工作目录：05.twofilter_data
# 构建configure文件
ls ../03.onefilter_data/*.fa > fa.config
ls ../04.bowtie/*.txt > includedID_text.config
cut -d "/" -f 3 fa.config | cut -d "_" -f 1 > sample_name.config
paste sample_name.config fa.config includedID_text.config > config

# 过滤未比对到水稻参考基因组上的序列
cat config | while read id 
do
	arr=($id)
	sample=${arr[0]}
	seq=${arr[1]}
	text=${arr[2]}
	# echo $sample $seq $text
	echo "seqkit grep -f $text $seq > ${sample}_twofilter.fa" >> seqkit.commandlines
done

# Run ParaFly
ParaFly -c seqkit.commandlines -CPU 3 &

# 结果统计
seqkit stats *.fa > twofilter.table.txt
```



### 5、Sequence Alignment against Osa hairpin sequences

使用bowtie参数，

```shell
-a # 保留所有比对结果
```



```shell
# 工作目录：06.bowtie
# 构建configure文件
ls ../05.twofilter_data/*.fa > fa.config

cat fa.config | while read id
do
	DT=`date +"%y%m%d%H%M"`
	seq=$(basename $id)
	name=${seq%_*}
	echo "bowtie -p 4 -v 0 -a -f /home/changqing/Documents/miRBase/osa_hairpin_u_t_bowtie/hairpin $id ./${name}.${DT}.sam 2>${name}.${DT}.log" >> osa_hairpin.bowtie.commandlines
done

# Run ParaFly
ParaFly -c osa_hairpin.bowtie.commandlines -CPU 3 &
```



### 6、Filtering miRNA based on hairpinRNA - mature miRNA position

1）获取对应水稻microRNA前体ID

```shell
# 工作目录：06.bowtie
# 构建configure文件
ls *.sam > sam_id.config
# 生成ID
cat sam_id.config | while read id
do
	sample=${id%%.*}
	# echo $sample
	awk '{print $3}' $id > ${sample}.microRNAID.txt
done
```

2）根据`.pos`文件进行microRNA过滤以及gene expression matrxi的生成

使用脚本：`pandas_organizebowtie.py`

> Note：`pandas_organizebowtie.py`尚未加入`sys.argv`，需进入脚本手动修改，
>
> ```python
> pos_df = ParsePos('osa_hairpin_mature.pos')
> bowtieoutput_list = BowtieOutputNameGet('bowtie_output')
> # bowtie_output，为06.bowtie文件夹下最终的bowtie比对结果
> MultipleBowtieOuputOrganize(bowtieoutput_list, pos=pos_df)
> ```

结果文件：`pandas_all_samples.csv`



## bowtie比对注意事项

### 1）参数设置

- 同一品种/物种，错配数为0
- 同一物种不同品种允许错配数为1，考虑到SNP的存在

### 2）输出文件格式

默认情况下，输出格式为bwt，即bowtie默认的输出格式。

若要生成`sam`文件格式，需使用`-S`参数。

> Note：本文档中的结果均为bwt格式，但以SAM结尾



# Downstream Analysis

## （1）差异表达分析

分析使用脚本：`DESeq.R`

可视化使用脚本：

- 1）热图：`pheatmap.R`，`DESeq.R`

To be continued



## （2）raw counts to RPM

使用脚本：`RPM.py`

使用数据：

- `pandas_all_samples.csv`
- `04.bowtie`目录下包含的SAM文件



## （3）靶基因预测分析

使用软件：`psRobot`

使用数据：水稻CDS序列，来源为`Ensembl Plants`

```shell
# 数据下载
wget -c http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/oryza_sativa/cds/Oryza_sativa.IRGSP-1.0.cds.all.fa.gz

# 靶基因预测
# 也可使用网页版默认参数即可
psRobot_tar -s input.fa -t target.fa -o target.gTP
```





# Code Availability

https://github.com/Youpu-Chen/microRNome/tree/dev1