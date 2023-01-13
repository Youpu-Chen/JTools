# Shell

- `run_RNAseq-full.sh`：shell流程。完整度较snakemake略差。



# snakemake

！！！流程搭建只完成了95%，剩下5%为流程分析步骤之间的耦合性以及高级性能的加入。

- `snakemake.dir.treeview`：目录文件，即输入文件需要如何存放等
- `Snakefile`：snakemake流程文件

```shell
# snakemake-demo
snakemake -c16 --use-conda 03.featureCounts/wtssRNA_seq.mat     # 得到最终的表达矩阵
snakemake -c16 --use-conda 01.qc_initial/wtssRNA_seq/wtssRNA_seq_fastqc.{zip,html}    # 得到raw data的质控结果
snakemake -c16 --use-conda 01.qc_after/wtssRNA_seq/wtssRNA_seq_fastqc.{zip,html}     # 得到过滤后数据的质控结果
snakemake -c16 --use-conda 02.aligments/wtssRNA_seq.sorted.rev.sam 02.aligments/wtssRNA_seq.sorted.fwd.sam 03.featureCounts/wtssRNA_seq.mat 00.genome_index/genome.fa.fai 00.genome_index/genome.txt 00.genome_index/genome.bed   # 流程的完成结果

# single-run
snakemake -c16 --use-conda 01.clean_data/wtssRNA_seq_trimmed.fq 01.clean_data/wtssRNA_seq.fastq_trimming_report.txt
snakemake -c16 --use-conda 02.aligments/wtssRNA_seq.sam 
snakemake -c16 --use-conda 03.featureCounts/wtssRNA_seq.mat
snakemake -c16 --use-conda 04.variants/wtssRNA_seq.vcf
```

