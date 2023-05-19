# <p align='center'>hisat2基本操作</p>

Hisat2 是一种用于 RNA-seq 数据分析的软件，主要用于将 RNA-seq reads 映射到参考基因组。以下是其基本使用方法：


## 软件安装
```bash
mamba search hisat2

mamba install -c bioconda hisat2
```


## 基本使用

```bash
# 第一步：构建索引
# 使用 HISAT2 构建参考基因组的索引。
# 输入文件为参考基因组的 fasta 文件，例如 `genome.fa`。
# 输出文件将是索引文件，例如 `genome_index`。
hisat2-build genome.fa genome_index

# 第二步：对reads进行比对
# 使用 HISAT2 对 RNA-seq reads 进行比对。
# 输入文件为两个 fastq 文件，例如 `reads_1.fastq` 和 `reads_2.fastq`。
# 输出文件将是比对结果的 SAM 文件，例如 `aligned.sam`。
hisat2 -x genome.index -1 R_1.fq.gz -2 R_2.fq.gz | samtools sort -O bam -@ 14 -o - > *.bam 

```











