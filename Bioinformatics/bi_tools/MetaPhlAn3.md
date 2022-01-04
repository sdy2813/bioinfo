# MetaPhlAn3基本操作

## 安装
```bash
# 安装软件
mamba search MetaPhlAn
mamba create -n MetaPhlAn MetaPhlAn=3.0.13

# 安装数据库
metaphlan --install --bowtie2db ~/db/metaphlan
```


```bash
# 基本使用
# --bowtie2db 指定数据库路径，，就是你上一步搁数据库的地方
# --bowtie2out 输出比对结果
# --input_type 指定输入文件格式，第一次都是fastq，后面可以选择bowtie2out以节约比对的时间
 metaphlan I10-6_R1.fq.gz,I10-6_R2.fq.gz --bowtie2db ~/db/metaphlan/ --bowtie2out ~/MP/mpa/I10-6.bowtie2.bz2 --nproc 15 --input_type fastq -o ~/MP/mpa/I10-6.mpa

```
