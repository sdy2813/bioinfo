# spades 基本操作

## 软件安装
```bash
mamba search spades

mamba install -c bioconda spades
```

## 基本使用
```
# PE150
spades.py -k 21,33,55,77 --careful -1 Sample04_R1.fastq.gz -2 Sample04_R2.fastq.gz -o ../assemble/sample04/spades_output -t 12

# PE250
spades.py -k 21,33,55,77,99,127 --careful -1 Sample04_R1.fastq.gz -2 Sample04_R2.fastq.gz -o ../assemble/sample04/spades_output -t 12
```







