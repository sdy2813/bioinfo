# fastuniq基本操作

FastUniq 是一种用于移除重复的测序 reads 的工具，特别是对于成对端测序数据。以下是其基本使用方法：

首先，你需要创建一个输入文件列表（例如，input_list.txt），在这个文件中，每行包含一对配对的 FASTQ 文件。文件应该是这样的：

```
# input_list.txt
reads_1.fastq
reads_2.fastq
```
```
# 运行 FastUniq
# `-i` 参数指定输入文件列表。
# `-t q` 参数指定输入文件的格式是 FASTQ 格式。
# `-o` 参数指定输出的第一对 FASTQ 文件。
# `-p` 参数指定输出的第二对 FASTQ 文件。
fastuniq -i input_list.txt -t q -o out_reads_1.fastq -p out_reads_2.fastq
```


FastUniq 主要是用于移除 PCR 重复，因此通常在测序 reads 的质控步骤中使用。但请注意，如果你的研究目标是定量的（例如，RNA-seq 或 ChIP-seq），移除 PCR 重复可能会影响后续的定量分析。在这些情况下，你可能需要采用不同的策略来处理 PCR 重复。


---
参考资料：

1. [去除PCR冗余](http://www.excel-jiqiao.com/subject/cwagbftx.html)

1. [FastUniq去除paired reads的duplicates](http://blog.sina.com.cn/s/blog_670445240101lqat.html)

1. [FastUniq去除paired reads的duplicates](https://qinqianshan.com/bioinformatics/trim/fastuniq/)
