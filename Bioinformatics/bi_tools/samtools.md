# samtools基本操作


## 软件安装
```bash
conda install samtools
```
## samtools的具体使用
```bash
   >>$ samtools

Program: samtools (Tools for alignments in the SAM format)
Version: 1.11 (using htslib 1.11)

Usage:   samtools <command> [options]

Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     fqidx          index/extract FASTQ
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates
     ampliconclip   clip oligos from the end of reads

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA

  -- Statistics
     bedcov         read depth per BED region
     coverage       alignment depth and percent coverage
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)
     ampliconstats  generate amplicon specific stats

  -- Viewing
     flags          explain BAM flags
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM


```

### samtools view
```bash
   >>$ samtools view

Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

Options:
  -b       output BAM 
  # 默认下输出是 SAM 格式文件，该参数设置输出 BAM 格式
  -C       output CRAM (requires -T)
  -1       use fast BAM compression (implies -b)
  -u       uncompressed BAM output (implies -b)
  -h       include header in SAM output
  -H       print SAM header only (no alignments)
  -c       print only the count of matching records
  -o FILE  output file name [stdout]
  -U FILE  output reads not selected by filters to FILE [null]
  -t FILE  FILE listing reference names and lengths (see long help) [null]
  -X       include customized index file
  -L FILE  only include reads overlapping this BED FILE [null]
  -r STR   only include reads in read group STR [null]
  -R FILE  only include reads with read group listed in FILE [null]
  -d STR:STR
           only include reads with tag STR and associated value STR [null]
  -D STR:FILE
           only include reads with tag STR and associated values listed in
           FILE [null]
  -q INT   only include reads with mapping quality >= INT [0]
  -l STR   only include reads in library STR [null]
  -m INT   only include reads with number of CIGAR operations consuming
           query sequence >= INT [0]
  -f INT   only include reads with all  of the FLAGs in INT present [0]
  -F INT   only include reads with none of the FLAGS in INT present [0]
  -G INT   only EXCLUDE reads with all  of the FLAGs in INT present [0]
  -s FLOAT subsample reads (given INT.FRAC option value, 0.FRAC is the
           fraction of templates/read pairs to keep; INT part sets seed)
  -M       use the multi-region iterator (increases the speed, removes
           duplicates and outputs the reads as they are ordered in the file)
  -x STR   read tag to strip (repeatable) [null]
  -B       collapse the backward CIGAR operation
  -?       print long help, including note about region specification
  -S       ignored (input format is auto-detected)
  --no-PG  do not add a PG line
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
  -T, --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
      --write-index
               Automatically index the output files [off]
      --verbosity INT
               Set level of verbosity
```

常用命令
```bash
# 使用samtools查看sam文件的header部分
samtools view -H test.sam

# 将sam文件转换成bam文件
samtools view -bS abc.sam > abc.bam
# 或
samtools view -b -S abc.sam -o abc.bam

# 提取比对到参考序列上的比对结果
samtools view -bF 4 abc.bam > abc.F.bam

# 提取paired reads中两条reads都比对到参考序列上的比对结果，只需要把两个4+8的值12作为过滤参数即可
samtools view -bF 12 abc.bam > abc.F12.bam

# 提取没有比对到参考序列上的比对结果
samtools view -bf 4 abc.bam > abc.f.bam

# 提取bam文件中比对到caffold1上的比对结果，并保存到sam文件格式
samtools view abc.bam scaffold1 > scaffold1.sam

# 提取scaffold1上能比对到30k到100k区域的比对结果
samtools view abc.bam scaffold1:30000-100000 $gt; scaffold1_30k-100k.sam

# 根据fasta文件，将 header 加入到 sam 或 bam 文件中
samtools view -T genome.fasta -h scaffold1.sam > scaffold1.h.sam
```


### samtools sort

```bash
   >>$ samtools sort
Usage: samtools sort [options...] [in.bam]
Options:
  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
  -u         Output uncompressed data (equivalent to -l 0)
  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
  -M         Use minimiser for clustering unaligned/unplaced reads
  -K INT     Kmer size to use for minimiser [20]
  -n         Sort by read name (not compatible with samtools index command)
  -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
  -o FILE    Write final output to FILE rather than standard output
  -T PREFIX  Write temporary files to PREFIX.nnnn.bam
  --no-PG    do not add a PG line
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
      --verbosity INT
               Set level of verbosity

```


实例
```bash
#  tmp.bam 按照序列名称排序，并将结果输出到tmp.sort.bam  
samtools sort -n tmp.bam  -o tmp.sorted.bam   
```

### samtools index
```bash
   >>$ samtools index
Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
Options:
  -b       Generate BAI-format index for BAM files [default]
  -c       Generate CSI-format index for BAM files
  -m INT   Set minimum interval size for CSI indices to 2^INT [14]
  -@ INT   Sets the number of threads [none]
```

对排序后的序列建立索引，并输出为bai文件，用于快速随机处理。在很多情况下，特别是需要显示比对序列的时候，bai文件是必不可少的，例如之后的tview命令。
```
samtools index abc.sort.bam
```


### samtools depth
```
   >>$ samtools depth

Usage: samtools depth [options] in1.bam [in2.bam [...]]
Options:
   -a                  output all positions (including zero depth)
   -a -a (or -aa)      output absolutely all positions, including unused ref. sequences
   -b <bed>            list of positions or regions
   -X                  use customized index files
   -f <list>           list of input BAM filenames, one per line [null]
   -H                  print a file header
   -l <int>            read length threshold (ignore reads shorter than <int>) [0]
   -d/-m <int>         maximum coverage depth [8000]. If 0, depth is set to the maximum
                       integer value, effectively removing any depth limit.
   -o FILE             where to write output to [stdout]
   -q <int>            base quality threshold [0]
   -Q <int>            mapping quality threshold [0]
   -r <chr:from-to>    region
   -g <flags>          remove the specified flags from the set used to filter out reads
   -G <flags>          add the specified flags to the set used to filter out reads
                       The default set is UNMAP,SECONDARY,QCFAIL,DUP or 0x704
   -J                  include reads with deletions in depth computation
   -s                  for the overlapping section of a read pair, count only the bases
                       of a single read. This option requires raising the base quality
                       threshold to 1.
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
      --verbosity INT
               Set level of verbosity

The output is a simple tab-separated table with three columns: reference name,
position, and coverage depth.  Note that positions with zero coverage may be
omitted by default; see the -a option.

```

得到每个碱基位点的测序深度，并输出到标准输出。输入的bam文件必须先做samtools index
```
samtools depth tmp.index.bam  >  tmp.depth.bam
```

[[SAMtools] 常用指令总结](https://www.cnblogs.com/xiaofeiIDO/p/6805373.html)

[samtools常用命令详解](https://www.cnblogs.com/emanlee/p/4316581.html)

[SAM/BAM相关的进阶知识](https://ming-lian.github.io/2019/02/07/Advanced-knowledge-of-SAM/)

[宏基因组/转录组去除宿主污染](http://blog.sciencenet.cn/blog-2379401-1268509.html)




