# fastqc基本操作

1.质控
```bash
# 基本格式

# fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN

# 主要是包括前面的各种选项和最后面的可以加入N个文件
# -o --outdir FastQC生成的报告文件的储存路径，生成的报告的文件名是根据输入来定的
# --extract 生成的报告默认会打包成1个压缩文件，使用这个参数是让程序不打包
# -t --threads 选择程序运行的线程数，每个线程会占用250MB内存，越多越快咯
# -c --contaminants 污染物选项，输入的是一个文件，格式是Name [Tab] Sequence，里面是可能的污染序列，如果有这个选项，FastQC会在计算时候评估污染的情况，并在统计的时候进行分析，一般用不到
# -a --adapters 也是输入一个文件，文件的格式Name [Tab] Sequence，储存的是测序的adpater序列信息，如果不输入，目前版本的FastQC就按照通用引物来评估序列时候有adapter的残留
# -q --quiet 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况。
```

2.fastqc结果查看
1. 产生两个结果文件：
+ html：网页版结果
+ zip：本地结果压缩文件


网页版结果
网页版结果页面左上角是一个summary,就是整个报告的目录，整个报告分成若干个部分。合格会有个绿色的对勾，警告是个“!”，不合格是个红色的叉子。  
![summary](https://upload-images.jianshu.io/upload_images/7493830-979e62961335827f?imageMogr2/auto-orient/strip|imageView2/2/w/274/format/webp)  

需要重点关注的结果：
+ Basic Statistics：**对数据量的概览**
	+ ![summary](https://pic2.zhimg.com/b19f28a107a90f7c4f7b02510c7884ab_r.jpg)  
	+ Encoding指测序平台的版本和相应的编码版本号，这个在计算Phred反推error P的时候有用，如果不明白可以参考之前的文章。
	+ Total Sequences记录了输入文本的reads的数量
	+ Sequence length 是测序的长度
	+ %GC 是我们需要重点关注的一个指标，这个值表示的是整体序列中的GC含量，这个数值一般是物种特意的，比如人类细胞就是42%左右。
+ Per base sequence quality：reads每个位置测序质量最直接的展示
	+ ![Per base sequence quality](https://picb.zhimg.com/38670ee6d5f373e326e3fe3e23ba4f9b_r.jpg)  
	+ 此图中的横轴是测序序列第1个碱基到第101个碱基
	+ 纵轴是质量得分，Q = -10*log10（error P）即20表示1%的错误率，30表示0.1%
	+ 图中每1个boxplot，都是该位置的所有序列的测序质量的一个统计，上面的bar是90%分位数，下面的bar是10%分位数，箱子的中间的横线是50%分位数，箱子的上边是75%分位数，下边是25%分位数
	+ 图中蓝色的细线是各个位置的平均值的连线
	+ 一般要求此图中，所有位置的10%分位数大于20,也就是我们常说的Q20过滤
	+ 所以上面的这个测序结果，需要把后面的87bp以后的序列切除，从而保证后续分析的正确性
	+ Warning 报警 如果任何碱基质量低于10,或者是任何中位数低于25
	+ Failure 报错 如果任何碱基质量低于5,或者是任何中位数低于20
+ Per sequence quality scores：总体reads测序质量趋势
+ Per base sequence content：ATGC含量估计测序是否存在偏差
+ Sequence Duplication Levels]：影响测序的因素太多，查看是否存在污染，数据处理时是否需要去冗余；现在数据量都可以满足需求，因此前期数据处理时，尽量高标准，严格质控；。












	
