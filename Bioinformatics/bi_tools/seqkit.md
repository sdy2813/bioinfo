# <center>seqkit基本操作</center>

## 序列操作
```
# 取反向序列
seqkit seq *.fa -r > out 
# 取互补序列
seqkit seq *.fa -p > out 
# 取反向互补序列
seqkit seq *.fa -r -p > out 
# DNA序列转换为RNA序列
seqkit seq *.fa --nda2rna > out
# RNA序列转换为DNA序列
seqkit seq *.fa rna2dna  > out 
# 将序列以小写字母的形式输出
seqkit seq *.fa -l > out 
# 将序列以大写字母的形式输出
seqkit seq *.fa -u > out 
# (指定序列的长度为10)指定每行序列的输出长度(为0的话，代表为一整行，默认的输出 长度是60个碱基)
seqkit seq *.fa -w 10 > out 
 # 将多行序列转换为一行序列
seqkit seq *.fa -w 0 > out
# 只输出序列
seqkit seq *.fa -s -w 0 > out 
# 将只输出的序列的，指定每行输出的碱基数
seqkit seq *.fa -s -w 40 > out 
```

## 文件转换
```
# 将fataq文件转化为fasta格式.
seqkit seq fq2fa test.fq -o test.fa

# 将fasta格式转化为tab格式
seqkit fx2tab test.fa > test.tab

# 将tab格式转为fasta格式
seqkit tab2fx test.tab > test.fa
```

## 序列信息统计
```
# 序列碱基含量
seqkit fx2tab -l -g -n -i -H *fa 
# 序列长度的整体分布统计
seqkit stat *.fa 
```
## 根据ID或特定的motif筛选提取序列
```
# 选取有起始密码子的序列
seqkit grep -s -r -i -p ^atg cds.fa
# 根据ID提取序列
seqkit grep -f list test.fa > new.fa
# 简并碱基使用。S 代表C or G.
seqkit grep -s -d -i -p TTSAA
# 匹配限定到某区域
seqkit grep -s -R 1:30 -i -r -p GCTGG
```

## 多个序列文件比较寻找相同的序列或者ID相同的序列
```
# By ID (default,>后面，空格之前的名字)输出ID名字相同的。
seqkit common test1.fa test2.fa -o common.fasta 
# By full name（整个序列的名字，包含description部分）。输出序列名字相同的。
seqkit common test1.fa test2.fa -n -o common.fasta 
# 输出要比较的文件中序列相同的序列
seqkit common test1.fa test2.fa -s -i -o common.fasta 
# 输出要比较的文件中序列相同的序列 (for large sequences)
seqkit common test1.fa test2.fa -s -i -o common.fasta --md5 
```

## 提取序列
```
# 随机抽取序列
seqkit sample -n 10000  -s 11  test1_1.fq -o sample.fq
seqkit sample -p 0.1 -s 11  test1_1.fq -o sample.fq

#  排序输出命令
seqkit sort -l test.fa

#文件切割
seqkit split hairpin.fa.gz -p 4
```
