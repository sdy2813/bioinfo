#!/bin/bash

# 定义帮助函数
helpFunction()
{
   echo ""
   echo "Usage: bam2cns.sh -i input.bam -o output.fasta -d depth -f reference.fa"
   echo -e "\t-i Specify the input BAM file."
   echo -e "\t-o Specify the output FASTA file."
   echo -e "\t-d Specify the depth for consensus calling."
   echo -e "\t-f Specify the reference genome file in fasta format."
   exit 1 # Exit script after printing help
}

# 解析命令行参数
while getopts "i:o:d:f:h" flag
do
    case "${flag}" in
        i) bamfile=${OPTARG};;
        o) fastafile=${OPTARG};;
        d) depth=${OPTARG};;
        f) referencefile=${OPTARG};;
        h) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# 检查是否已经设置深度和参考基因组文件
if [ -z "$depth" ] || [ -z "$referencefile" ]
then
   echo "Some or all of the parameters are empty.";
   helpFunction
fi

# 输出参数
echo "Input File: $bamfile";
echo "Output File: $fastafile";
echo "Depth: $depth";
echo "Reference genome: $referencefile";

# 获取一致性序列
samtools mpileup -Q 0 -d $depth -uf $referencefile $bamfile | bcftools call -c - | vcfutils.pl vcf2fq -d $depth | seqtk seq -aQ64 -q 20 -n N > $fastafile
