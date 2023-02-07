"""
Created on 2023/2/7
@author: chatgpt
"""

from Bio import SeqIO

# 读取大的GenBank文件
records = SeqIO.parse("big_file.gbk", "genbank")

# 遍历每一个记录
for i, record in enumerate(records):
    # 将每个记录写入一个独立的文件
    SeqIO.write(record, "record_{}.gbk".format(i), "genbank")
