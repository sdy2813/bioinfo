#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: xum
"""

from Bio import SeqIO
#导入模块，初始传递命令
import argparse

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 gbk2faa.py -i [gbk] -o [faa]\n用于提取 genbank 中并转为 faa')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[gbk]', help = '输入文件，gbk 文件', required = True)
required.add_argument('-o', '--output', metavar = '[faa]', help = '输出文件，faa 文件', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#gbk_filename = "0803-1.gbk"
#faa_filename = "0803-1.faa"
input_handle  = open(args.input, "r")
output_handle = open(args.output, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                                     seq_feature.qualifiers['locus_tag'][0],
                                     seq_record.name,
                                     seq_feature.qualifiers['translation'][0]))
print("Fished")
output_handle.close()
input_handle.close()
