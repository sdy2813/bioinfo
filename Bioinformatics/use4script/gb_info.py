#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:37:13 2023

@author: xum
"""

from Bio import SeqIO

# 读取GenBank文件
records = SeqIO.parse("KX981911.gb", "genbank")

# 打开输出文件
with open("gb.info.tsv", "w") as f:
    # 写入标题行
    f.write("Accession\tLocation\tCollection Date\tHost\n")
    
    # 遍历所有记录
    for record in records:
        accession = record.id
        for feature in record.features:
            if feature.type == "source":
                # 提取地区信息
                location = feature.qualifiers.get("country", [None])[0]
                # 提取采样日期信息
                collection_date = feature.qualifiers.get("collection_date", [None])[0]
                # 提取host信息
                host= feature.qualifiers.get("host", [None])[0]
                # 写入记录
                f.write("{}\t{}\t{}\t{}\n".format(accession, location, collection_date, host))
