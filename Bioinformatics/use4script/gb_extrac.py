#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 18:45:31 2022

@author: xum
"""


#导入模块，初始传递命令
import argparse

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 gb_extrac.py -i [gb] -o [tsv]\n用于提取 genbank 中并转为 tsv')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[gb]', help = '输入文件，gb 文件', required = True)
required.add_argument('-o', '--output', metavar = '[tsv]', help = '输出文件，tsv 文件', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

tsv = open(args.output, "w") #核苷酸gb文件路径
gb = open(args.input,"r")   #输出文件路径
date_str = gb.read()
a = []
a = date_str.split("//")
first_line = "GenID"+"\t"+"Host"+"\t"+"Country"+"\t"+"Collection_date"+"\n"
tsv.write(first_line)
for n in range(len(a)-1):
    b = []
    b = a[n].split("\n")
    seqID = ""
    host = "N/A"
    country = "N/A"
    collection_date = "N/A"
    for m in range(len(b)):
        if b[m].find("LOCUS") != -1:
            seqID = b[m][10:22]
            seqID = seqID.strip(" ")
        if b[m].find("/host") != -1:
            host_str = b[m].strip()
            host_str = host_str.replace('"','')
            host = host_str[6:]
        if b[m].find("/country")!=-1:
            country_str = b[m].strip()
            country_str = country_str.replace('"','')
            country = country_str[9:]
        if b[m].find("/collection_date")!=-1:
            collection_date_str = b[m].strip()
            collection_date_str = collection_date_str.replace('"','')
            collection_date = collection_date_str[17:]

    S = seqID+"\t"+ host+"\t"+country+"\t"+collection_date
    tsv.write(S+"\n")
    
print("Fished")
gb.close()
tsv.close()
