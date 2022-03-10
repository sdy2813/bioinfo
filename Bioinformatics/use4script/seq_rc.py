#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 19:52:58 2022

@author: xum
"""

from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 seq_rc.py -i [fa] -o [fa]\n用于获取反向序列和反向互补序列')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i'， '--input'， metavar = '[fa]'， help = '输入文件，fa 文件'， required = True)
required.add_argument('-o'， '--output'， metavar = '[fa]'， help = '输出文件，fa 文件'， required = True)
optional.add_argument('-c'， '--complement'， help = '互补序列'， action = "store_true")
optional.add_argument('-rc'， '--reverse_complement'， help = '反向互补序列'， action = "store_true")
optional.add_argument('-h'， '--help', action = 'help'， help = '帮助信息')

args = parser.parse_args()
seq=open(args.input,"r")
out=open(args.output, "w")

a=[]
for line in seq:
    line=line.strip()
    a.append(line)
for i in range (len(a)):
    # print(Seq(a[1]).reverse_complement())
    if i % 2==0:
        out.write(a[i]+"\n")
    if i % 2==1: 
        if args.reverse_complement:
            a[i]=Seq(a[i]).reverse_complement()
            out.write(str(a[i]) + "\n")
        if args.complement:
            a[i]=Seq(a[i]).complement()
            out.write(str(a[i]) + "\n")



