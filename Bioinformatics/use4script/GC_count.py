##python3.6 或以上版本执行
import re

#读取基因组序列，将碱基信息存储到字典中
genome_dict = {}
genome = open('/data/home/xum/MP/circos/p1_2/bam/ref.fa', 'r')
for line in genome:
    line = line.strip()
    if line[0] == '>':
        seq = line.split('>')[1]
        genome_dict[seq] = ''
    else:
        genome_dict[seq] += line    
genome.close()  

 
#读取滑窗位置信息，并统计每段滑窗的 GC 含量并追加在滑窗测序深度统计文件后方
output = open('/data/home/xum/MP/circos/p1_2/bam/0825-2_depth_gc.txt', 'w')
print('seq\tstart\tend\tDepth\tGC\tA\tC\tG\tT', file = output)

depth = open("/data/home/xum/MP/circos/p1_2/bam/0825-2_depth.txt", 'r')
for l in depth:
    l= l.strip().split()
    GC = len(re.findall('[GCgc]', genome_dict[l[0]][int(l[1]):int(l[2])]))/(int(l[2])-int(l[1]))
    A = len(re.findall('[Aa]', genome_dict[l[0]][int(l[1]):int(l[2])]))/(int(l[2])-int(l[1]))
    C = len(re.findall('[Cc]', genome_dict[l[0]][int(l[1]):int(l[2])]))/(int(l[2])-int(l[1]))
    G = len(re.findall('[Gg]', genome_dict[l[0]][int(l[1]):int(l[2])]))/(int(l[2])-int(l[1]))
    T = len(re.findall('[Tt]', genome_dict[l[0]][int(l[1]):int(l[2])]))/(int(l[2])-int(l[1]))
    print(f'{l[0]}\t{int(l[1])+1}\t{l[2]}\t{l[3]}\t{GC}\t{A}\t{C}\t{G}\t{T}', file = output)
    
depth.close()
output.close()

