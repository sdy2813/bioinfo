# fastp

```
for i in {1..6}
do
fastp -i $read1 I20210910-${i}_L5_1.fq.gz -o ./NO_adapter/I10-${i}_R1.fq.gz -I I20210910-${i}_L5_2.fq.gz -O ./NO_adapter/I10-${i}_R2.fq.gz -z 9 -q 20 -l 20 -w 16 -c -h ./NO_adapter/I10-${i}.html 2>./NO_adapter/I10-${i}_log.txt
done
```
