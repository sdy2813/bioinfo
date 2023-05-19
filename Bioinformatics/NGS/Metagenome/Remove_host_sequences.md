# Computational removal of host sequences (paired-end reads)

## 1. Using bowtie2 option:  `--un-conc`
quick solution to get the paired reads that do not map to the host reference genome (both reads unmapped).

### download ready to use bowtie2 database of human host genome GRCh38 (hg38)
```
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

unzip GRCh38_noalt_as.zip
```
### run bowtie2 mapping (using --un-conc-gz to get gzip compressed output files;  8 processors)
```
bowtie2 -p 8 -x GRCh38_noalt_as \

  -1 SAMPLE_R1.fastq.gz \

  -2 SAMPLE_R2.fastq.gz \

  --un-conc-gz \

  SAMPLE_host_removed \

  > SAMPLE_mapped_and_unmapped.sam
```


### bowtie2 results (gz files without gz ending)
```
ls

  SAMPLE_host_removed.1

  SAMPLE_host_removed.2

 rename host-sequence free samples

mv SAMPLE_host_removed.1 SAMPLE_host_removed_R1.fastq.gz

mv SAMPLE_host_removed.2 SAMPLE_host_removed_R2.fastq.gz
```

Option `--un-conc` shows results like samtools options -F 2 (excluding reads "mapped in proper pair").

Paired reads that do not map both to the host sequence might still be included in the "host removed" output.

For better control about read filtering options, see workflow below.

If multi-processor option -p is used, output reads might have a different order compared to input files.

Use option --reorder to keep the original read order.

(read order refers to .sam output but might effect also host-removed read output files .1 .2)

http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#performance-options


## 2. Using bowtie2 together with samtools and bedtools

complex solution that gives better control over the rejected reads by using SAM-flags

How to filter out host reads from paired-end fastq files?

2.1 bowtie2 mapping against host genome: write all (mapped and unmapped) reads to a single .bam file

2.2 samtools view: use filter-flags to extract unmapped reads

2.3 samtools fastq: split paired-end reads into separated R1 and R2 fastq files

### 2.1 bowtie2 mapping against host sequence
``` 
# create or download ready to use bowtie2 index databases (host_DB)

# create bowtie2 database using a reference genome (fasta file)

  bowtie2-build host_genome_sequence.fasta host_DB

       → Download reference genomes of common hosts (human, mouse, ..)

       → How to create a bowtie2 index database

   # or, download ready to use bowtie2 database of human genome GRCh38 (hg38)

  wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

  unzip GRCh38_noalt_as.zip

       # move all files into your working directory (or into your predefined $BOWTIE2_INDEXES location)

       # use  "GRCh38_noalt_as" for your bowtie2 database name (instead of "host_DB")

        → see all available bowtie2 databases  (host species list is shown on the right)

# bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
  bowtie2 -p 8 -x host_DB -1 SAMPLE_R1.fastq.gz -2 SAMPLE_R2.fastq.gz | samtools sort -O bam -@ 10 -o - > 4.bam
```


### 2.2 filter required unmapped reads
```
# SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)

 samtools view -b -f 12 -F 256 \

   SAMPLE_mapped_and_unmapped.bam \

   > SAMPLE_bothReadsUnmapped.bam 

    -f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>

    -F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>
```
see meaning of [SAM-flags](http://broadinstitute.github.io/picard/explain-flags.html)

### 2.3 split paired-end reads into separated fastq files .._R1 .._R2
```
# sort bam file by read name ( -n ) to have paired reads next to each other (2 parallel threads, each using up to 5G memory)

samtools sort -n -m 5G -@ 2 SAMPLE_bothReadsUnmapped.bam -o SAMPLE_bothReadsUnmapped_sorted.bam

samtools fastq -@ 8 SAMPLE_bothReadsUnmapped_sorted.bam \

  -1 SAMPLE_host_removed_R1.fastq.gz \

  -2 SAMPLE_host_removed_R2.fastq.gz \

  -0 /dev/null -s /dev/null -n
  
 或者
 
 bedtools bamtofastq -i SAMPLE_bothReadsUnmapped_sorted.bam \
	-fq sample_remove_host_1.fastq \
	-fq2 sample_remove_host_2.fastq
```

http://www.htslib.org/doc/samtools-fasta.html

see other options for  → converting .bam to .fastq

Result

Two files of paired-end reads, containing non-host sequences

SAMPLE_host_removed_R1.fastq.gz

SAMPLE_host_removed_R2.fastq.gz



---
1. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#output-options)






