Trimmomatic - Trim Reads
----------------------------

Remove low quality ends of reads and then filter reads that are too short.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/paired_end/reads/preprocess/trimmed/trimmomatic.snake
- *Rule name:* trimmomatic__trim_reads

**Input(s):**

- *r1:* Left side of sequenced fragments in gzipped fastq format
- *r2:* Right side of sequenced fragments in gzipped fastq format

**Output(s):**

- *r1:* Left side of subsampled fragments without overlap in gzipped fastq format
- *r2:* Right side of subsampled fragments without overlap in gzipped fastq format

