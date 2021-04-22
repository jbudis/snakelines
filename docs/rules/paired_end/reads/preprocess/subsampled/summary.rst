Seqtk - Subsample Reads
---------------------------

Randomly select user-configured number of reads from fastq files

**Location**

- *Filepath:* <SnakeLines_dir>/rules/paired_end/reads/preprocess/subsampled/seqtk.snake
- *Rule name:* seqtk__subsample_reads

**Input(s):**

- *r1:* Left side of sequenced fragments in gzipped fastq format
- *r2:* Right side of sequenced fragments in gzipped fastq format

**Output(s):**

- *r1:* Left side of subsampled fragments without overlap in gzipped fastq format
- *r2:* Right side of subsampled fragments without overlap in gzipped fastq format

