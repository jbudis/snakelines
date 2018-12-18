Seqtk - Subsample Reads
---------------------------

Randomly select user-configured number of reads from fastq files.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/preprocess/subsampled/seqtk.snake
- *Rule name:* seqtk__subsample_reads

**Input(s):**

- *r1:* filename with R1 reads
- *r2:* filename with R2 reads

**Output(s):**

- *r1:* filename with subsampled R1 reads
- *r2:* filename with subsampled R2 reads

**Param(s):**

- *seed:* int - seed of the random generator for subsampling
- *n_reads:* int - number of reads to keep in subsampled set

