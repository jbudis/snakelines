Samtools - Sort Mapped Reads
--------------------------------

Sort aligned reads according to mapped position on reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/mapping/postprocess/sorted/samtools.snake
- *Rule name:* samtools__sort_mapped_reads

**Input(s):**

- *ref:* Reference genomic sequences in fasta format
- *bam:* Unordered mapped reads in bam format

**Output(s):**

- *bam:* Ordered mapped reads according to their location on reference genome

