Samtools - Convert Sam To Bam
---------------------------------

Convert .bam to .sam.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/mapping/index/samtools.snake
- *Rule name:* samtools__convert_sam_to_bam

**Input(s):**

- *bam:* Mapped reads in bam format

**Output(s):**

- *sam:* Mapped reads in sam format

Samtools - Bam Index
------------------------

Generate .bai index to .bam files to quick recover reads from genomic location of interest.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/mapping/index/samtools.snake
- *Rule name:* samtools__bam_index

**Input(s):**

- *bam:* Mapped reads in bam format

**Output(s):**

- *bai:* Index to mapped reads for enable fast read retrieval from desired genomic region

