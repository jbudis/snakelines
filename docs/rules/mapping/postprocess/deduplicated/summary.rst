Picard - Mark Duplicates
----------------------------

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as
originating from a single fragment of DNA.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/deduplicated/picard.snake
- *Rule name:* picard__mark_duplicates

**Input(s):**

- *bam:* Mapped reads in bam format
- *bai:* Index to mapped reads for enable fast read retrieval from desired genomic region

**Output(s):**

- *bam:* Mapped reads with additional attribute that marks if read is PCR duplicate

