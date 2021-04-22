Custom - Remove Overlapping Reads
-------------------------------------

Removes overlapping mapped reads, keeping only a single copy.
This ensures that the each base is covered with at most one read.
The tool is suitable for specific use cases.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/mapping/postprocess/deoverlapped/custom.snake
- *Rule name:* custom__remove_overlapping_reads

**Input(s):**

- *bam:* Mapped reads in bam format
- *bai:* Index to mapped reads for enable fast read retrieval from desired genomic region

**Output(s):**

- *bam:* Mapped reads without overlapping reads.

