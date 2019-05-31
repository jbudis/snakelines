Custom - Infer Read Groups
------------------------------

Infer sample name, flow cell, barcode and lanes from the input BAM file from Illumina sequencing.
Write meta information as read groups to the output bam file.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/read_group/custom.snake
- *Rule name:* custom__infer_read_groups

**Input(s):**

- *bam:* Mapped reads in BAM format

**Output(s):**

- *bam:* Mapped reads with inferred read groups in BAM format

**Param(s):**

- *infer_script:* File path to script that do real job here

