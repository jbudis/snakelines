Bamtools - Merge Bams Into Single Bam File
----------------------------------------------

Merge mapped reads in BAM files into a single BAM file

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/merged/bamtools.snake
- *Rule name:* bamtools__merge_bams_into_single_bam_file

**Input(s):**

- *bams:* BAM files with mapped reads

**Output(s):**

- *bam:* single BAM file with all reads from input BAM files
- *list:* List of files that were used for merging (one line per file)

