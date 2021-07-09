Bowtie2 - Filter Reads From Reference
-----------------------------------------

Remove reads that do not map to the reference, and so may be caused by contamination in lab processing.
Alternatively, using keep: False configuration removes all fragments that belongs to reference, and so is suitable
to remove contamination caused by host with known genome, e.g. human fragments.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/single_end/reads/preprocess/decontaminated/bowtie2.snake
- *Rule name:* bowtie2__filter_reads_from_reference

**Input(s):**

- *reads:* Left side of sequenced fragments in gzipped fastq format
- *index_files:* List of reference bowtie2 databases

**Output(s):**

- *reads:* Left side of filtered fragments in gzipped fastq format

