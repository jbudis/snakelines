Samtools - Filter Fragments
-------------------------------

 Filter reads that fall into selected regions

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/filtered/samtools.snake
- *Rule name:* samtools__filter_fragments

**Input(s):**

- *bam:* Mapped reads in BAM format
- *bed:* Regions to filter in BED format

**Output(s):**

- *bam:* Filtered reads in BAM format

**Param(s):**

- *min_map_quality:* Minimal quality of mapping
- *regions:* Name of the directory with bed file under reference/{reference}/annotation/{regions}/regions.bed

Bamtools - Filter Fragments
-------------------------------

 Remove fragments that do not meet user specified conditions. All values are inclusive (i.e. <= or =>)

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/filtered/bamtools.snake
- *Rule name:* bamtools__filter_fragments

**Input(s):**

- *bam:* Mapped reads in BAM format

**Output(s):**

- *bam:* Filtered reads in BAM format

**Param(s):**

- *min_map_quality:* Minimal quality of mapping
- *drop_improper_pairs:* Eliminate reads that do not pass paired-end resolution

