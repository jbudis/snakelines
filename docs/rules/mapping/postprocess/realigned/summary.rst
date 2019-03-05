Gatk - Identify Regions To Realign
--------------------------------------

Find places on genomes with indels for realignment.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/realigned/gatk.snake
- *Rule name:* gatk__identify_regions_to_realign

**Input(s):**

- *ref:* Reference genomic sequences in FASTA format
- *bam:* Mapped reads for realignment in BAM format

**Output(s):**

- *intervals:* Regions that need to be realigned

Gatk - Realign Regions
--------------------------

Realign regions locally in places with indels that are problematic for standard mappers.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/realigned/gatk.snake
- *Rule name:* gatk__realign_regions

**Input(s):**

- *ref:* Reference genomic sequences in FASTA format
- *bam:* Mapped reads for realignment in BAM format
- *intervals:* Regions that need to be realigned

**Output(s):**

- *bam:* Realigned reads in BAM format

