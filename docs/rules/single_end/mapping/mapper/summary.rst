Minimap2 - Map Reads To Reference
-------------------------------------

Align reads to the reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/single_end/mapping/mapper/minimap2.snake
- *Rule name:* minimap2__map_reads_to_reference

**Input(s):**

- *reads:* gzipped fastq file with reads, e.g. 'reads/%s/{sample}.fastq.gz'
- *ref:* Reference genomic sequences in fasta format

**Output(s):**

- *sam:* Ordered mapped reads according to their location on reference genome

Bowtie2 - Map Reads To Reference
------------------------------------

For input preprocessed reads bowtie2 finds the most similar genomic region in the provided reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/single_end/mapping/mapper/bowtie2.snake
- *Rule name:* bowtie2__map_reads_to_reference

**Input(s):**

- *read:* gzipped fastq file with reads, e.g. 'reads/%s/{sample}_R1.fastq.gz'
- *index:* reference index (created by rule bowtie2__prepare_index), e.g. 'reference/{reference}/bowtie2_index/{reference}.1.bt2'
- *ref:* reference genome, e.g. 'reference/{reference}/{reference}.fa'

**Output(s):**

- *bam:* mapped read in BAM file, e.g. 'mapping/{reference}/original/{sample}.bam'

**Param(s):**

- *index:* name of reference, technically filename's path prefix, e.g. 'reference/{reference}/bowtie2_index/{reference}'
- *additional:* additional params

