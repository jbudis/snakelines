Bowtie2 - Map Reads To Reference
------------------------------------

For input preprocessed reads bowtie2 finds the most similar genomic region in the provided reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/bowtie2.snake
- *Rule name:* bowtie2__map_reads_to_reference

**Input(s):**

- *r1:* gzipped fastq file with left reads, e.g. 'reads/%s/{sample}_R1.fastq.gz'
- *r2:* gzipped fastq file with right reads, e.g. 'reads/%s/{sample}_R2.fastq.gz'
- *index:* reference index (created by rule bowtie2__prepare_index), e.g. 'reference/{reference}/bowtie2_index/{reference}.1.bt2'
- *ref:* reference genome, e.g. 'reference/{reference}/{reference}.fa'

**Output(s):**

- *bam:* mapped read in BAM file, e.g. 'mapping/{reference}/original/{sample}.bam'

**Param(s):**

- *index:* name of reference, technically filename's path prefix, e.g. 'reference/{reference}/bowtie2_index/{reference}'
- *additional:* additional params
- *concordant:* only concordant reads

Bwa - Map Reads To Reference
--------------------------------

For input preprocessed reads bwa finds the most similar genomic region in the provided reference genome, using bwa mem algorithm.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/bwa.snake
- *Rule name:* bwa__map_reads_to_reference

**Input(s):**

- *r1:* gzipped fastq file with left reads, e.g. 'reads/%s/{sample}_R1.fastq.gz'
- *r2:* gzipped fastq file with right reads, e.g. 'reads/%s/{sample}_R2.fastq.gz'
- *index:* reference index (created by rule bwa__prepare_index), e.g. 'reference/{reference}/bwa_index/{reference}.bwt'

**Output(s):**

- *bam:* mapped read in BAM file, e.g. 'mapping/{reference}/original/{sample}.bam'

**Param(s):**

- *index:* name of reference, technically filename's path prefix, e.g. 'reference/{reference}/bwa_index/{reference}'
- *additional:* additional params
- *concordant:* more strict conditions for scoring options should help to map only concordant reads

Bismark - Map Methyl Seq Reads To Reference
-----------------------------------------------

For input preprocessed reads treated by Bisulfide Bismark finds the most similar genomic region in the provided reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/bismark.snake
- *Rule name:* bismark__map_methyl_seq_reads_to_reference

**Input(s):**

- *r1:* gzipped fastq file with left reads, e.g. 'reads/%s/{sample}_R1.fastq.gz'
- *r2:* gzipped fastq file with right reads, e.g. 'reads/%s/{sample}_R2.fastq.gz'
- *ct_index:* CT reference index (created by rule bismark__prepare_index)
- *ga_index:* GA reference index (created by rule bismark__prepare_index)
- *ref:* reference genome, i.e. 'reference/{reference}/{reference}.fa'

**Output(s):**

- *bam:* mapped reads in BAM format
- *alignment_report:* Control report with basic mapping statistics generated by Bismark

**Param(s):**

- *aux_bam:* bam file that have preset name by Bismark. Has to be renamed to match downstream analysis.
- *aux_alignment_report:* alignment report that have preset name by Bismark. Better to move to separate folder to keep bam directory clean
- *bam_dir:* directory with the output bam file
- *reference_dir:* directory with reference fasta, i.e. 'reference/{reference}'

