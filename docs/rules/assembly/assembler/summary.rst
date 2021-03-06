Spades - Assemble Reads Into Contigs
----------------------------------------

Assemble preprocessed reads into the larger genomic sequences, contigs. Also generate contig overlap graphs and
initial scaffolds. Use 'spades' program.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/assembly/assembler/spades.snake
- *Rule name:* spades__assemble_reads_into_contigs

**Input(s):**

- *r1:* gzipped fastq file with left reads
- *r2:* gzipped fastq file with right reads

**Output(s):**

- *fastg:* assembled graph in .fastg format
- *gfa:* assembly graph in .gfa format
- *contigs:* contigs in .fa file
- *scaffolds:* scaffolds in .fa file

**Param(s):**

- *outdir:* output directory (do not change)
- *contigs:* contig file generated by spades (do not change)
- *scaffolds:* scaffolds file generated by spades (do not change)
- *mode:* mode of operation of spades (do not change), extracted from config file
- *careful:* whether to use --careful parameter for spades (do not change), extracted from config file

Unicycler - Assemble Reads Into Contigs
-------------------------------------------

Assemble preprocessed reads into the larger genomic sequences, contigs. Also generate contig overlap graphs and
initial scaffolds. Use 'unicycler' program.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/assembly/assembler/unicycler.snake
- *Rule name:* unicycler__assemble_reads_into_contigs

**Input(s):**

- *r1:* gzipped fastq file with left reads
- *r2:* gzipped fastq file with right reads

**Output(s):**

- *gfa:* assembly graph in .gfa format as needed by downstream analysis
- *contigs:* assembled contigs in .fa file

**Param(s):**

- *outdir:* output directory (do not change)
- *contigs:* contig file generated by spades (do not change)
- *gfa:* assembly graph in .gfa format as generated by unicycler

