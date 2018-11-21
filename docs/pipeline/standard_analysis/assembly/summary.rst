Assemble Reads
-----------------------------

Pipeline filter and revise fastq files from the reads/original directory. All preprocess steps, such as trimming,
deduplication and filtering of contamination reads are defined in the configuration file, part 'preprocess'.
Preprocess steps are executed gradually, i.e. output reads of one step are input to the following step.

**Location**

- *Filepath:* <SnakeLines_dir>/pipeline/standard_analysis/assembly/Snakefile
- *Rule name:* finalise__assemble_reads

**Input(s):**

- *contigs:* assembled contigs in Fasta format
- *overlap_graphs:* contig overlap graph
- *quality_reports:* quality control metrics of assembled contigs

**Output(s):**

- *overlap_graphs:* contig overlap graph (in the report directory)
- *quality_reports:* quality control metrics of assembled contigs (in the report directory)

