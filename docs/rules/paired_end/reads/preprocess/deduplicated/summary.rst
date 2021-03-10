Fastuniq - Deduplicate Reads
--------------------------------

Remove fragments with identical sequences that are usually consequence of extensive PCR multiplication

**Location**

- *Filepath:* <SnakeLines_dir>/rules/paired_end/reads/preprocess/deduplicated/fastuniq.snake
- *Rule name:* fastuniq__deduplicate_reads

**Input(s):**

- *r1:* Left side of sequenced fragments in gzipped fastq format
- *r2:* Right side of sequenced fragments in gzipped fastq format

**Output(s):**

- *r1:* Left side of fragments without duplications in gzipped fastq format
- *r2:* Right side of fragments without duplications in gzipped fastq format
- *pair_description:* Auxiliary file with names for paired files
- *unzipped_r1:* Auxiliary file with unzipped left side of sequenced fragments
- *unzipped_r2:* Auxiliary file with unzipped right side of sequenced fragments

