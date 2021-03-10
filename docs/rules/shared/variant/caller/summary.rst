Vardict - Create Wgs Bed File
---------------------------------

Creates bed file with whole genomic regions in reference fasta file. Simple way to unite WGS and panel analysis.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__create_wgs_bed_file

**Input(s):**

Vardict - Prepare Bed File
------------------------------

Vardict throws error for bed files with more than 4 columns, rule therefore cut other columns

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__prepare_bed_file

**Input(s):**

Vardict - Call Germline Variants
------------------------------------

Identify small variation (SNP and indels) from the mapped reads.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__call_germline_variants

**Input(s):**

Vardict - Test Strand Bias
------------------------------

Perform statistical testing of strand bias.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__test_strand_bias

**Input(s):**

Vardict - Tsb To Vcf
------------------------

Convert intermediate files of Vardict to VCF format

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__tsb_to_vcf

**Input(s):**

Clair - Variant Call Mapped Reads
-------------------------------------

Identify small variants from reads mapped to a reference.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/clair.snake
- *Rule name:* clair__variant_call_mapped_reads

**Input(s):**

- *ref:* Reference genomic sequences in fasta format
- *bam:* Unordered mapped reads in bam format

**Output(s):**

- *bam:* Ordered mapped reads according to their location on reference genome

Vcfcat - Merge Vcf Files
----------------------------

Concatente all vcf files into one

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/clair.snake
- *Rule name:* vcfcat__merge_vcf_files

**Input(s):**

Medaka - Variant Call Mapped Reads
--------------------------------------

Identify small variants from reads mapped to a reference.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/medaka.snake
- *Rule name:* medaka__variant_call_mapped_reads

**Input(s):**

- *ref:* Reference genomic sequences in fasta format
- *bam:* Unordered mapped reads in bam format

**Output(s):**

- *bam:* Ordered mapped reads according to their location on reference genome

Freebayes - Call Mapped Reads
---------------------------------

Identify small variation (SNP and indels) from the mapped reads.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/freebayes.snake
- *Rule name:* freebayes__call_mapped_reads

**Input(s):**

