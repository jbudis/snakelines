Vardict - Create Wgs Bed File
---------------------------------

Creates bed file with whole genomic regions in reference fasta file. Simple way to unite WGS and panel analysis.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/caller/vardict.snake
- *Rule name:* vardict__create_wgs_bed_file

**Input(s):**

Vardict - Prepare Bed File
------------------------------

Vardict throws error for bed files with more than 4 columns, rule therefore cut other columns

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/caller/vardict.snake
- *Rule name:* vardict__prepare_bed_file

**Input(s):**

Vardict - Call Germline Variants
------------------------------------

Identify small variation (SNP and indels) from the mapped reads.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/caller/vardict.snake
- *Rule name:* vardict__call_germline_variants

**Input(s):**

Vardict - Test Strand Bias
------------------------------

Perform statistical testing of strand bias.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/caller/vardict.snake
- *Rule name:* vardict__test_strand_bias

**Input(s):**

Vardict - Tsb To Vcf
------------------------

Convert intermediate files of Vardict to VCF format

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/caller/vardict.snake
- *Rule name:* vardict__tsb_to_vcf

**Input(s):**

