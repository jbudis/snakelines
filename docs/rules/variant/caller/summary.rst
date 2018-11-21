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


Undocumented rules
------------------
WARNING: found  1 undocumented rules:

- rule vardict__tsb_to_vcf is UNDOCUMENTED
