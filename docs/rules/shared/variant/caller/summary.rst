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

Config parameters:
    hard_filter:
        min_nonref_allele_freq
        min_alternate_count
        min_map_quality

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/vardict.snake
- *Rule name:* vardict__call_germline_variants

**Input(s):**

- *bam:* Unordered mapped reads in bam format
- *bai:* Ordered and indexed mapped reads in bam format
- *bed:* Bed file containing whole genomic regions from the reference fasta
- *reffasta:* Reference genomic sequences in fasta format
- *reffaidx:* Indexed reference genomic sequences in fasta.fai format

**Output(s):**

- *vcf:* Variants in the Variant Call Format (VCF)

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

Config parameters:
    hard_filter:
        min_nonref_allele_freq
        min_coverage

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/clair.snake
- *Rule name:* clair__variant_call_mapped_reads

**Input(s):**

- *bam:* Unordered mapped reads in bam format
- *bai:* Ordered and indexed mapped reads in bam format
- *reffasta:* Reference genomic sequences in fasta format
- *reffaidx:* Indexed reference genomic sequences in fasta.fai format

**Output(s):**

- *vcf:* Variants in the Variant Call Format (VCF)

Vcfcat - Merge Vcf Files
----------------------------

Concatente all vcf files into one

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/clair.snake
- *Rule name:* vcfcat__merge_vcf_files

**Input(s):**

Medaka - Variant Call Mapped Reads
--------------------------------------

Identify small variants from reads mapped to a reference. Mainly for diploid samples.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/medaka.snake
- *Rule name:* medaka__variant_call_mapped_reads

**Input(s):**

- *bam:* Unordered mapped reads in bam format
- *reffasta:* Reference genomic sequences in fasta format
- *reffaidx:* Indexed reference genomic sequences in fasta.fai format

**Output(s):**

- *vcf:* Variants in the Variant Call Format (VCF)

Freebayes - Variant Call Mapped Reads
-----------------------------------------

Identify small variation (SNP and indels) from the mapped reads.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/variant/caller/freebayes.snake
- *Rule name:* freebayes__variant_call_mapped_reads

**Input(s):**

