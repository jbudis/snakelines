Gatk - Fix Vcf Header
-------------------------

Adds missing sequence dictionary to VCF header. This job also generates VCF index (.idx).

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/report/calling/gatk.snake
- *Rule name:* gatk__fix_vcf_header

**Input(s):**

- *vcf:* Raw VCF from variant caller vardict.
- *fasta:* Reference sequence.

**Output(s):**

- *vcf:* Fixed VCF.
- *vcf_index:* Fixed VCF index.

Tabix - Index Vcf
---------------------

Create tabix intex on BGZF (bgzipped) VCF file.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/report/calling/gatk.snake
- *Rule name:* tabix__index_vcf

**Input(s):**


Picard - Bed To Interval List
---------------------------------

Conversion of BED file to GATK specific interval_list.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/report/calling/gatk.snake
- *Rule name:* picard__bed_to_interval_list

**Input(s):**

- *bed:* BED file
- *seq_dict:* sequence dictionary

**Output(s):**

- *intervals:* interval list

Gatk - Collect Variant Calling Metrics
------------------------------------------

GATK tool for generating

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/report/calling/gatk.snake
- *Rule name:* gatk__collect_variant_calling_metrics

**Input(s):**

- *vcf:* called variants
- *vcf_index:* index of called variants
- *dbsnp:* DBSNP in BGZF format.
- *dbsnp_index:* DBSNP index of BGZF format.
- *intervals:* Genomic regions of interest.

**Output(s):**

- *filename:* Text file with summary of calling metrics.

