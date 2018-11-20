rule gatk__fix_vcf_header
-------------------------
located in: '<SnakeLines_dir>/rules/variant/report/calling/gatk.snake'

Adds missing sequence dictionary to VCF header. This job also generates VCF index (.idx).

Input(s):
	vcf: Raw VCF from variant caller vardict.
	fasta: Reference sequence.
Output(s):
	vcf: Fixed VCF.
	vcf_index: Fixed VCF index.

rule tabix__index_vcf
---------------------
located in: '<SnakeLines_dir>/rules/variant/report/calling/gatk.snake'

Create tabix intex on BGZF (bgzipped) VCF file.

Input(s):

rule picard__bed_to_interval_list
---------------------------------
located in: '<SnakeLines_dir>/rules/variant/report/calling/gatk.snake'

Conversion of BED file to GATK specific interval_list.

Input(s):
	bed: BED file
	seq_dict: sequence dictionary
Output(s):
	intervals: interval list

rule gatk__collect_variant_calling_metrics
------------------------------------------
located in: '<SnakeLines_dir>/rules/variant/report/calling/gatk.snake'

GATK tool for generating

Input(s):
	vcf: called variants
	vcf_index: index of called variants
	dbsnp: DBSNP in BGZF format.
	dbsnp_index: DBSNP index of BGZF format.
	intervals: Genomic regions of interest.
Output(s):
	filename: Text file with summary of calling metrics.

