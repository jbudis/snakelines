rule finalise__call_germline_variants
-------------------------------------
located in: '<SnakeLines_dir>/pipeline/variant_calling/germline/Snakefile'

Pipeline identify SNPs and small insertions/deletions from mapped reads.

Input(s):
	variant: SNP and indel variants in VCF file
	report: Summary PDF report
Output(s):
	variant: SNP and indel variants in VCF file (in the report directory)
	report: Summary PDF report (in the report directory)

