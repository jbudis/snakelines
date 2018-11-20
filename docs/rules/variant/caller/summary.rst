rule vardict__create_wgs_bed_file
---------------------------------
located in: '<SnakeLines_dir>/rules/variant/caller/vardict.snake'

Creates bed file with whole genomic regions in reference fasta file. Simple way to unite WGS and panel analysis.

Input(s):

rule vardict__prepare_bed_file
------------------------------
located in: '<SnakeLines_dir>/rules/variant/caller/vardict.snake'

Vardict throws error for bed files with more than 4 columns, rule therefore cut other columns

Input(s):

rule vardict__call_germline_variants
------------------------------------
located in: '<SnakeLines_dir>/rules/variant/caller/vardict.snake'

Identify small variation (SNP and indels) from the mapped reads.

Input(s):

rule vardict__test_strand_bias
------------------------------
located in: '<SnakeLines_dir>/rules/variant/caller/vardict.snake'

Perform statistical testing of strand bias.

Input(s):

undocumented rules
------------------
WARNING: found  1 undocumented rules:
	- rule vardict__tsb_to_vcf is UNDOCUMENTED
