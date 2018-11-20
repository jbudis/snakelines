rule qualimap__mapping_quality_report_accross_reference
-------------------------------------------------------
located in: '<SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake'

Generate summary statistics for mapped reads stored in BAM files. Statistics are calculated across
whole reference genome.

Input(s):

rule qualimap__mapping_quality_report_accross_panel
---------------------------------------------------
located in: '<SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake'

Generate summary statistics for mapped reads stored in BAM files. Statistics are calculated across
genomic regions specified in the BED file.

Input(s):

rule qualimap__summarize_quality_reports
----------------------------------------
located in: '<SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake'

Aggregate results from individual bamqc results to a single summary report

Input(s):

