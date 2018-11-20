rule finalise__identify_viruses
-------------------------------
located in: '<SnakeLines_dir>/pipeline/viral_identification/assembly_based/Snakefile'

Pipeline annotates assembled contigs from preprocessed reads. Annotations, such as homologies with
reference databases, coverage and results of prediction tools are presented in sortable, filterable HTML table.

Input(s):
	summary_html: summary HTML report with annotations of assembled contigs
Output(s):
	summary_html: summary HTML report with annotations of assembled contigs (in the report directory)

