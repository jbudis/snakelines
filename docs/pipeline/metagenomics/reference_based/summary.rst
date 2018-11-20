rule finalise__classify_reads_reference_based
---------------------------------------------
located in: '<SnakeLines_dir>/pipeline/metagenomics/reference_based/Snakefile'

Pipeline classify preprocess reads against chosen reference sequences. According to their taxonomic labels,
number of reads per taxonomic units are summarized into graphical and tabular reports.

Input(s):
	kronas: interactive multi-level pie plots of taxonomic counts
	tables: summary tables with number and proportions of reads mapped to taxonomic units
	barplots: barplots with number and proportions of reads mapped to taxonomic units
Output(s):
	kronas: interactive multi-level pie plots of taxonomic counts (in the report directory)
	tables: summary tables with number and proportions of reads mapped to taxonomic units (in the report directory)
	barplots: barplots with number and proportions of reads mapped to taxonomic units (in the report directory)

