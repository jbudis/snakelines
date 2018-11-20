rule finalise__preprocess_reads
-------------------------------
located in: '<SnakeLines_dir>/pipeline/standard_analysis/preprocess/Snakefile'

Pipeline filter and revise fastq files from the reads/original directory. All preprocess steps, such as trimming,
deduplication and filtering of contamination reads are defined in the configuration file, part 'preprocess'.
Preprocess steps are executed gradually, i.e. output reads of one step are input to the following step.

Input(s):
	reads: fastq files of all preprocess types that are not configured as temporary

