rule finalise__map_reads_to_reference
-------------------------------------
located in: '<SnakeLines_dir>/pipeline/standard_analysis/mapping/Snakefile'

Pipeline maps preprocessed reads to the selected reference sequences. Mapped reads in bam files
are further processed - sorting, deduplication, realignment... Final mapping statistics and quality
reports are generated in the end.

Input(s):
	alignments: mapped sorted reads in BAM file (.bam) (for each of the sample references and map types)
	bam_indices: BAM file indices (.bam.bai)
	quality_reports: mapping statistics retrieved from BAM file (.pdf)
	summary_report: aggregated mapping statistics (.pdf)
Output(s):
	quality_reports: mapping statistics retrieved from BAM file (.pdf) (in the report directory)
	summary_report: aggregated mapping statistics (.pdf) (in the report directory)

rule finalise__preprocess_reads
-------------------------------
located in: '<SnakeLines_dir>/pipeline/standard_analysis/mapping/Snakefile.preprocess'

Pipeline filter and revise fastq files from the reads/original directory. All preprocess steps, such as trimming,
deduplication and filtering of contamination reads are defined in the configuration file, part 'preprocess'.
Preprocess steps are executed gradually, i.e. output reads of one step are input to the following step.

Input(s):
	reads: fastq files of all preprocess types that are not configured as temporary

rule finalise__quality_report
-----------------------------
located in: '<SnakeLines_dir>/pipeline/standard_analysis/mapping/Snakefile.read_quality_report'

Make FastQC quality reports for all reads of defined types. Should be used as the first assessment of
read quality of sequencing run. After that, preprocess pipeline may be utilize to eliminate sequencing artefacts.

Input(s):
	fastqcs: fastQC quality reports of all preprocess types that are configured to be generated
	reports: Summary table of fastQC reports
Output(s):
	fastqcs: fastQC quality reports of all preprocess types that are configured to be generated (in the report directory)
	reports: Summary table of fastQC reports (in the report directory)

