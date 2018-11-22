Map Reads To Reference
-------------------------------------

Pipeline maps preprocessed reads to the selected reference sequences. Mapped reads in bam files
are further processed - sorting, deduplication, realignment... Final mapping statistics and quality
reports are generated in the end.

**Location**

- *Filepath:* <SnakeLines_dir>/pipeline/standard_analysis/mapping/Snakefile.rule
- *Rule name:* finalise__map_reads_to_reference

**Input(s):**

- *alignments:* mapped sorted reads in BAM file (.bam) (for each of the sample references and map types)
- *bam_indices:* BAM file indices (.bam.bai)
- *quality_reports:* mapping statistics retrieved from BAM file (.pdf)
- *summary_report:* aggregated mapping statistics (.pdf)

**Output(s):**

- *quality_reports:* mapping statistics retrieved from BAM file (.pdf) (in the report directory)
- *summary_report:* aggregated mapping statistics (.pdf) (in the report directory)

