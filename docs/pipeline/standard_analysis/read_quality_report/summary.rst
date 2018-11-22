Quality Report
-----------------------------

Make FastQC quality reports for all reads of defined types. Should be used as the first assessment of
read quality of sequencing run. After that, preprocess pipeline may be utilized to eliminate sequencing artefacts.

**Location**

- *Filepath:* <SnakeLines_dir>/pipeline/standard_analysis/read_quality_report/Snakefile
- *Rule name:* finalise__quality_report

**Input(s):**

- *fastqcs:* fastQC quality reports of all preprocess types that are configured to be generated
- *reports:* Summary table of fastQC reports

**Output(s):**

- *fastqcs:* fastQC quality reports of all preprocess types that are configured to be generated (in the report directory)
- *reports:* Summary table of fastQC reports (in the report directory)

