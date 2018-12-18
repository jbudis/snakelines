Fastqc - Html Summary For Joined Reads
------------------------------------------

Aggregate quality control results from all FastQC reports and generate summary HTML table with
PASS/WARN/FAIL icons

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/report/quality_report/fastqc.snake
- *Rule name:* fastqc__html_summary_for_joined_reads

**Input(s):**

- *fastq_files:* All joined Fastq files to quality report

**Output(s):**

- *summary:* HTML report with summary statistics of multiple fastq files

Fastqc - Html Summary For Paired Reads
------------------------------------------

Aggregate quality control results from all FastQC reports and generate summary HTML table with
PASS/WARN/FAIL icons

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/report/quality_report/fastqc.snake
- *Rule name:* fastqc__html_summary_for_paired_reads

**Input(s):**

- *fastq_files:* All Fastq files to quality report

**Output(s):**

- *summary:* HTML report with summary statistics of multiple fastq files

Fastqc - Quality Report
---------------------------

Generate HTML report with plots assessing various quality control aspects of NGS reads

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/report/quality_report/fastqc.snake
- *Rule name:* fastqc__quality_report

**Input(s):**

- *reads:* Sequenced reads in fastq format

**Output(s):**

- *html:* Quality report in HTML format
- *txt:* Quality report in text format suitable for automated processing

