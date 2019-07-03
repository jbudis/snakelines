Qualimap - Mapping Quality Report Across Reference
------------------------------------------------------

Generate summary statistics for mapped reads stored in BAM files. Statistics are calculated across
whole reference genome.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__mapping_quality_report_across_reference

**Input(s):**

- *bam:* Mapped reads in bam format
- *bai:* Index to mapped reads for enable fast read retrieval from desired genomic region

**Output(s):**

- *html:* Quality report of mapped reads in HTML format
- *pdf:* Quality report of mapped reads in PDF format
- *text:* Quality report of mapped reads in format suitable for automated processing

Qualimap - Mapping Quality Report Across Panel
--------------------------------------------------

Generate summary statistics for mapped reads stored in BAM files. Statistics are calculated across
genomic regions specified in the BED file.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__mapping_quality_report_across_panel

**Input(s):**

- *bam:* Mapped reads in bam format
- *bai:* Index to mapped reads for enable fast read retrieval from desired genomic region
- *bed:* Genomic regions of interest in bed format

**Output(s):**

- *html:* Quality report of mapped reads in HTML format
- *pdf:* Quality report of mapped reads in PDF format
- *text:* Quality report of mapped reads in format suitable for automated processing

Qualimap - Summarize Quality Reports
----------------------------------------

Aggregate results from individual bamqc results to a single summary report

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__summarize_quality_reports

**Input(s):**

- *reports:* List of pdf qualimap reports to aggregate

**Output(s):**

- *html:* Aggregated quality report of mapped reads in HTML format
- *pdf:* Aggregated quality report of mapped reads in PDF format

