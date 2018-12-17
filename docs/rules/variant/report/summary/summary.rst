Custom - Summary Report
---------------------------

Summarizing key metrics from sequencing (FASTQC),
mapping (Qualimap) and variant calling (GATK) in a single report file.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/variant/report/summary/custom.snake
- *Rule name:* custom__summary_report

**Input(s):**

- *fastqc_original:* FASTQC report on original reads.
- *fastqc_trimmed:* FASTQC report on trimmed reads.
- *qualimap_txt:* Plain text Qualimap report on deduplicated mapped reads.
- *qualimap_html:* HTML Qualimap report on deduplicated mapped reads.
- *calling_stats:* GATK variant calling metrics.

**Output(s):**

- *xml:* Report for further computer processing.
- *html:* Formatted report which is used to generate PDF report.
- *pdf:* Human readable report in form of PDF document.

