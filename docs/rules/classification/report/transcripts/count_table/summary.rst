Custom - Summarize Transcriptomic Counts Into Tsv Table
-----------------------------------------------------------

Take transcriptomic counts from several samples and merge them together into a single table.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/report/transcripts/count_table/custom.snake
- *Rule name:* custom__summarize_transcriptomic_counts_into_tsv_table

**Input(s):**

- *krns:* Count files for each analysed sample in standardised format suitable for Krona graph generation

**Output(s):**

- *summary_xlsx:* Aggregated counts for each transcript in Excel format
- *summary_tsv:* Aggregated counts for each transcript in tabular format

