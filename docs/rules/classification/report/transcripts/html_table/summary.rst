Custom - Visualise Transcriptomic Counts In Html Table
----------------------------------------------------------

Take transcriptomic counts from several samples and merge them together into a single table.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/report/transcripts/html_table/custom.snake
- *Rule name:* custom__visualise_transcriptomic_counts_in_html_table

**Input(s):**

- *table:* TSV table with statistical evaluation of change in expression
- *desc:* Description of reference sequences
- *template:* HTML template with basic report outline
- *annotations:* TSV files with attributes for annotated transcripts

**Output(s):**

- *html:* HTML page with sortable, filterable table of transcriptomic results

