Edger - Identify Transcripts With Changed Expression
--------------------------------------------------------

Compare counts of selected groups of samples and find features with significant change across the groups

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/differential_analysis/edger.snake
- *Rule name:* edger__identify_transcripts_with_changed_expression

**Input(s):**

- *counts:* Aggregated counts for each transcript and sample in tabular format
- *metadata:* File with samples attributes - use for separating samples into compared groups

**Output(s):**

- *norm_counts:* TSV table with normalised number of reads (by sample read depth)
- *table:* TSV table with statistical evaluation of change of expression
- *design:* TSV table describing separation of samples into groups for differential expression

**Param(s):**

- *group_by:* attribute (in the metadata file header) that would split samples into two categories
- *batch:* attribute (in the metadata file header) that would split samples categories with similar batch bias effect (e.g. samples from same sequencing run)

