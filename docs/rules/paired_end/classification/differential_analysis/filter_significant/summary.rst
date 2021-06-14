Custom - Filter Significantly Expressed Transcripts
-------------------------------------------------------

Select transcripts with significant change in expression according to user specified filters.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/paired_end/classification/differential_analysis/filter_significant/custom.snake
- *Rule name:* custom__filter_significantly_expressed_transcripts

**Input(s):**

- *table:* TSV table with statistical evaluation of change in expression
- *design:* TSV table describing separation of samples into groups for differential expression

**Output(s):**

- *table:* TSV table with transcripts with significant change in expression

**Param(s):**

- *max_fdr:* Maximal value of fold discovery change for transcript to be reported
- *min_fold_change:* Minimal value of fold change for transcript to be reported
- *reproducible_expression:* At least one read must be mapped to transcript in all samples from over-expressed group to be reported

