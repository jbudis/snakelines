Custom - Summarize Taxonomic Counts Into Tsv Table
------------------------------------------------------

Take taxonomic counts and proportions from several samples and merge them together into a single table.
Also calculate number of reads per sample for each taxonomic unit.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/report/taxonomic_counts/count_table/custom.snake
- *Rule name:* custom__summarize_taxonomic_counts_into_tsv_table

**Input(s):**

- *krns:* Taxonomy files for each analysed sample in standardised format suitable for Krona graph generation

**Output(s):**

- *summary_xlsx:* Aggregated counts for each discovered taxonomy in Excel format
- *summary_tsv:* Aggregated counts for each discovered taxonomy in tabular format
- *subtax_xlsx:* Aggregated counts for each discovered taxonomy and its subtaxes in Excel format
- *subtax_tsv:* Aggregated counts for each discovered taxonomy and its subtaxes in tabular format
- *subtax_xlsx:* Normalised (counts for sample sum to 1) aggregated counts for each discovered taxonomy and its subtaxes in Excel format
- *subtax_tsv:* Normalised (counts for sample sum to 1) aggregated counts for each discovered taxonomy and its subtaxes in tabular format

Custom - Extract Taxonomic Level From Taxonomic Table
---------------------------------------------------------

Extract read counts for organisms at selected taxonomic level, such as genus, species.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/report/taxonomic_counts/count_table/custom.snake
- *Rule name:* custom__extract_taxonomic_level_from_taxonomic_table

**Input(s):**

- *subtax_tsv:* Aggregated counts for each discovered taxonomy and its subtaxes in tabular format

**Output(s):**

- *tax_level_abs:* Aggregated counts for single taxonomic level, e.g. species, order
- *tax_level_rel:* Normalised (counts for a sample sum to 1) aggregated counts for single taxonomic level, e.g. species, order

