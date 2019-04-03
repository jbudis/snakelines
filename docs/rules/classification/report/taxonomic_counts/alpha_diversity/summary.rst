Custom - Alpha Diversity
----------------------------

Compute aplha diversity at selected taxonomic level, such as genus, species.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/report/taxonomic_counts/alpha_diversity/custom.snake
- *Rule name:* custom__alpha_diversity

**Input(s):**

- *tax_level_abs:* Aggregated counts for single taxonomic level, e.g. species, order

**Output(s):**

- *alpha_div:* Alpha diversities for single taxonomic level, e.g. species, order; either '_norm'-alized to the same counts or kept '_pure'.

