Qiime - Feature Table Summarize
-----------------------------------

Create qiime visualization of clustered feature table, calls 'qiime feature-table summarize'

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/report/prepare/cluster_feature_table/feature_table_summarize.snake
- *Rule name:* qiime__feature_table_summarize

**Input(s):**

- *table:* qiime artifact of type FeatureTable[Frequency] with similarity, e.g. reads/qiime/joined-table-dn-99.qza
- *metadata:* TSV metadata filepath, e.g. description/joined-sample-metadata.tsv

**Output(s):**

- *qzv:* qiime visualization file, can be exported to HTML, e.g. reads/qiime/joined-table-dn-99.qzv

