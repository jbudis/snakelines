Qiime - Vsearch Cluster Features De Novo
--------------------------------------------

Creates clusters of sequences (OTUs) with {similarity}% similarity

Calls 'qiime vsearch cluster-features-de-novo'.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/clustered/vsearch_cluster_features_de_novo.snake
- *Rule name:* qiime__vsearch_cluster_features_de_novo

**Input(s):**

- *table:* qiime artifact of type FeatureTable[Frequency], e.g. reads/qiime/joined-table.qza
- *rep_seqs:* qiime artifact of type FeatureData[Sequence], e.g. reads/qiime/joined-rep-seqs.qza

**Output(s):**

- *table:* qiime artifact of type FeatureTable[Frequency] with similarity, e.g. reads/qiime/joined-table-dn-99.qza
- *rep_seqs:* qiime artifact of type FeatureData[Sequence] with similarity, e.g. reads/qiime/joined-rep-seqs-dn-99.qza

