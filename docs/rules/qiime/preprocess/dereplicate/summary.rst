Qiime - Vsearch Dereplicate Sequences
-----------------------------------------

Dereplicate sequence data and create a feature table and feature representative sequences.

Calls 'qiime vsearch dereplicate-sequences'.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/dereplicate/vsearch_dereplicate_sequences.snake
- *Rule name:* qiime__vsearch_dereplicate_sequences

**Input(s):**

- *qza:* qiime artifact of type SampleData[JoinedSequencesWithQuality] | SampleData[SequencesWithQuality] | SampleData[Sequences], e.g. reads/qiime/joined-table.qza

**Output(s):**

- *table:* qiime artifact of type FeatureTable[Frequency], e.g. reads/qiime/joined-table.qza
- *rep_seqs:* qiime artifact of type FeatureData[Sequence], e.g. reads/qiime/joined-rep-seqs.qza

