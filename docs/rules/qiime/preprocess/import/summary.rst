Qiime - Import Manifest
---------------------------

Imports samples written down in metadata file, sequences should NOT be pair-end. Calls 'qiime tools import'

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/import/import_manifest.snake
- *Rule name:* qiime__import_manifest

**Input(s):**

- *manifest:* TSV manifest file in format 'sample-id,absolute-filepath,direction', e.g. reads/qiime/joined.manifest

**Output(s):**

- *qza:* qiime artifact of type SampleData[{params.data_type}], e.g. reads/qiime/joined-imported.qza

Qiime - Import Pair End Manifest
------------------------------------

Imports samples written down in metadata file, sequences MUST BE pair-end. Calls 'qiime tools import'

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/import/import_pair_end_manifest.snake
- *Rule name:* qiime__import_pair_end_manifest

**Input(s):**

- *manifest:* TSV manifest file in format 'sample-id,absolute-filepath,direction', e.g. reads/qiime/joined.manifest

**Output(s):**

- *qza:* qiime artifact of type SampleData[PairedEndSequencesWithQuality], e.g. reads/qiime/joined-imported.qza

