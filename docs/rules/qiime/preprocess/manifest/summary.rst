Qiime - Create Joined Manifest
----------------------------------

Creates manifest file for joined samples, that it required later for importing those samples into qiime.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/manifest/create_joined_manifest.snake
- *Rule name:* qiime__create_joined_manifest

**Input(s):**

- *rm:* relative path to (preprocessed) joined fastQ samples, e.g. reads/joined/*_RM.fastq.gz

**Output(s):**

- *manifest:* TSV manifest file in format 'sample-id,absolute-filepath,direction', e.g. reads/qiime/joined.manifest

