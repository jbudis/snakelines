Qiime - Create Metadata File
--------------------------------

Creates generic metadata file from samples filenames. Each metadata information is separated in filename by {sep}. Headers are generic.
E.g.:
    filename: 2017-dna-replicaA_RM.fastq.gz
    metadata:
        sample-id   col1    col2    col3
        #q2:types   categorical categorical categorical
        2017-dna-replicaA   2017    dna replicaA

Note: for several analysis, mainly for emporer plotting, it is better to create your own metadata file manually

**Location**

- *Filepath:* <SnakeLines_dir>/rules/qiime/preprocess/metadata/create_metadata_file.snake
- *Rule name:* qiime__create_metadata_file

**Input(s):**

- *rm:* relative path to (preprocessed) joined fastQ samples, e.g. reads/joined/*_RM.fastq.gz

**Output(s):**

- *metadata:* samples' metadata filepath in TSV format (tab-separated value), e.g. description/joined-sample-metadata.tsv

