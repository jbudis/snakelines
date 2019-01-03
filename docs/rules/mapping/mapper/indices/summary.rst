Bismark - Prepare Index
---------------------------

Generate mapping index that is utilized by Bismark algorithm to map bisulfide treated reads to a reference genome

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/indices/bismark_index.snake
- *Rule name:* bismark__prepare_index

**Input(s):**

- *fa:* fasta reference genome

**Output(s):**

- *ct_index:* CT reference index (created by rule bismark__prepare_index)
- *ga_index:* GA reference index (created by rule bismark__prepare_index)

**Param(s):**

- *fadir:* root directory of the reference sequence

Bowtie2 - Prepare Index
---------------------------

Generate mapping index that is utilized by Bowtie2 algorithm to map reads to a reference genome

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/indices/bowtie2_index.snake
- *Rule name:* bowtie2__prepare_index

**Input(s):**

- *fa:* fasta reference genome

**Output(s):**

- *indecis:* 6 files with indeces. e.g '{fadir}/bowtie2_index/{sequence}.[1-4].bt2', '{fadir}/bowtie2_index/{sequence}.rev.[1,2].bt2'

**Param(s):**

- *index:* name of output reference, technically filename's path prefix, e.g. '{fadir}/bowtie2_index/{sequence}'

Bwa - Prepare Index
-----------------------

Generate mapping index that is utilized by BWA algorithm to map reads to a reference genome

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/mapper/indices/bwa_index.snake
- *Rule name:* bwa__prepare_index

**Input(s):**

- *fa:* fasta reference genome

**Output(s):**

- *indecis:* 5 files with indices. e.g '{fadir}/bwa_index/{sequence}.[amb,ann,bwt,pac,sa]'

**Param(s):**

- *index:* output filename's path prefix, e.g. '{fadir}/bwa_index/{sequence}'

