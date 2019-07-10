Custom - Prepare Description File
-------------------------------------

Extract ID and description from input reference sequences in FASTA format.
Assuming that id is the first word after '>'.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reference/attributes/desc/custom.snake
- *Rule name:* custom__prepare_description_file

**Input(s):**

- *reference:* Genomic sequences in FASTA format

**Output(s):**

- *desc:* Description TSV file - 1st column is id of sequence, 2nd is description

