Blast - Find Homologues
---------------------------

Find homologues for input sequences stored in fasta file in query reference database.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/contig_based/blast.snake
- *Rule name:* blast__find_homologues

**Input(s):**

- *blast_db:* Blast index of reference sequences, made by makeblastdb command
- *fasta:* Genomic sequences to be classified

**Output(s):**

- *tsv:* Homologues for input sequences in tabular format

**Param(s):**

- *nohead:* Auxiliary file for immediate results
- *header:* Auxiliary file for immediate results
- *max_target_seqs:* Maximal number of homologue sequences to be reported for each query sequence
- *blast_binary:* Blast tool to use - this differ in types of input files (nucleotide/protein)
- *blast_database:* Prefix name of the input blast database files

Blast - Annotate With Taxonomy
----------------------------------

Append additional column to the blast tabular output with the taxonomy of homology sequence

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/contig_based/blast.snake
- *Rule name:* blast__annotate_with_taxonomy

**Input(s):**

- *blast:* Blast tabular output
- *taxes:* Preprocessed taxonomies in .pickle

**Output(s):**

- *blast:* Blast tabular output with additional column with the taxonomy of homology sequence

Blast - Prepare Reference Index For Nucleotide
--------------------------------------------------

Creates blast indices from genomic sequences.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/contig_based/blast.snake
- *Rule name:* blast__prepare_reference_index_for_nucleotide

**Input(s):**

- *fa:* Genomic nucleotide sequences in fasta format

**Output(s):**

- *nhr:* Part of blast index
- *nin:* Part of blast index
- *nsq:* Part of blast index

Blast - Prepare Reference Index For Protein
-----------------------------------------------

Creates blast indices from protein sequences.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/classification/contig_based/blast.snake
- *Rule name:* blast__prepare_reference_index_for_protein

**Input(s):**

- *fa:* Genomic protein sequences in fasta format

**Output(s):**

- *phr:* Part of blast index
- *pin:* Part of blast index
- *psq:* Part of blast index

