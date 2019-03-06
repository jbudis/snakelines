Entrez - Download Sequences By Genbank Id
---------------------------------------------

Download sequences from NCBI Genbank database according to the list of enumerated genbank ids

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reference/download/entrez.snake
- *Rule name:* entrez__download_sequences_by_genbank_id

**Input(s):**

**Output(s):**

- *reference:* Downloaded sequences in FASTA format
- *taxonomy:* Taxonomies for downloaded sequences

**Param(s):**

- *method_config:* Configuration of entrez. It has to contain list of Genbank identifiers for sequences to download
- *email:* Inform NCBI who you are to contact you in case of excessive use. Otherwise they may block your access directly.

