Pear - Join Read Pairs
--------------------------

Join paired reads into single end reads based on sequence overlap

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/preprocess/joined/pear.snake
- *Rule name:* pear__join_read_pairs

**Input(s):**

- *r1:* Left side of sequenced fragments in fastq format
- *r2:* Right side of sequenced fragments in fastq format

**Output(s):**

- *r1:* Left side of fragments without overlap in fastq format
- *r2:* Right side of fragments without overlap in fastq format
- *rm:* Paired reads with overlap joined into merged sequences in fastq format

Pear - Concat Joined With Single
------------------------------------

Unite joined reads, and reads that could not be joined into a single fastq file.
Sequence between reads without overlap would be filled with N symbol.

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/preprocess/joined/pear.snake
- *Rule name:* pear__concat_joined_with_single

**Input(s):**

- *r1:* Left side of fragments without overlap in fastq format
- *r2:* Right side of fragments without overlap in fastq format
- *rm:* Paired reads with overlap joined into merged sequences in fastq format

**Output(s):**

- *rc:* Concatenated file, first goes joined reads and then reads without overlap with N symbols between reads

