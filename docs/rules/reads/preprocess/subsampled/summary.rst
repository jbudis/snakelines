Seqtk - Subsample Reads
---------------------------

seqtk sample -s{params.seed} {input.r1} {params.n_reads} | gzip > {output.r1}
seqtk sample -s{params.seed} {input.r2} {params.n_reads} | gzip > {output.r2}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/preprocess/subsampled/seqtk.snake
- *Rule name:* seqtk__subsample_reads

**Input(s):**


