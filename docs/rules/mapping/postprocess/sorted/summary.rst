Samtools - Sort Mapped Reads
--------------------------------

samtools sort \
    -o {output.bam} \
    --threads {threads} \
    --output-fmt BAM \
    --reference {input.ref} \
    {input.bam} \
>  {log.out} \
2> {log.err}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/sorted/samtools.snake
- *Rule name:* samtools__sort_mapped_reads

**Input(s):**


