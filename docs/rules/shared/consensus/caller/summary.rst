Bcftools - Build Consensus Sequence
---------------------------------------

input:
    fa = 'reference/{reference}/{reference}.fa.fai',
    vcf = 'variant/'
output:
    bed = 'reference/{reference}/annotation/{panel}/regions.vardict.bed',
shell:

**Location**

- *Filepath:* <SnakeLines_dir>/rules/shared/consensus/caller/bcftools.snake
- *Rule name:* bcftools__build_consensus_sequence

**Input(s):**

