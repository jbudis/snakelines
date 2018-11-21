rule bowtie2__filter_reads_from_reference
-----------------------------------------
located in: '<SnakeLines_dir>/rules/reads/preprocess/decontaminated/bowtie2.snake'

Remove reads that do not map to the reference, and so may be caused by contamination in lab processing.
Alternatively, using keep: False configuration removes all fragments that belongs to reference, and so is suitable
to remove contamination caused by host with known genome, e.g. human fragments.

Input(s):
