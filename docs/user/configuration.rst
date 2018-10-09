Configuration of pipeline
=========================

Supported pipelines are stored at the pipeline/ directory of SnakeLines scripts.
Pipelines typically contains a master file (Snakefile), several sub-workflows (Snakefile.<sub-workflow>) and a configuration file (config.yaml).
The master file and sub-workflows define, what files would be generated.
Rules, that define, how they would be generated are stored in the rules/ directory.
These rules may be easily configured adjusting the config.yaml configuration of the pipeline.

Configurations typically consists of several steps.
At first, define set of samples and reference genomes to analyse.
Next, define where to store summary reports of the analysis.
Finally, adjust parameters of rules in the pipeline.
See `example <#example-configuration>`_ at the end of this chapter.


Define set of samples
---------------------

Fastq files stored in the reads/original directory can be analysed together.
To analyse only a subset of samples configuration must start with `samples:` attribute.
The attribute may contain several group categories, each with its own reference and target panel set.
For example, you may specify to all sample with name started with `one` to be analysed against hg19 genome, one panel.
For `sureselect6` samples you may use different genome and panel.

.. code-block:: yaml

   samples:                  # List of sample categories to be analysed
      # Category with one panel
      - name: one.*          # Regex expression of sample names to be analysed (reads/original/one.*_R1.fastq.gz)
        reference: hg19      # Reference genome for reads in the category (reference/hg19/hg19.fa)
        panel: one           # Bed file with target region coordinates (reference/hg19/annotation/one/regions.bed)

      # Another category with sureselect6 panel
      - name: sureselect6.*  # Regex expression of sample names to be analysed (reads/original/sureselect6.*_R1.fastq.gz)
        reference: hg38      # Reference genome for reads in the category (reference/hg19/hg19.fa)
        panel: sureselect6   # Bed file with target region coordinates (reference/hg19/annotation/sureselect6/regions.bed)




Example configuration
---------------------

Example for the read mapping pipeline.

.. literalinclude:: ../_static/configuration/mapping.yaml
   :language: yaml