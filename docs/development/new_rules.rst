Developing new rules
====================

All rules are stored in the `rules` directory.
They should contain only logic, how to generate output files.
Here we present conventions and best practices for developing new rules for SnakeLines.

Rules folder structure
----------------------

Rules are imported automatically, based on the `configuration <../user/configuration.html#adjust-rules-parameters>`_.
Therefore we recommend to place your source code into semantic folder structure, which is predictable and understandable by both user and developer.
For example it is easier to understand what `seqtk` tool does, if you put it into
::

    rules/reads/preprocess/subsampled/seqtk.snake

rather than
::

    rules/tools/seqtk.snake

The former path would be defined in the configuration as

.. code-block:: yaml

    reads:
        preprocess:
           subsampled:                 # Randomly select subset of reads
               method: seqtk           # Supported values: seqtk
               n_reads: 10             # Number of reads to select


Output and dependencies
-----------------------

SnakeLines internally requires information about the required output of each rule,
whether the outputs of the rule should be copied into the report directory, and (optionally) dependency of each rule.
This information is written in `src/dependency.yaml` file in human readable yaml format.
This way, every output of every rule is in one place.

The structure of this file directly copies the structure from the configuration file.
However, it is independent on the exact methods (for example `bowtie2` and `bwa` tools for mapping should produce files with the same names) and on the methods' parameters.

A small example part of the dependency file:

.. code-block:: yaml

    mapping:
        mapper:
            output:
                alignments: expand('mapping/{sr.reference}/original/{sr.sample}.bam', sr=pipeline.sample_references)

        report:
            quality_report:
                output:
                    quality_reports:
                        from: expand('mapping/{sr.reference}/{map_type}/stats-wgs/samples/{sr.sample}/report.pdf', sr=pipeline.sample_references, map_type=pipeline.postprocessed_map_type)
                        to: expand('{report_dir}/{sr.reference}/samples/{sr.sample}/mapping_quality.pdf', report_dir=config['report_dir'], sr=pipeline.sample_references)
                    summary_report:
                        from: expand('mapping/{reference}/{map_type}/stats-wgs/summary/report.pdf', reference=multisample_references, map_type=pipeline.postprocessed_map_type)
                        to: expand('{report_dir}/{reference}/summary/mapping_quality.pdf', report_dir=config['report_dir'], reference=multisample_references)
                depends:
                    - mapping/mapper


In this example, `mapping/mapper` pipeline generates alignments as BAM files for each sample from an internal object (this object stores all samples defined in `configuration <../user/configuration.html#define-set-of-samples>`_). Finally, there is `mapping/report/quality_report` pipeline defined that generates two types of .pdf reports (`quality_report` and `summary_report`), which are copied into the report directory (the data is generated to path defined in `from:` directive and if everything went correctly, it is copied to the path defined in `to:`).
This pipeline depends on the previous (`mapping/mapper`) pipeline - if the mapping/mapper pipeline is not defined, SnakeLines generates a warning, but tries to proceed (although it will likely fail, since this pipeline requires .bam files - outputs of `mapping/mapper` pipeline). Although it is allowed to define empty output of a pipeline (i.e. `output: []`), we would like to discourage such practise. If there is a rule that is used as a preliminary computation (whose output would never be needed in the end), this rule should be included directly in the snake file with a pipeline, that generates desired output (or use the snakemake `include:` directive if you need it on more locations - see how the `reads/conversion/unzip` rule is done).

If you create a directory for a new rule file (a subdirectory somewhere in the `rules` directory), you need to add the corresponding configuration into `src/dependency.yaml`.
For example, when we created new rules in file 'rules/mapping/report/quality_report/qualimap.snake', we added the 'mapping/report/quality_report' configuration block into the `src/dependency.yaml` file.


Naming of rules
---------------

Always name the rule by {tool}__{what_tool_does} convention.
This way, user always know, what tool is generating files and what is his purpose.
If you use your own script without name, use `custom` as tool name.
For example
::

   rule samtools__sort_mapped_reads:
   rule bowtie2__map_reads_to_reference:
   rule trimmomatic__trim_reads:
   rule custom__visualise_taxonomic_counts_as_barplot:


Rules of thumb
--------------

Be sure your rule contains

* docstring after the name of the rule
* named input files
* named output files
* log files for both error and standard output stream (stored in the log/ directory inside output files directory)
* threads (if applicable for the rule) - use value from configuration, such as in the example below
* bash command or python script

For example:

.. code-block:: yaml

    rule samtools__sort_mapped_reads:
    """
    Sort aligned reads according to mapped position on reference genome.
    :input ref: Reference genomic sequences in fasta format
    :input bam: Unordered mapped reads in bam format
    :output bam: Ordered mapped reads according to their location on reference genome
    """
    input:
        ref = 'reference/{reference}/{reference}.fa',
        bam = 'mapping/{{reference}}/{map_type}/{{sample}}.bam'.format(map_type=method_config['input_map_type'])
    output:
        bam = 'mapping/{reference}/sorted/{sample}.bam'
    log:
        out = 'mapping/{reference}/sorted/log/{sample}.log',
        err = 'mapping/{reference}/sorted/log/{sample}.err'
    threads:
        int(config['threads'])
    shell:
        """
        samtools sort \
            -o {output.bam} \
            --threads {threads} \
            --output-fmt BAM \
            --reference {input.ref} \
            {input.bam} \
        >  {log.out} \
        2> {log.err}
        """

When using bash script, be sure you use full parameter names, where applicable.
For example, --output-fmt is more informative than -O.

Code style guidelines
---------------------

Before pull request, check if your code meets these guidelines

* source code meets PEP-8 style guidelines - PyCharm will underline code that is problematic
* remove unused imports
* prefer docstrings before comments ('#')
* avoid long lines (> 80 chars) - when possible
* only fasta file should be in example references - avoid indices and other files that can be generated automatically
* keep only a single reference for any example (if possible)
* use method_config instead of config variable in snake files
* use format instead of %
* check grammar (especially misspelled words - they are detected and underlined by PyCharm)
* documentation pages for new or modified rules are updated

Method configuration
--------------------

Configuration for a rule in config.yaml would be accessible from the rule source code in the form of `method_config` dictionary.
For example,

.. code-block:: yaml

   reads:                           # Prepare reads and quality reports for downstream analysis
      preprocess:                   # Pre-process of reads, eliminate sequencing artifacts, contamination ...
         trimmed:                   # Remove low quality parts of reads
            method: trimmomatic     # Supported values: trimmomatic
            temporary: False        # If True, generated files would be removed after successful analysis
            crop: 500               # Maximal number of bases in read to keep. Longer reads would be truncated.
            quality: 20             # Minimal average quality of read bases to keep (inside sliding window of length 5)
            headcrop: 20            # Number of bases to remove from the start of read
            minlen: 35              # Minimal length of trimmed read. Shorter reads would be removed.

In the /rules/reads/preprocess/trimmed/trimmomatic.snake you may use dictionary `method_config` with these values:

.. code-block:: python

   method_config = {'temporary': False,
                    'crop': 500,
                    'quality': 20,
                    'headcrop': 20,
                    'minlen': 35}
