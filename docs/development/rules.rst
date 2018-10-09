Developing new rules
====================

All rules are stored at the `/rules` directory.
They should contain only logic, how to generate output files.
Here we present conventions and best practices for developing new rules for SnakeLines.

Location of rules
-----------------

Rules are imported automatically, based on configuration.
See example at the `configuration page <../user/configuration.html#adjust-rules-parameters>`_.
Therefore you want to place your source code to subdirectories, that are predictable and understandable by both user and developer.
For example it is easier to understand what `seqtk` tool does, if you put it into
::

	/rules/reads/preprocess/subsampled/seqtk.snake

than
::

	/rules/tools/seqtk.snake

The first one would be defined in the configuration as

.. code-block:: yaml

   reads:
       preprocess:
           subsampled:                 # Randomly select subset of reads
               method: seqtk           # Supported values: seqtk
               n_reads: 10             # Number of reads to select

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

* comment after the name of the rule
* named input files
* named output files
* log files for both error and standard output stream (stored in the log/ directory inside output files directory)
* threads (if applicable for the rule) - use value from configuration, such as in the example below
* bash command or python script

For example:

.. code-block:: yaml

    rule samtools__sort_mapped_reads:
    '''
    Sort aligned reads according to mapped position on reference genome.
    '''
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
        '''
        samtools sort \
            -o {output.bam} \
            --threads {threads} \
            --output-fmt BAM \
            --reference {input.ref} \
            {input.bam} \
        >  {log.out} \
        2> {log.err}

When using bash script, be sure you use full parameter names, where applicable.
For example, --output-fmt is more informative than -O.

Method config
-------------

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
