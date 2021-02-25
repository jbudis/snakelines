Useful bash aliases
===================

Snakemake is very flexible in workflow execution, see `its documentation <https://snakemake.readthedocs.io/en/stable/executable.html#all-options>`_ for detailed description.
Here you may find some useful aliases for multi-threaded execution, running on SGE cluster or visualization of individual steps of the pipeline.

* dsnake - shows rules to be executed, but do not start analysis (--dryrun)
* vsnake - visually examine rules to be executed and their dependencies
* snake  - runs standard SnakeLines pipeline in the current working directory
* qsnake - distribute snake jobs on the cluster - use for multi-threaded rules
* fsnake - distribute snake jobs on the cluster - use for single-threaded rules

Running SnakeLines
------------------

Edit your ~/.bash_aliases file, to load aliases at login to the system.

.. code-block:: bash

   vim ~/.bash_aliases

Add aliases for using pipelines locally. You can set up using conda with `USE_CONDA` variable and location for downloaded tools from Anaconda repository.

.. code-block:: bash

   THREADS_LOCAL=4     # Number of used threads, when running pipelines locally
   SNAKELINES_DIR=/data/snakelines/$USER   # Path to snakelines source files
   USE_CONDA=true
   CONDA_DIR=/data/snakelines/snakemake_repos
   if [ "$USE_CONDA" = true ] ; then
     alias basesnake='snakemake -d `pwd` --jobname {rulename}.{jobid} --reason --printshellcmds --snakefile $SNAKELINES_DIR/snakelines.snake' --use-conda --conda-prefix=$CONDA_DIR ; else
     alias basesnake='snakemake -d `pwd` --jobname {rulename}.{jobid} --reason --printshellcmds --snakefile $SNAKELINES_DIR/snakelines.snake'
   fi

   alias snake='basesnake --config threads=$THREADS_LOCAL --cores $THREADS_LOCAL'
   alias dsnake='basesnake --config threads=$THREADS_LOCAL --dryrun'

   function vsnake {
      dsnake $@ --rulegraph | dot | display
   }

You may also run scripts on computational cluster

.. code-block:: bash

   # SGE job scheduling system specific
   SGE_PRIORITY=0   # Priority of tasks on cluster
   SGE_THREADS=16   # Number of used threads, when running pipelines on cluster
   SGE_NODES=10     # Number of computational nodes in the cluster
   SGE_CORES=160    # Total number of cpus on all nodes in the cluster
   SGE_PARAMS="qsub -cwd -o log/ -e log/ -p $PRIORITY -r yes" # Additional parameters for the SGE scheduler
   SGE_QUEUE_MAIN=main.q
   SGE_CLUSTER_PARAMS="qsub -cwd -o log/ -e log/ -p $SGE_PRIORITY -r yes"

   function clustersnake {
       NOW=`date "+%y-%m-%d_%H-%M-%S"`
       LOGDIR=log/$NOW
       mkdir -p $LOGDIR
       rm -f log/last
       ln -f -s $NOW log/last
       basesnake --config threads=$1 --cluster "$SGE_CLUSTER_PARAMS -q $3 -l thr=$1 -o $LOGDIR -e $LOGDIR" --jobs $2 ${@:4}
   }

   alias qsnake='clustersnake 16 10 $SGE_QUEUE_MAIN'
   alias fsnake='clustersnake 1 160 $SGE_QUEUE_MAIN'


Next time, you will log into system, aliases would be ready.
When you do not want to re-login, you may reload aliases

.. code-block:: bash

   source ~/.bash_aliases

Logging
-------

Logs are stored in the logs/ directory, when running pipeline on the computational cluster using qsnake or fsnake alias.
Each pipeline run has its own directory named by time of the execution.
The last/ directory is link to the last executed analysis.
For example:
::

   logs/
      |-- 18-10-08_16-19-48/
      |-- 18-10-08_17-16-23/
      |-- last -> 18-10-08_17-16-23/

Summary logs of the last run may be examined using log alias.
You must be inside of project directory to find which logs to display.

.. code-block:: bash

   function log {

      TYPE=e
      if [ "$1" == "o" ]; then
         TYPE=o
      elif [ "$1" == "a" ]; then
         TYPE=
      fi

      cat `python -c 'import os; cwd = os.getcwd(); print("/".join(cwd.split("/")[:4]))'`/log/last/*.$TYPE* | less
   }

Command ``log e`` will display only error messages, ``log o`` messages on standard stream.

Example run
-----------

Assuming you have input file in the SnakeLines compatible project structure, you may start analysis using these aliases

.. code-block:: bash

   # Go to screen - analysis would not terminate, if connection fails
   screen

   # Run always in the project root directory
   cd /data/projects/example

   # First try dryrun, check if pipeline is correct
   dsnake --configfile config_variant_calling.yaml

   # Optionally visualise pipeline - but only on rules with small number of samples
   vsnake --configfile config_variant_calling.yaml

   # Run test analysis with one, small sample
   snake --configfile config_variant_calling.yaml

   # Distribute tasks for all samples on cluster
   ## For multi-threaded analysis
   qsnake --configfile config_variant_calling.yaml
   ## For single-threaded analysis
   fsnake --configfile config_variant_calling.yaml

Changing SGE queue
------------------

SGE engine supports organising computational nodes into groups, called queues.
You may specify, which queue you want to use in the `clustersnake` command.
Alternately, you may prepare your own aliases for each cluster queue.

For example, assume computational cluster with 8 computational nodes organised in groups:

* main.q - nodes 1, 2, 3, 4, 5, 6, 7, 8
* pat.q - nodes 1, 2, 3, 4
* mat.q - nodes 5, 6, 7, 8

Aliases for SnakeLines calls may be specified as:

.. code-block:: bash

   SGE_QUEUE_MAIN=main.q
   alias qsnake='clustersnake 16 10 $SGE_QUEUE_MAIN'
   alias fsnake='clustersnake 1 160 $SGE_QUEUE_MAIN'

   # Cluster specific
   SGE_QUEUE_MAT='mat.q'
   SGE_QUEUE_PAT='pat.q'

   alias qsnake.pat='clustersnake 16 4 $SGE_QUEUE_PAT'
   alias fsnake.pat='clustersnake 1 64 $SGE_QUEUE_PAT'
   alias qsnake.mat='clustersnake 16 4 $SGE_QUEUE_MAT'
   alias fsnake.mat='clustersnake 1 64 $SGE_QUEUE_MAT'

Call of aliases would generate:

.. code-block:: bash

   # Distribute tasks on all cluster nodes - 1-8
   qsnake --configfile config_variant_calling.yaml

   # Distribute tasks on 4 cluster nodes - 1-4
   qsnake.pat --configfile config_variant_calling.yaml

   # Distribute tasks on 4 cluster nodes - 5-8
   qsnake.mat --configfile config_variant_calling.yaml
