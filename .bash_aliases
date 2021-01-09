THREADS_LOCAL=3     # Number of used threads, when running pipelines locally
SNAKELINES_DIR=/home/adog/snakelines/   # Path to snakelines source files
USE_CONDA=true

if [ $USE_CONDA = true ] ; then
  alias basesnake='snakemake -d `pwd` --jobname {rulename}.{jobid} --reason --printshellcmds  --snakefile $SNAKELINES_DIR/snakelines.snake --use-conda --conda-prefix=$CONDA_DIR'
else
  alias basesnake='snakemake -d `pwd` --jobname {rulename}.{jobid} --reason --printshellcmds --snakefile $SNAKELINES_DIR/snakelines.snake'
fi

alias snake='basesnake --config threads=$THREADS_LOCAL --cores $THREADS_LOCAL'
alias dsnake='basesnake --config threads=$THREADS_LOCAL --dryrun'

function vsnake {
   dsnake $@ --rulegraph | dot | display
}

export CLAIR='/home/adog/Clair/clair.py'

