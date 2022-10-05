eval "$(conda shell.bash hook)"
conda activate freyja-env
freyja update --outdir $1
