eval "$(conda shell.bash hook)"
conda activate freyja-env
freyja variants $1 --variants $2 --depths $3 --ref $4
