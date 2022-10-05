eval "$(conda shell.bash hook)"
conda activate freyja-env
freyja demix $1 $2 --output $3 $4 $5
