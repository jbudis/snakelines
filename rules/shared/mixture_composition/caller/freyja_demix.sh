eval "$(conda shell.bash hook)"
conda activate freyja-env
freyja demix $1 $2 --output $3 --meta /home_pfs/data/narcos/production/scripts/snakelines/dev/rules/shared/mixture_composition/caller/lineages/curated_lineages.json --barcodes /home_pfs/data/narcos/production/scripts/snakelines/dev/rules/shared/mixture_composition/caller/lineages/usher_barcodes.csv
