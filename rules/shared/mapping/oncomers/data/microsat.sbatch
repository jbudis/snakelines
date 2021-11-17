#!/bin/bash
#SBATCH --job-name=microsat_${SID}
#SBATCH --partition=compute,viz,mic,gpu
#SBATCH --output=/home_pfs/data/projects/population/microsat/log/%j.log
#SBATCH --error=/home_pfs/data/projects/population/microsat/log/%j.err
#SBATCH --cpus-per-task=1

set -e
set -u

SIDS=$1

SIDFILE=/home_pfs/data/projects/population/scripts/microsat/sids/${SIDS}.sids
SID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SIDFILE})

hostname
python /home_pfs/data/projects/population/scripts/microsat/analysis.py ${SID}
