#!/bin/bash
#SBATCH --job-name=prody_anm
#SBATCH --output=prody_%j.out
#SBATCH --error=prody_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# Adjust these lines to match your HPC system.
# Example options:
#   module load Python/3.12
#   module load Miniconda3
#
# If you created a venv:
#   source ~/venvs/prody-hpc/bin/activate
#
# If you created a conda env:
#   source "$(conda info --base)/etc/profile.d/conda.sh"
#   conda activate prody-hpc

source ~/venvs/prody-hpc/bin/activate

python compute_prody_frequencies.py \
  --structures-dir structures \
  --outdir results \
  --no-download \
  --n-modes all
