#!/bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=Proposed
#SBATCH --output=Proposed.o%j
#SBATCH --time=0-40:59:59
# SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user xilu@ksu.edu
#SBATCH -J Proposed
module load R
R --no-save -q < TP_LADBLVC.R


