#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=sigmat_nac
#SBATCH --output=/home/tl3031/project/git/sigmat/workflow/sigmat_nac_%j.out
#SBATCH --error=/home/tl3031/project/git/sigmat/workflow/sigmat_nac_%j.err
#SBATCH -p CSG
#SBATCH --mail-type=FAIL
#SBATCH --mail-user tl3031@cumc.columbia.edu

source ~/mamba_activate.sh

cd /home/tl3031/project/git/sigmat/workflow/
Rscript sigmat_nac.R
