#!/bin/sh
#SBATCH --job-name=PSPFullYear
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3900
#SBATCH --mail-user=yuting.mou@uclouvain.be
#SBATCH --mail-type=ALL

# load modules required for execution, can put in your bash_profile
# module load scip
# module load julia

# run script
julia --depwarn=no PSPMIP.jl
