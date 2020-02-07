#!/bin/sh
#SBATCH --job-name=DDCP
#SBATCH --time=4:30:00
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=yuting.mou@uclouvain.be
#SBATCH --mail-type=ALL

# load modules required for execution
# can put in your bash_profile
# run script
export SLURM_NODEFILE=`generate_pbs_nodefile`
julia --depwarn=no --machinefile $SLURM_NODEFILE ./source/Main.jl
